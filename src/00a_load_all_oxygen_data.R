# 
# Purpose: To load all oxygen data
# Author: Jake Diamond
# Date: 22 November 2022
# 

# Load libraries
library(plotly)
library(mgcv)
library(lubridate)
library(tidyverse)
library(tidytable)

# Load all hourly DO data---------------------------------------------------------------
# Data directory
data_dir <- file.path("data", "01_EDF", "raw", "dissolved_oxygen")

# files
files <- fs::dir_ls(data_dir, regexp = "\\.txt$")

# Load DO data
df <- files |> 
  map_dfr(read.table, skip = 2, header = T,
          col.names = c("datetime", "DO_mgL", "QC", "niv_mod",
                        "no_v_data", "validity", "niv_coh"),
          .id = "filename") |>
  mutate(datetime = dmy_hm(datetime))

# Get names of each site and position (upstream or downstream)
df <- df |>
  mutate(pos = if_else(str_detect(filename, "3_"), "down", "up"),
         site = case_when(str_detect(filename, "511") ~ "civaux",
                          str_detect(filename, "581") ~ "belleville",
                          str_detect(filename, "582") ~ "dampierre",
                          str_detect(filename, "583") ~ "stlaurent",
                          str_detect(filename, "584") ~ "chinon"))

# There's extra years of hourly data for 1992-1994 for dampierre up
# Already cleaned by Florentina for her thesis
df_dam_hr_1992 <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                      "F.M thesis-1990-1995.xlsx"),
                            sheet = 2) |>
  select(datetime = Date, DO_mgL = `Oxy`) |>
  mutate(site = "dampierre", pos = "up") |>
  filter(between(datetime, ymd(19920101), ymd_h(1994123123)))

# Same for Belleville 1990-1991, 1992 is missing completely
df_bel_hr_1990 <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                              "F.M thesis-1990-1995.xlsx"),
                                    sheet = 1) |>
  select(datetime = Date, DO_mgL = `Oxy`) |>
  mutate(site = "belleville", pos = "up") |>
  filter(between(datetime, ymd(19900101), ymd_h(1991123123))) |>
  mutate(datetime = round_date(datetime, "hour")) #for some reason the hours are off by one second??

# Remove that data from the raw data and add this in
df <- filter(df, !(site == "dampierre" & pos == "up" & 
                     datetime < ymd_h(1994123123))) |>
  filter(!(site == "belleville" & pos == "up" & 
                 datetime < ymd_h(1991123123))) |>
  bind_rows(df_dam_hr_1992) |>
  bind_rows(df_bel_hr_1990)

# Some simple data cleaning, DO can't be less than 0 or more than 30
df <- df |>
  mutate(DO_mgL = if_else(DO_mgL < 0 | DO_mgL > 30, NA_real_, DO_mgL))

# Summary of the data
df |>
  dplyr::group_by(site, pos) |>
  dplyr::summarize(nas = sum(is.na(DO_mgL)),
                   dat = sum(!is.na(DO_mgL)),
                   tot = sum(dat + nas),
                   per = nas / (dat + nas),
                   mindat = min(datetime),
                   maxdat = max(datetime))

# Plot all sites upstream and downstream raw data
p_DO_raw <- plot_ly(data = dplyr::filter(df, site != "stlaurent"), 
                    x = ~datetime,
                    y = ~DO_mgL,
                    color = ~site,
                    linetype = ~pos,
                    colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) |>
  add_trace(type = "scatter", mode='lines') |>
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)),
    title = "raw DO data")

# Load extra daily max min data -----------------------
# For Chinon, fill in from 1993–1998
df_chi <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                      "Chinon_1993_2005_FM.xlsx"),
                            sheet = 4) |>
  rename(date = DATE, min = `Mini Amont`, max = `Max amont`) |>
  mutate(site = "chinon", pos = "up", date = as.Date(date)) |>
  filter(date < ymd(19980101))

# For Dampierre, from 1990-1991
df_dam <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                      "F.M thesis-1990-1995.xlsx"),
                            sheet = 3) |>
  select(date = Date, min = `Oxy_min`, max = `Oxy_max`) |>
  mutate(site = "dampierre", pos = "up", date = as.Date(date)) |>
  filter(date < ymd(19920101))


# Get hour of daily max and min for hourly data by day of year--------
df_hr <- filter(df, site %in% c("dampierre", "chinon"), pos == "up") |>
  mutate(doy = yday(datetime),
         year = year(datetime)) |>
  group_by(site, pos, doy, year) |>
  filter(DO_mgL == max(DO_mgL, na.rm = T) | 
           DO_mgL == min(DO_mgL, na.rm = T)) |>
  mutate(hour = hour(datetime),
         minmax = if_else(DO_mgL == max(DO_mgL, na.rm = T), "max", "min"))

ggplot(data = drop_na(df_hr, site),
       aes(x = doy,
           y = hour,
           color = minmax)) +
  # facet_wrap(~site) +
  stat_summary() +
  stat_smooth()

# summary dataframe for fitting
df_hrsum <- ungroup(df_hr) |>
  group_by(doy, minmax) |>
  summarize(meanhr = mean(hour, na.rm = T))

# Get the gam fit
mod_max <- gam(meanhr ~ s(doy, bs = "cs"), data = filter(df_hrsum, minmax == "max"))
mod_min <- gam(meanhr ~ s(doy, bs = "cs"), data = filter(df_hrsum, minmax == "min"))

# Get the predicted max and min hours
predmaxmin <- tibble(doy = 1:365) |>
  mutate(maxhr = predict(mod_max, newdata = .),
         minhr = predict(mod_min, newdata = .))

# create hourly time series to fill in with max and min at estimated hour
# Chinon
df_tofill_chi <- tibble(datetime = seq(floor_date(min(df_chi$date), "hours"),
                                   ceiling_date(max(df_chi$date), "hours"),
                                   by = "hours")) |>
  mutate(doy = yday(datetime),
         date = date(datetime)) |>
  left_join(mutate(predmaxmin, across(everything(), round))) |>
  left_join(df_chi) |>
  mutate(DO_mgL = case_when(hour(datetime) == maxhr ~ max, 
                        hour(datetime) == minhr ~ min,
                        TRUE ~ NA_real_))

# Dampierre
df_tofill_dam <- tibble(datetime = seq(floor_date(min(df_dam$date), "hours"),
                                       ceiling_date(max(df_dam$date), "hours"),
                                       by = "hours")) |>
  mutate(doy = yday(datetime),
         date = date(datetime)) |>
  left_join(mutate(predmaxmin, across(everything(), round))) |>
  left_join(df_dam) |>
  mutate(DO_mgL = case_when(hour(datetime) == maxhr ~ max, 
                            hour(datetime) == minhr ~ min,
                            TRUE ~ NA_real_))

# Both
df_tofill <- bind_rows(df_tofill_dam, df_tofill_chi)

# Impute the max/min to hourly data --------------------------------------------------
# We know that late fall/winter time data there is rarely a diel signal (e.g. <1 mgL)
# So we can make imputation easier for these time by filling in the hours
# with a smooth average
df_tofill <- df_tofill |>
  group_by(date) |>
  mutate(DO_mgL = if_else(((doy < 60 | doy > 298) & abs(max-min) < 1),
                          mean(DO_mgL, na.rm = T),
                          DO_mgL)) |>
  ungroup()

# dataframe to fill, add in known data
df_tofill <- bind_rows(df_tofill, filter(df, site %in% c("chinon", "dampierre"), 
                                 pos == "up")) |>
  arrange(site, pos, datetime)

# Fill in missing NAs for dampierre using seasonal kalman filtering
# Get into time series format, then impute, only for first 8 years,
# less time to caclulate
df_imp <- df_tofill |>
  group_by(site, pos) |>
  nest() |>
  mutate(imp = map(data, imputeTS::na_kalman)) |>
  select(imp) |>
  unnest()

# Plot to review
z <- tibble(do = df_tofill_chi$DO_mgL,
            tst = filter(df_imp, site == "chinon", 
                         datetime <= ymd_h(1997123101)) |> select(DO_mgL),
            d = df_tofill_chi$datetime)

plot_ly(data = z, 
        x = ~d) |>
  add_trace(y = ~tst$DO_mgL, type = "scatter", mode='lines') |>
  add_trace(y = ~do, type = "scatter", mode='points', color = I("red"))

# Get all data together ---------------------------------------------------
df_all <- ungroup(df_imp) |>
  filter(!(site == "chinon" & datetime >= ymd_h(1998010101))) |>
  filter(!(site == "dampierre" & datetime >= ymd_h(1992010101))) |>
  bind_rows(df) |> 
  arrange(site, pos, datetime) |>
  select(site, pos, datetime, DO_mgL, QC)

# Save
saveRDS(df_all, file.path("data", "01_EDF", "do_timeseries.RDS"))
