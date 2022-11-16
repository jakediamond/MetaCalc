# -------------------------------------
# Author: Jake Diamond
# Purpose: Final cleaning of data with random forest filling of NAs
# Date: 27 October 2022
# -------------------------------------

# Load libraries
library(lubridate)
library(tidyverse)
library(missForest)
library(plotly)
library(tidytable)

# Set seed
set.seed(42)

# Load physical data
df_phys <- readRDS(file.path("data", "05_hourly_data_clean", "temp_discharge_rad_data.RDS"))

# Load pH data
df_pH <- read_csv(file.path("data", "01_EDF", "hourly EDF data", "hourly_pH.csv")) %>%
  rename(dampierre = DAM, belleville = BEL, chinon = CHB) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "pH")

# Load DO data
df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))

# Combine data and clear workspace
df <- left_join(df_DO, df_phys) %>%
  ungroup() %>%
  left_join(df_pH)

rm(df_DO, df_phys, df_pH)

# Remove entire years of missing data (more than a half year)
# Remove DO data again when NA is longer than 3 days, initial na_kalman didn't work
df <- df %>%
  group_by(site, pos, year) %>%
  filter.(!(sum(is.na(DO_mgL)) > 180 * 24)) %>%
  ungroup() %>%
  group_by(site, year, pos, narun) %>%
  mutate(DO_use = if_else((namax > 24 * 3)  & is.na(DO), NA_real_, DO_use)) %>%
  ungroup()

# A bit more cleaning -----------------------------------------------------
# Another cleaning function
do_clean_fun <- function(data) {
  # Fill in NAs as best as possible with kalman filtering
  do_ts <- ts(data$DO_use, deltat = 1/(365*24))
  do_ts_clean <- imputeTS::na_seasplit(do_ts, algorithm = "kalman")
  data$DO_use <- as.numeric(do_ts_clean)
  data
}
# Apply functions ---------------------------------------------------------
# Clean the data and use the lowpass filter
df_n <- df %>%
  group_by(site, pos) %>%
  nest() %>%
  mutate(clean = map(data, possibly(do_clean_fun, NA_real_)))

# Get all in one dataframe
df_DO <- filter.(df_n, !(is.na(clean))) %>%
  unnest(clean)

a %>%
  group_by(year) %>%
  summarize(x = sum(is.na(DO_use)) / 24)

do_ts <- ts(a$DO_use, deltat = 1/(365*24))
do_ts_clean <- imputeTS::na_seasplit(do_ts, algorithm = "kalman")
data$DO_use <- as.numeric(do_ts_clean)
data
# Clean up for some known periods of bad data
# Belleville up from 1997 july 30 00:00 to 1997 august 21 13:00
clean1 <- filter(df, site == "belleville", pos == "down",
                 between(datetime, ymd_h(1997073000), ymd_h(1997082113))) %>%
  distinct(site, pos, date, datetime, DO_use) %>%
  pivot_wider(names_from = pos, values_from = DO_use)


# Belleville down from 2001 aug 22 04:00 to 2001 august 27 09:00, make NA and
# fill in with naseades
clean2 <- df %>%
  mutate(DO_use = if_else(site == "belleville", pos == "down",
                 between(datetime, ymd_h(2001082204), ymd_h(2001082709)),
                 NA_real_, DO_use))



# Data for metabolism -----------------------------------------------------
# Calculate DO saturation
df_DOsat <- select(df_DO, site, datetime = soltime, pos, wtr = temp_C) %>% 
  arrange(site, pos, datetime) %>%
  ungroup() %>%
  select(datetime, wtr) %>%
  as.data.frame() %>%
  LakeMetabolizer::o2.at.sat() %>%
  select(-datetime)

# save data for metabolism
savefile <- select(df_DO, site, date, solar.time = soltime,
                   period, pos, DO.obs = DO_use, 
                   temp.water = temp_C, rad_Wm2, depth = depth_m,
                   discharge = Q_m3s) %>%
  ungroup() %>%
  mutate(light = streamMetabolizer::convert_SW_to_PAR(rad_Wm2),
         DO.sat = df_DOsat$do.sat) %>%
  select(-rad_Wm2, -discharge, -date)

# All data
saveRDS(savefile, file.path("data", "02_metabolism", "hourly_inputs.RDS"))

# All discharge data
qfile <- select(savefile, site, pos, period, discharge, date) %>%
  group_by(site, pos, period, date) %>%
  summarize(discharge.daily = mean(discharge)) %>%
  ungroup()

# All discharge data
saveRDS(qfile, file.path("data", "02_metabolism", "daily_Q_inputs.RDS"))

# Split data
select(savefile, -discharge, -date) %>%
  mutate(filename = paste(site, pos, period, sep = "_")) %>%
  group_split(site, pos, period) -> list_of_dfs

# name of each datafile
list_of_dfs %>%
  map(~pull(., filename)) %>%
  map(~unique(.)) -> names(list_of_dfs)

imap(list_of_dfs, ~saveRDS(.x, file = file.path("data", "02_metabolism", 
                                                paste0(.y, "_for_meta", ".RDS"))))

# daily discharge, too
select(savefile, site, pos, period, discharge, date) %>%
  group_by(site, pos, period, date) %>%
  summarize(discharge.daily = mean(discharge)) %>%
  ungroup() %>%
  mutate(filename = paste(site, pos, period, sep = "_")) %>%
  group_split(site, pos, period) -> list_of_dfs_q

# name of each datafile
list_of_dfs_q %>%
  map(~pull(., filename)) %>%
  map(~unique(.)) -> names(list_of_dfs_q)

imap(list_of_dfs_q, ~saveRDS(.x, file = file.path("data", "02_metabolism", 
                                                paste0(.y, "_daily_Q",".RDS"))))
# More cleaning -----------------------------------------------------------


# Belleville up from 1997 july 30 00:00 to 1997 august 21 13:00
clean1 <- filter(df, site == "belleville", pos == "down",
                 between(datetime, ymd_h(1997073000), ymd_h(1997082113))) %>%
  distinct(site, pos, date, datetime, DO_use) %>%
  pivot_wider(names_from = pos, values_from = DO_use)

df_clean1 <- df %>%
  left_join(clean1) %>%
  mutate(DO_use = case_when(site == "belleville" & pos == "up" &
                              between(datetime, ymd_h(1997073000), ymd_h(1997082113)) ~
          down,
         TRUE ~ DO_use))

# Subsample just to test
x <- filter(df, site == "chinon", 
            between(date, ymd(20050101), ymd(20050301)), 
            pos == "down") %>%
  select(datetime, site, DO_use) %>%
  # pivot_wider(names_from = pos, values_from = DO) %>%
  left_join(filter(df, pos == "up") %>%
              select(datetime, site, temp_C, rad_Wm2, Q_m3s, pH, up = DO_use) %>%
              distinct(),
            by = c("datetime", "site")) %>%
  mutate(hour = hour(datetime),
         month = month(datetime))

summary(x)
library(nlme)
mod <- gls(DO_use ~ temp_C*rad_Wm2*Q_m3s*pH*up, data = drop_na(x), corr = corAR1())
summary(mod)
plot(mod)
acf(residuals(mod, type = "normalized"))
plot(residuals(mod, type = "normalized"))

modpred <- predict(mod, newdata = newd)
# plot(mod)
# y = predict(mod, newdata = x)
# plot(modpred)
# lines(x$DO_use)
plot(xpred$DO_use)
lines(modpred)
lines(newd$up, col = "red")
mean(sqrt((y-x$down)^2))

xpred <- filter(df, site == "chinon", between(date, ymd(20030129), ymd(20030228)), 
                pos == "down")
newd <- filter(df, site == "chinon", between(date, ymd(20030129), ymd(20030228)), 
               pos == "down") %>%
  select(datetime, site, DO_use) %>%
  # pivot_wider(names_from = pos, values_from = DO) %>%
  left_join(filter(df, pos == "up") %>%
              select(datetime, site, temp_C, rad_Wm2, Q_m3s, pH, up = DO_use) %>%
              distinct(),
            by = c("datetime", "site")) %>%
  mutate(hour = hour(datetime),
         month = month(datetime))
summary(newd)


a = filter(df, site == "dampierre", year == 1999)

xdf <- x %>%
  select(-site, -datetime) %>%
  as.data.frame()
summary(x)
x_im <- missForest(xdf, verbose = TRUE, ntree = 50)

x_im$OOBerror

y <- x_im$ximp




plot_ly(data = bind_cols(y, datetime = x$datetime) %>%
          pivot_longer(cols = c(up, down)), 
        x = ~datetime,
        y = ~value,
        color = ~name) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))



plot_ly(data = filter(df, site == "chinon", pos == "up", between(year, 2005, 2007)) %>%
          select(datetime, DO_mgL, DO_use) %>%
          left_join(select( bind_cols(y, datetime = x$datetime), 
                            datetime, DO)) %>%
          pivot_longer(cols = c(DO, DO_use, DO_mgL)), 
        x = ~datetime,
        y = ~value,
        color = ~name) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))


plot_ly(data = filter(df_clean1, site == "belleville", 
                      between(date, ymd(20010601), ymd(20010830))),
        x = ~datetime,
        y = ~DO_mgL,
        color = ~pos) %>%
  add_trace(type = "scatter", mode='lines')

plot_ly(data = filter.(df_DO, site == "dampierre"),
        x = ~datetime,
        y = ~DO_use,
        color = ~pos) %>%
  add_trace(type = "scatter", mode='lines')


# Save as xlsx with sheets
# Split one data frame per site
df_w %>%
  group_split(site) -> list_of_dfs
# name of each sheet will be the site
list_of_dfs %>%
  map(~pull(., site)) %>%
  map(~unique(.)) -> names(list_of_dfs)

list_of_dfs %>%
  writexl::write_xlsx(path = file.path("data", "05_hourly_clean_data", "all_data.xlsx")