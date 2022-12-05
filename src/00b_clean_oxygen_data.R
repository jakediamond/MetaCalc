# 
# Purpose: To clean DO time-series for middle Loire and Vienne, reduce outliers
# Author: Jake Diamond
# Date: 25 October 2022
# 

# Load libraries
library(zoo)
library(signal)
library(imputeTS)
library(plotly)
library(lubridate)
library(tidyverse)
library(tidytable)

# Load data ----------------------------------------
# first round of cleaned DO data
df <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))

# Load temp, rad, discharge data
df_rtq <- readRDS(file.path("data", "05_hourly_data_clean", 
                            "temp_discharge_rad_data.RDS")) %>%
  arrange(site, datetime)

# Calculate DO saturation to help account for temperature when deciding
# if data is an outlier or not
df_rtq <- select(df_rtq, site, datetime, wtr = temp_C) %>% 
  arrange(site, datetime) %>%
  select(datetime, wtr) %>%
  as.data.frame() %>%
  LakeMetabolizer::o2.at.sat() %>%
  select(-datetime) %>%
  bind_cols(df_rtq)

# get all data together
df <- left_join(df, df_rtq) %>%
  mutate(DOper = DO_use / do.sat * 100)

# Plot to inspect ---------------------------------------------------------
# Plot to inspect for manual cleaning
p <- plot_ly(data = dplyr::filter(df, site == "belleville"), 
                      x = ~datetime) %>%
  add_trace(y = ~DOper,
            color = ~pos,
            colors = c("#1E88E5", "#FFC107"), 
            type = "scatter", mode='lines') %>%
  add_trace(y = ~Q_m3s, type = "scatter", mode='lines',
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                title = "discharge (m3/s)"),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (% sat)",
                 range = list(0, 270)),
    title = "1st round clean DO data")

htmltools::browsable(p)

# Remove data to clean
rm(df_rtq)
# Add the flags from manual cleaning -------------------------------------
# flags
flags <- readxl::read_xlsx(file.path("data", "do_data_flags.xlsx"))

# Combine flags with data
# function to help speed up
fuz_fun <- function(data, flag) {
  data %>% fuzzyjoin::fuzzy_left_join(flag,
                             by = c("datetime" = "start",
                                    "datetime" = "end"),
                             match_fun = list(`>=`, `<=`))
}

# Much easier/faster to do this with a few columns and split by sites
df_flags <- df %>%
  dplyr::select(site, pos, datetime) %>% 
  dplyr::group_by(site, pos) %>%
  tidyr::nest() %>%
  dplyr::left_join(dplyr::group_by(flags, site, pos) %>% 
                     tidyr::nest(flag = c(start, end, flag))) %>%
  dplyr::summarize(map2_df(data, flag, fuz_fun))

# Do the cleaning based on manual flags
# Linear stuff first
df_clean_linear <- df_flags %>%
  dplyr::filter(flag == "linear") %>% #do linear first then rejoin other data
  left_join(df, by = c("site", "pos", "datetime")) %>%
  dplyr::group_by(site, pos, start, end) %>%
  tidyr::nest() %>%
  dplyr::mutate(DOper_clean = map(data, ~ .$DOper - 
                                    coef(lm(DOper~seq_along(datetime), .))[2] * #everything based on %sat
                                    seq_along(.$datetime))) %>%
  unnest(c(data, DOper_clean)) %>%
  ungroup()

# Linear with drop
df_clean_linear_drop <- df_flags %>%
  left_join(df, by = c("site", "pos", "datetime")) %>%
  dplyr::filter(str_detect(flag, "linear_drop")) %>% #do linear first then rejoin other data
  dplyr::mutate(drop = as.numeric(word(flag, 3, sep = "_"))) %>%
  dplyr::group_by(site, pos, start, end, drop) %>%
  tidyr::nest() %>%
  dplyr::mutate(DOper_clean = map2(data, drop, ~ .y + .x$DOper - 
                                    coef(lm(DOper~seq_along(datetime), .x))[2] * #everything based on %sat
                                    seq_along(.x$datetime))) %>%
  unnest(c(data, DOper_clean)) %>%
  ungroup()
  
# Everything else
df_clean <- bind_rows(df_clean_linear, df_clean_linear_drop) %>%
  bind_rows(left_join(df_flags, df, by = c("site", "pos", "datetime")) %>% 
              filter(!str_detect(flag, "linear") | is.na(flag))) %>%
  arrange(site, pos, datetime) %>%
  mutate(DOper_clean = case_when(flag == "remove" ~ NA_real_,
                                 str_detect(flag, "rm") ~ if_else(hour(datetime) == as.numeric(word(flag, 2, sep = "_")),
                                                                  NA_real_,
                                                                  DOper),
                                 str_detect(flag, "rd") ~ {if_else(hour(datetime) >= as.numeric(word(flag, 2, sep = "_")) |
                                                                     hour(datetime) <= as.numeric(word(flag, 3, sep = "_")),
                                                                   NA_real_,
                                                                   DOper + as.numeric(word(flag, 4, sep = "_")))},
                                 !is.na(as.numeric(flag)) ~ DOper + as.numeric(flag),
                                 str_detect(flag, "remove_") ~ {if_else(hour(datetime) >= as.numeric(word(flag, 2, sep = "_")) |
                                                                          hour(datetime) <= as.numeric(word(flag, 3, sep = "_")),
                                                                        NA_real_,
                                                                        DOper)},
                                 str_detect(flag,"linear") ~ DOper_clean,
                                 TRUE ~ DOper)) %>%
  mutate(DO_use_clean = DOper_clean / 100 * do.sat) #convert back to DO from %sat

# Plot
p2 <- plot_ly(data = dplyr::filter(df_clean, site == "belleville"),
             x = ~datetime) %>%
  add_trace(y = ~DOper_clean,
            color = ~pos,
            colors = c("#1E88E5", "#FFC107"),
            type = "scatter", mode='lines') %>%
  add_trace(y = ~Q_m3s, type = "scatter", mode='lines',
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = "discharge (m3/s)"),
         xaxis = list(title = ""),
         yaxis = list(title = "DO (mg/L)",
                      range = list(0, 270)),
         title = "1st round clean DO data")

htmltools::browsable(p2)

saveRDS(df_clean, file.path("data", "05_hourly_data_clean", 
                            "DO_cleaned_part2.RDS"))
rm(df, df_flags, flags)

# Function to add progress bar --------------------------------------------
map_progress <- function(.x, .f, ..., .id = NULL) {
  .f <- purrr::as_mapper(.f, ...)
  pb <- progress::progress_bar$new(total = length(.x), 
                                   format = " [:bar] :current/:total (:percent) eta: :eta",
                                   force = TRUE)
  
  f <- function(...) {
    pb$tick()
    .f(...)
  }
  purrr::map(.x, f, ..., .id = .id)
}
# Imputation --------------------------------------------------------------
df_clean <- readRDS(file.path("data", "05_hourly_data_clean", 
                            "DO_cleaned_part2.RDS"))

# Finally fill in with na_seadec, this will take a very long time
df_imp <- df_clean %>%
  arrange(site, pos, datetime) %>%
  nest_by(site, pos, period) %>%
  mutate(ts = map(data, ~ts(.$DOper_clean, frequency = 24)),
         imputed = map_progress(ts, possibly(~na_seasplit(.x, "kalman", 
                                               # find_frequency = TRUE,
                                               maxgap = 24*40),
                                     otherwise = NULL)))

saveRDS(df_imp, file.path("data", "05_hourly_data_clean",
                          "DO_cleaned_part3_split.RDS"))

x= select(df_imp, -ts) %>%
  # filter(!map_lgl(imputed, is.null)) %>%
  unnest(c(data, imputed))

p3 <- plot_ly(data = dplyr::filter(ungroup(x), site == "belleville"),
             x = ~datetime) %>%
  add_trace(y = ~imputed,
            color = ~pos,
            colors = c("#1E88E5", "#FFC107"),
            type = "scatter", mode='lines') %>%
  add_trace(y = ~Q_m3s, type = "scatter", mode='lines',
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = "discharge (m3/s)"),
         xaxis = list(title = ""),
         yaxis = list(title = "DO (mg/L)",
                      range = list(0, 270)),
         title = "1st round clean DO data")

htmltools::browsable(p3)
# One last clean ----------------------------------------------------------
source(file.path("src", "000_DO_cleaning_functions.R"))

df_imp <- readRDS(file.path("data", "05_hourly_data_clean",
                          "DO_cleaned_part3_split.RDS")) %>%
  select(-ts) %>%
  unnest(c(data, imputed))

# Clean the data and use the lowpass filter
df_final <- df_imp %>%
  select(site, pos, period, datetime, imputed, do.sat) %>%
  mutate(DO_mgL = imputed / 100 * do.sat) %>%
  select(-imputed, -do.sat) %>%
  nest_by(site, pos, period) %>%
  mutate(clean = map(data, lowpass_fun, variable = DO_mgL, cutoff_frequency = 0.15))

df_final <- unnest(df_final, clean)

p3 <- plot_ly(data = dplyr::filter(ungroup(df_final2), site == "chinon"),
              x = ~datetime) %>%
  add_trace(y = ~DO,
            color = ~pos,
            colors = c("#1E88E5", "#FFC107"),
            type = "scatter", mode='lines') %>%
  # add_trace(y = ~Q_m3s, type = "scatter", mode='lines',
  #           yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = "discharge (m3/s)"),
         xaxis = list(title = ""),
         yaxis = list(title = "DO (mg/L)",
                      range = list(0, 25)),
         title = "1st round clean DO data")

htmltools::browsable(p3)

# fill back for NAs
df_final2 <- df_final %>%
  dplyr::group_by(site, pos, period,
                  narun = {narun = rle(is.na(DO_mgL)); rep(seq_along(narun$lengths), 
                                                           narun$lengths)}) %>%
  dplyr::mutate(namax = length(is.na(DO_mgL))) %>%
  dplyr::mutate(DO = if_else((is.na(DO_mgL) & namax > 40*24), NA_real_, filtered)) %>%
  ungroup() %>%
  select(site, pos, period, datetime, DO)


saveRDS(df_final2, file.path("data", "05_hourly_data_clean",
                          "DO_cleaned_final.RDS"))
