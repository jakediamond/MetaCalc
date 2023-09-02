# -------------------------------------
# Author: Jake Diamond
# Purpose: Fit regressions for Alkalinity and Ca to create timeseries
# Date: 2023-04-28
# -------------------------------------

library(tidyverse)

# Clean up input data for modeling ----------------------------------------
# Load all the clean K600 data
df_k <- df_met <- readRDS(file.path("data", "dampierre_metabolism_clean.RDS")) |>
  select(date, K600 = K600_daily_mean)

# Hourly O2 data
df_o2 <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS")) |>
  filter(site == "dampierre", pos == "up") |>
  mutate(datetime = floor_date(solar.time, "hours")) |>
  select(datetime, O2_mgL = DO.obs, O2sat_mgL = DO.sat, depth, temp = temp.water) |>
  mutate(O2_uM = O2_mgL * 1000 / 32, O2sat_uM = O2sat_mgL * 1000 / 32,
         date = date(datetime))

# Calculate hourly NEP
df_o2 <- df_o2 |>
  left_join(df_k) |>
  mutate(Sc_O2 = 1801 - 120.1 * temp + 3.782 * temp^2 - 0.0476 * temp^3,
         KO2 = K600 * (600 / Sc_O2)^0.5) |>
  arrange(datetime) |>
  # smooth the K and depth data, 3 day moving average
  mutate(K600_smooth = slider::slide_dbl(K600, mean, .before = 3*24, na.rm = T),
         KO2_smooth = slider::slide_dbl(KO2, mean, .before = 3*24, na.rm = T),
         depth_smooth = slider::slide_dbl(depth, mean, .before = 3*24, na.rm = T),
         NEP = ((O2_uM - lag(O2_uM)) - (KO2_smooth / 24) * 
                  (O2sat_uM - O2_uM)) * depth_smooth)

# Save this data
saveRDS(df_o2, file.path("data", "hourly_o2_system.RDS"))
