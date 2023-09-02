# 
# Purpose: Clean metabolism data
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(tidyverse)

# Load all the most recent metabolism 
df_met <- readRDS(file.path("data", "02_metabolism", "update_metabolism_results_07dec2022.RDS")) |>
  select(-`...1`) |>
  mutate(source = str_extract(source, "([^/]+$)"),
         source = str_remove(source, ".csv")) |>
  separate(col = source, into = c("type", "site", "pos", "period"), sep = "_") |>
  select(-type) |>
  filter(site == "dampierre", pos == "up")

# Count number of days with GPP and positive ER
sum(df_met$ER_mean > 0, na.rm = T) #456, = 3.9%
sum(df_met$ER_mean > 0 & df_met$ER_2.5pct < 0, na.rm = T) #319, 70.0%
sum(df_met$GPP_mean < 0, na.rm = T) # 1406 = 11.9%
sum(df_met$GPP_mean < 0 & df_met$GPP_97.5pct > 0, na.rm = T) #1251, 89.0%

# Step 1 Remove negative GPP and positive ER
# Set to 0 if biologically impossible, but 95% CI contains 0, otherwise NA
df_met_clean <- df_met |>
  mutate(GPP_mean = if_else(GPP_mean < 0 & df_met$GPP_97.5pct > 0, 0, GPP_mean),
         GPP_mean = if_else(GPP_mean < 0 & df_met$GPP_97.5pct < 0, NA_real_, GPP_mean),
         ER_mean = if_else(ER_mean > 0 & df_met$ER_2.5pct < 0, 0, ER_mean),
         ER_mean = if_else(ER_mean > 0 & df_met$ER_2.5pct > NA_real_, 0, ER_mean)) |>
  mutate(across(contains("GPP"), ~ifelse(is.na(GPP_mean), NA_real_, .)),
         across(contains("ER"), ~ifelse(is.na(ER_mean), NA_real_, .)))

# Step 2 is to remove really low K600 values and interpolate them
# Based on observations, anything less than 0.5 1/d is questionable
df_k <- select(df_met_clean, date, contains("K")) |>
  mutate(K600_daily_mean = if_else(K600_daily_mean < 0.5, NA_real_, K600_daily_mean)) |>
  mutate(across(contains("K"), ~ifelse(is.na(K600_daily_mean), NA_real_, .))) |>
  imputeTS::na_kalman()

sum(df_met_clean$K600_daily_mean < 0.5, na.rm = T) #594 or 5% of data

# Replace these new K values and recalculate flux
df_met_clean <- select(df_met_clean, -contains("K")) |>
  left_join(df_k)

saveRDS(df_met_clean, file.path("data", "dampierre_metabolism_clean.RDS"))
