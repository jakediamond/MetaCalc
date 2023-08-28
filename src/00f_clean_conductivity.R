# -------------------------------------
# Author: Jake Diamond
# Purpose: Clean conductivity data for Dampierre
# Date: 2023-05-01
# -------------------------------------
# Load libraries
library(plotly)
library(tidyverse)

# Load functions
source(file.path("src", "000_lowpass_function.R"))

# Load raw hourly data --------------------------------------------------------
# All time series from EDF, just want conductivity
df_raw <- readRDS(file.path("data", "01_EDF", "raw", 
                            "Raw_cond_oxy_pH_temp_1992_2022.rds")) |>
  select(site, position, datetime, SpC = Conductivity)

# Clean data for pH, temp, and DO
df_hr <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete.RDS")) |>
  select(datetime, discharge, O2, pH, temp)

# Get all data together
df <- filter(df_raw,
             site == "Dampierre") |>
  pivot_wider(names_from = position, values_from = SpC) |>
  left_join(df_hr) |>
  mutate(year = year(datetime), month = month(datetime), hr = hour(datetime))

# Only for Dampierre upstream
df_damup <- filter(df_raw,
                   site == "Dampierre",
                   position == "upstream") |>
  select(datetime, SpC)

# Predictability of conductivity based on other variables -----------------
# Quick regression just to see
mod <- lm(upstream ~ downstream*log(discharge)*month*temp,
           data = df)
summary(mod)

# plot(lm(upstream ~ downstream * log(discharge) * O2 * pH, 
#         data =slice_sample(df, prop = 0.01)))

# Cleaning of data --------------------------------------------------------
# Fill missing data with model
pred <- predict(mod, newdata = semi_join(df, df_damup, by = "datetime"))
df_damup$SpC <- if_else(is.na(df_damup$SpC), pred, df_damup$SpC)

# Low pass filter to get rid of noise
data <- lowpass_fun(df_damup, cutoff_frequency = 8) |>
    left_join(rename(df_damup, SpC_raw = SpC))

# Some values to control cleaning
prob <- 0.975
mult <- 1.5
minSpC <- 150
maxSpC <- 375

# Calculate upper 97.5% average delta DO (and DO) by month 
# Use this to set a bound on what to expect for big jumps
del_SpC <- data |>
  filter(between(filtered, minSpC, maxSpC)) |>
  mutate(month = month(datetime),
         year = year(datetime),
         dSpC = filtered - lag(filtered)) |>
  group_by(month) |>
  summarize(dSpC_lim = quantile(abs(dSpC), 
                                probs = prob, 
                                na.rm = TRUE))

# minimum and maximum monthly limits for conductivity
minmax_SpC <- data |>
  filter(between(filtered, minSpC, maxSpC)) |>
  mutate(month = month(datetime),
         year = year(datetime),
         date = date(datetime)) |>
  group_by(month, date) |>
  summarize(minSpC = min(filtered, na.rm = T),
            maxSpC = max(filtered, na.rm = T)) |>
  ungroup() |>
  group_by(month) |>
  summarize(SpC_limmax = quantile(maxSpC, 
                                  probs = prob, 
                                  na.rm = TRUE),
            SpC_limmin = quantile(minSpC, 
                                  probs = 1 - prob, 
                                  na.rm = TRUE))
  
# Make the data NA where there are big jumps, or just wrong data (anomalies)
# If the data jump more than 2x the 95% value for that time period
# Finally make DO NA if an anomaly is detected up to 2 hours ahead or behind
# this aids in smoothing the data later
data_natests <- data |> 
  mutate(month = month(datetime),
         year = year(datetime)) |>
  left_join(del_SpC) |>
  left_join(minmax_SpC) |>
  mutate(dSpC = filtered - lag(filtered),
         dSpCtest = if_else(abs(dSpC) > 2 * dSpC_lim | 
                              is.na(dSpC), "fail", "pass"),
         minSpCtest = if_else(filtered <= minSpC, "fail", "pass"),
         minSpCmonthtest = if_else(filtered <= SpC_limmin, "fail", "pass"),
         maxSpCmonthtest = if_else(filtered >= SpC_limmax, "fail", "pass")) |>
  mutate(passfail = paste(dSpCtest, minSpCtest, minSpCmonthtest, maxSpCmonthtest),
         SpC_test = if_else(str_detect(passfail, "fail"),  NA_real_, filtered )) |>
  mutate(SpC_test = if_else((is.na(lead(SpC_test)) & lead(dSpCtest) == "fail") |
                              (is.na(lag(SpC_test)) & lag(dSpCtest) == "fail"),
                            NA_real_,
                            SpC_test)) |>
  mutate(SpC_test = if_else((is.na(lead(SpC_test)) | is.na(lead(SpC_test, 2)) |
                        is.na(lag(SpC_test)) | is.na(lag(SpC_test, 2))),
                      NA_real_,
                      SpC_test))

# Check to see where original kalman filter was too filling too much (12 h)
data_nalengths <- data_natests |>
  dplyr::group_by(narun = {narun = rle(is.na(SpC_raw)); rep(seq_along(narun$lengths), 
                                                            narun$lengths)}) |>
  dplyr::mutate(namax = length(is.na(SpC_raw))) |>
  dplyr::mutate(SpC_real = if_else((is.na(SpC_raw) & namax > 12), NA_real_, SpC_test)) |>
  ungroup()

# Fill in NAs as best as possible with kalman filtering
# First scale the data so algorithm works
spc_scaled <- scales:::rescale(data_nalengths$SpC_real, c(0, 1))
ts_scaled <- ts(spc_scaled, frequency = 24)

# Do the kalman filter
spc_ts_clean <- imputeTS::na_kalman(ts_scaled, maxgap = 48)

# Rescale
ts_kalman <- scales:::rescale(as.numeric(spc_ts_clean),
                              c(min(data_nalengths$SpC_real, na.rm = T),
                                max(data_nalengths$SpC_real, na.rm = T)))

data_final <- data_nalengths
data_final$SpC_use <- ts_kalman


p_test <- plot_ly(data = data_final,
                  x  = ~datetime) |>
  add_trace(y = ~ SpC_use, type = "scatter", mode='lines') |> 
  add_trace(y = ~ SpC_raw, type = "scatter", mode='lines', color = I("gold")) |>
  # color = ~position, showlegend = FALSE) |>
  layout(xaxis = list(title = ""),
         yaxis = list(range = list(0, 450),
                      title = TeX("\\color{#1E88E5}{upstream~(mu~S~cm^{-1})}")),
         title = TeX("\\text{Specific Conductance}"),
         margin = list(r = 50)) |>
  config(mathjax = "cdn")

htmltools::browsable(p_test)

savefile <- select(data_final, datetime, SpC = SpC_use)
saveRDS(savefile, file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))
