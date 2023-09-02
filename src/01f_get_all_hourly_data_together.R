# -------------------------------------
# Author: Jake Diamond
# Purpose: Get all hourly data together
# Date: 6 March 2023
# -------------------------------------

# Load libraries
library(tidyverse)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_co2_chem_enhance_functions.R"))
source(file.path("src", "000_carbonate_functions.R"))
source(file.path("src", "000_calcite_kinetics_function.R"))
# Load data ---------------------------------------------------------------
# Load hourly fco2 and carbonate system with daily metabolism + uncertainty
df_nepco2 <- readRDS(file.path("data", "nep_fco2_uncertainty.RDS"))

# Load hourly O2 data
df_o2 <- readRDS(file.path("data", "hourly_o2_system.RDS")) |>
  mutate(exO2  = O2_uM - O2sat_uM, 
         O2per = O2_uM / O2sat_uM * 100)

# Get daily discharge and radiation data
df_q <- readRDS(file.path("data", "05_hourly_data_clean", 
                          "temp_discharge_rad_data.RDS")) |>
  filter(site == "dampierre") |>
  select(datetime, Q_m3s, rad_Wm2)

# Some quicks calcs and combine -------------------------------------------
df_all <- df_nepco2 |>
  left_join(df_o2) |>
  left_join(df_q) |>
  mutate(year = year(datetime),
         month = month(datetime),
         date = date(datetime),
         hr = hour(datetime)) |>
  relocate(year, month, date, datetime, hr, Q_m3s, depth_m = depth, 
           temp, pH, AT_mM = AT, Ca_mM = Ca)

# Estimate calcification from half change in alkalinity (Escoffier et al 2023)
df_all <- df_all |>
  mutate(calc_mMm2hr = if_else(abs(SI) > 0.4,
                        -0.5 * (AT_mM - lag(AT_mM)) * depth_m * 1000,
                        0)) |>
  mutate(
    # calcification shouldn't be positive if SI is < 0
    calc_mMm2hr = if_else(SI < 0 & calc_mMm2hr > 0, 0, calc_mMm2hr),
    # And the inverse
    calc_mMm2hr = if_else(SI > 0 & calc_mMm2hr < 0, 0, calc_mMm2hr)
         ) |>
  mutate(
    # Add the chemical enhancement
    FCO2_enh = FCO2 * enh,
    dFCO2_enh = dFCO2 * enh
  )

# Save this data
saveRDS(df_all, file.path("data", "hourly_data.RDS"))
