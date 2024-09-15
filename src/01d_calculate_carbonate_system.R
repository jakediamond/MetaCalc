# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate carbonate system from Alkalinity and pH
# Date: 6 March 2023
# -------------------------------------

# Load libraries
library(tidyverse)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_carbonate_functions.R"))

# Load hourly pH data, recalculate carbonate system
# Need to do this for new data based on hourly alkalinity
# The first file is all pH data cleaned by An from 1990-2021
# The second file has all the 2022 data from EDF
df_pH <- readxl::read_xlsx(file.path("data", "05_hourly_data_clean", "2.Dampierre_hourly_newDepth.xlsx")) |>
  select(datetime, pH) |>
  bind_rows(read_csv(file.path("data", "01_EDF", "hourly EDF data", "hourly_pH.csv")) |>
              select(datetime, pH = DAM) |>
              filter(year(datetime) > 2021))
  
# Load hourly Alkalinity data, calculated by me from regression with SpC, discharge, and Ca
# Also, Ca, temperature, and specific conductance
df_AT <- readRDS(file.path("data", "03_CO2", "hourly_SpC_AT_Ca_system_v2.RDS")) |>
  select(datetime, AT,  SpC, Ca, temp)

# Join these and clean a bit
df_carb <- left_join(df_AT, df_pH) |>
  arrange(datetime) |>
  zoo::na.trim() |> # there are like 15 hours of NA pH at the end
  imputeTS::na_kalman() # and 20 hours of NA in between 2021 and 2022

# Now calculate the carbonate system
df_carb <- df_carb |>
  mutate(dic_sys = DIC(TK = temp + 273.15, AT = AT, pH, SpC)) |>
  unnest(cols = dic_sys) |>
  mutate(SI = SI(temp + 273.15, SpC, Ca, CO3_uM / 1000))

# Estimate the potential carbonate system at a given alkalinity and atm CO2
# This takes a while because it's an iterative solver
df_pot <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS")) %>%# this file has the atm. data
  janitor::clean_names() %>%
  select(datetime, pCO2_atm = atmosphere_co2_ppmv) %>%
  right_join(select(df_AT, datetime, AT, SpC, temp)) %>%
  # Last observation carried forward for filling in NAs at end of file for atm CO2
  imputeTS::na_locf() |>
  mutate(carbeq = seacarb::carb(flag = 24, 
                                var1 = pCO2_atm, 
                                var2 = AT / rhow(temp), # needs to be mol/kg
                                T = temp,
                                S = 1.6*10^-5 * SpC * 59,
                                pH = "F"))

# Get the error estimates in CO2
df_err <- df_carb |>
  drop_na() |>
  mutate(epH = if_else(year(datetime) < 2008, 0.1, 0.3), # from the manufacturer
         eAT = 133E-6, # from the regression model RMSE (mol/kg)
         ers = seacarb::errors(
           flag = 8,
           var1 = pH,
           var2 = AT / rhow(temp),
           evar1 = epH,
           evar2 = eAT,
           T = temp,
           S = 1.6*10^-5 * SpC * 59,
           eS = 0,
           pH = "F"))


# Get all that together
df <- df_pot |>
  unnest(cols = carbeq) |>
  mutate(DICeq_uM = DIC * rhow(temp) * 1000, 
         CO2eq_uM = CO2 * rhow(temp) * 1000) |>
  select(datetime, pCO2_atm, DICeq_uM, CO2eq_uM) |>
  right_join(df_carb) |>
  left_join(select(df_err, datetime, ers, temp) |> 
              unnest(ers) |>
              mutate(dCO2_uM = CO2 * rhow(temp) * 1000,
                     dDIC_uM = DIC * rhow(temp) * 1000) |>
              select(datetime, dCO2_uM, dDIC_uM)) |>
  rename(atmCO2_uatm = pCO2_atm)

# Save data
saveRDS(df, file.path("data", "hourly_carbonate_system_v2.RDS"))
