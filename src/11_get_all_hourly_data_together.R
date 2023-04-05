# -------------------------------------
# Author: Jake Diamond
# Purpose: Get all hourly data together
# Date: 6 March 2023
# -------------------------------------

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
# library(seacarb)
library(tidyverse)
library(tidytable)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_co2_chem_enhance_functions.R"))

# Load data ---------------------------------------------------------------
# Load clean DO data
df_do <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS")) %>%
  filter(site == "dampierre", pos == "up")

# Quick clean
df_do <- df_do %>%
  mutate(datetime = floor_date(solar.time, "hours"))

# Load hourly CO2 data
df_co2 <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS")) %>%
  janitor::clean_names()

# Get data together
df_o2co2 <- select(df_do, O2 = DO.obs, O2_sat = DO.sat, 
                   solartime = solar.time, datetime, light) %>%
  right_join(df_co2) %>%
  mutate(O2_mmol_m3 = O2 * 1000 / 32, O2_sat_mmol_m3 = O2_sat * 1000 / 32) %>%
  mutate(O2ex = O2_mmol_m3 - O2_sat_mmol_m3, 
         O2per = O2_mmol_m3 / O2_sat_mmol_m3 * 100)

# Load NEP:CO2 archetypes
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

# Load hourly calcite precipitation
df_calc <- readRDS(file.path("data", "03_CO2", "hourly_calcite_precip.RDS"))

# Some quicks calcs and combine -------------------------------------------
# Calculate hourly NEP (mmol C/m2/hr) from change in DO
df_o2co2 <- df_o2co2 %>%
  mutate(Sc_O2 = 1801 - 120.1 * temperature + 3.782 * temperature^2 - 0.0476 * temperature^3,
         Sc_CO2 = 1742 - 91.24 * temperature + 2.208 * temperature^2 - 0.0219 * temperature^3,
         KO2 = k_co2_meta_1_d * (Sc_O2 / Sc_CO2)^-0.5,
         NEP = ((O2_mmol_m3 - lag(O2_mmol_m3)) - (KO2 / 24) * 
                  (O2_sat_mmol_m3 - O2_mmol_m3)) * depth)

# get all data together
df <- select(df_o2co2, datetime, solartime, discharge, depth, O2_mmol_m3,
             O2, O2ex, O2per, light, pH = p_h, temp = temperature, 
             CO2_atm_uatm = atmosphere_co2_ppmv, 
             CO2_atm_mmolm3 = atmosphere_co2_mmol_m3,
             CO2_uM = co2_mmol_m3, KCO2 = k_co2_meta_1_d, KO2 = KO2,
             CO2ex = co2ex, NEP_mmolO2m3 = NEP, CO2_flux_mmol_m2_d = co2_flux_mmol_m2_d,
             dCO2_flux = co2_flux_std_mmol_m2_d) %>%
  left_join(select(df_calc, datetime, SpC = filtered, Alk_molkg,
                   HCO3, CO3, DIC, Ca, LSI, CCPP, alkdif_uM, calc_uM = calc)) %>%
  mutate(date = date(datetime)) %>%
  dplyr::left_join(select(df_nepco2, year, month, date, archetype, 
                   chem_enh_daily = filtered_enh_mean,
                   NEP_daily = filtered_NEP_mean, CO2_daily = filtered_CO2_meanenh)) %>%
  mutate(
    # the chemical enhancement factor
    enh = chem_enh_fun(temp, pH, KCO2, depth)
  ) %>%
  mutate(CO2_flux_enh_mmol_m2_hr = CO2_flux_mmol_m2_d * enh / 24, #mmol/m2/hr
         dCO2_flux_enh_mmol_m2_hr = dCO2_flux * enh / 24,
         hr = hour(solartime),
         calc_mmol_m2_hr = calc_uM * depth) #mmol/m2/hr

# Fill in alkalinity for 1990-1991
df$Alk_molkg = ifelse(is.na(df$Alk_molkg), df_o2co2$alkalinity/1000, df$Alk_molkg)

# Calculate exDIC ---------------------------------------------------------
# Potential carbonate system
df_pot <-  seacarb::carb(flag = 24, # input eq pCO2 with atm. and ALK
                         df$CO2_atm_uatm,
                         df$Alk_molkg,
                         S = 0,
                         T = df$temp)

# Add this to the dataframe
df$DIC_eq_molkg <- df_pot$DIC
df$CO2_eq_molkg <- df_pot$CO2

# Calculate exDIC and exCO2
df <- mutate(df,
             exDIC_uM = DIC - (DIC_eq_molkg * 1E6),
             exCO2_uM = CO2_uM - (CO2_eq_molkg * 1E6))

saveRDS(df, file.path("data", "03_CO2", "all_hourly_data_complete.RDS"))
