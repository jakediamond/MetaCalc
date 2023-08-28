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
# Load NEP:CO2 archetypes
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

# Load clean hourly O2, SpC, AT, and Ca data
df <- readRDS(file.path("data", "03_CO2", "hourly_SpC_AT_Ca_O2system.RDS")) |>
  mutate(exO2  = O2_uM - O2sat_uM, 
         O2per = O2_uM / O2sat_uM * 100,
         year = year(date),
         month = month(date),
         hr = hour(datetime))

# Load hourly pH data, recalculate carbonate system
# Need to do this for new data based on hourly alkalinity
df_carb <- read_csv(file.path("data", "01_EDF", "hourly EDF data", "hourly_pH.csv")) |> # this could be a significant change as I am using measured pH now and not pH bias corrected according to Liu et al.
  select(datetime, pH = DAM) |>
  right_join(df) |>
  mutate(dic_sys = DIC(TK = temp + 273.15, AT = AT, pH, SpC)) |>
  unnest(cols = dic_sys) |>
  mutate(SI = SI(temp + 273.15, SpC, Ca, CO3_uM / 1000))

# Estimate the potential carbonate system at a given alkalinity and atm CO2
# This takes a while because it's an interative solver
df_pot <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS")) |> # this file has the atm. data
  janitor::clean_names() |>
  select(datetime, pCO2_atm = atmosphere_co2_ppmv) |>
  right_join(select(df_carb, datetime, AT, SpC, temp)) |>
  # Last observation carried forward for filling in NAs
  imputeTS::na_locf() |>
  mutate(carbeq = seacarb::carb(flag = 24, 
                                var1 = pCO2_atm, 
                                var2 = AT / rhow(temp), # needs to be mol/kg
                                T = temp,
                                S = 1.6*10^-5 * SpC * 59,
                                pH = "F",
                                k1k2 = "m06"))

# Get all that together
df_all <- df_pot |>
  unnest(cols = carbeq) |>
  mutate(DICeq_uM = DIC * 1E6, CO2eq_uM = CO2 * 1E6) |>
  select(datetime, pCO2_atm, DICeq_uM, CO2eq_uM) |>
  right_join(df_carb)
  

 # # Load hourly CO2 data
# df_co2 <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS"))
#   
# 
# # Get data together
# df_o2co2 <- select(df_do, ) %>%
#   right_join(df_co2) %>%


# Some quicks calcs and combine -------------------------------------------
# Calculate hourly FCO2 with enhancement
df_all <- df_all %>%
  mutate(Sc_CO2 = Sc(temp = temp),
         KCO2 = K600 * (600 / Sc_CO2)^0.5,
         enh = chem_enh_fun(temp, pH+0.13, KCO2, depth),
         FCO2 = KCO2 * (CO2_uM - CO2eq_uM) * depth,
         exCO2 = CO2_uM - CO2eq_uM,
         exDIC = DIC_uM - DICeq_uM)

# Estimate calcification from half change in alkalinity
# Also estimate based on dissolution/precip kinetics
df_all <- df_all %>%
  mutate(calc_mMm2hr = if_else(abs(SI) > 0.4,
                        -0.5 * (AT - lag(AT)) * depth * 1000,
                        0))
         # calc_molcm2hr = calc_rate(temp, SpC, SI, 
         #                           pH, CO2_uM / 1000))

# Add the trophlux states
df_all <- left_join(df_all,
                    select(df_nepco2, date, troph, sourcesink,
                           trophlux = archetype))

# Save this data
saveRDS(df_all, file.path("data", "hourly_data_final.RDS"))


# # get all data together
# df <- select(df_o2co2, datetime, solartime, discharge, depth, O2_mmol_m3,
#              O2, O2ex, O2per, light, pH = p_h, temp = temperature, 
#              CO2_atm_uatm = atmosphere_co2_ppmv, 
#              CO2_atm_mmolm3 = atmosphere_co2_mmol_m3,
#              CO2_uM = co2_mmol_m3, KCO2 = k_co2_meta_1_d, KO2 = KO2,
#              CO2ex = co2ex, NEP_mmolO2m3 = NEP, CO2_flux_mmol_m2_d = co2_flux_mmol_m2_d,
#              dCO2_flux = co2_flux_std_mmol_m2_d) %>%
#   left_join(select(df_calc, datetime, SpC = filtered, Alk_molkg,
#                    HCO3, CO3, DIC, Ca, LSI, CCPP, alkdif_uM, calc_uM = calc)) %>%
#   mutate(date = date(datetime)) %>%
#   dplyr::left_join(select(df_nepco2, year, month, date, archetype, 
#                    chem_enh_daily = filtered_enh_mean,
#                    NEP_daily = filtered_NEP_mean, CO2_daily = filtered_CO2_meanenh)) %>%
#   mutate(
#     # the chemical enhancement factor
#     enh = chem_enh_fun(temp, pH, KCO2, depth)
#   ) %>%
#   mutate(CO2_flux_enh_mmol_m2_hr = CO2_flux_mmol_m2_d * enh / 24, #mmol/m2/hr
#          dCO2_flux_enh_mmol_m2_hr = dCO2_flux * enh / 24,
#          hr = hour(solartime),
#          calc_mmol_m2_hr = calc_uM * depth) #mmol/m2/hr

# Fill in alkalinity for 1990-1991
# df$Alk_molkg = ifelse(is.na(df$Alk_molkg), df_o2co2$alkalinity/1000, df$Alk_molkg)

# # Calculate exDIC ---------------------------------------------------------
# # Potential carbonate system
# df_pot <-  seacarb::carb(flag = 24, # input eq pCO2 with atm. and ALK
#                          df$CO2_atm_uatm,
#                          df$Alk_molkg,
#                          S = 0,
#                          T = df$temp)
# 
# # Add this to the dataframe
# df$DIC_eq_molkg <- df_pot$DIC
# df$CO2_eq_molkg <- df_pot$CO2
# 
# # Calculate exDIC and exCO2
# df <- mutate(df,
#              exDIC_uM = DIC - (DIC_eq_molkg * 1E6),
#              exCO2_uM = CO2_uM - (CO2_eq_molkg * 1E6))

# saveRDS(df, file.path("data", "03_CO2", "all_hourly_data_complete.RDS"))
