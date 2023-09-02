# 
# Purpose: Propagate error from metabolism to NEP, CO2 flux
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(tidyverse)

# Load error propagation function
source(file.path("src", "000_error_propagation_function.R"))
source(file.path("src", "000_carbonate_functions.R"))

# Load data and clean---------------------------------------------------------------
# Load all the metabolism data for all sites, clean for site info, then get Dampierre
df_met <- readRDS(file.path("data", "dampierre_metabolism_clean.RDS"))

# Load discharge, and K600 from Raymond equations data
df_kray <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                       "DAM_K600_Flux_all_Eq.xlsx"))

# Load carbonate system data with uncertainty for CO2
df_carb <- readRDS(file.path("data", "hourly_carbonate_system.RDS"))

# Calculate Schmidt number CO2 (2 eqns in Raymond et al. 2012), and 1 in Wanninkoff 2014
df_sc <- df_carb |>
  select(datetime, temp) |>
  mutate(Sc_CO2_1 = Sc(temp = temp), #this is the Wanninkoff, I made a function for this
         Sc_CO2_2 = 1742 - 91.24 * temp + 2.208 * temp^2 - 0.0219 * temp^3,
         Sc_CO2_3 = 1911 - 118.11 * temp + 3.453 *temp^2 - 0.0413 * temp^3) |>
  rowwise() |> 
  # Mean and uncertainty
  mutate(Sc_CO2_mean = mean(c(Sc_CO2_1, Sc_CO2_2, Sc_CO2_3)),
         dSc_CO2_mean = sd(c(Sc_CO2_1, Sc_CO2_2, Sc_CO2_3))) |>
  ungroup() |>
  mutate(date = date(datetime))

# Count number of days with GPP and positive ER
sum(df_met$ER_mean > 0, na.rm = T) #456, = 3.9%
sum(df_met$ER_mean > 0 & df_met$ER_2.5pct < 0, na.rm = T) #319, 70.0%
sum(df_met$GPP_mean < 0, na.rm = T) # 1406 = 11.9%
sum(df_met$GPP_mean < 0 & df_met$GPP_97.5pct > 0, na.rm = T) #1251, 89.0%

# Correlation of GPP and ER to account for in error propagation
cor_ge <- cor(df_met$GPP_mean, y = df_met$ER_mean, use = "pairwise.complete.obs")

# Data cleaning -----------------------------------------------------------

  
# Need to use a sampling of K600 for its standard deviation
# Otherwise we can have -K600, Which is impossible
df_met_clean <- df_met_clean |>
  select(date, 
         GPP_mean, dGPP_mean = GPP_daily_sd, 
         GPP_2.5 = GPP_2.5pct, GPP_97.5 = GPP_97.5pct,
         ER_mean, dER_mean = ER_daily_sd, 
         ER_2.5 = ER_2.5pct, ER_97.5 = ER_97.5pct,
         K600_mean = K600_daily_mean, dK600_mean = K600_daily_sd,
         K600_2.5 = K600_daily_2.5pct, K600_97.5 = K600_daily_97.5pct)

# Propagate errors ----------------------------------------------
# Propagate error from GPP and ER to NEP
# Assume symmetric error around the mean (mostly true) and account for covariance
# between GPP and ER
df_met_err <- df_met_clean |>
  mutate(NEP_mean = GPP_mean + ER_mean,
         dNEP_mean = sqrt(dGPP_mean^2 + dER_mean^2 + 2 * cor_ge * dGPP_mean * dER_mean),
         NEP_2.5 = NEP_mean - 1.96*dNEP_mean, # estimate 95% credible interval for NEP
         NEP_97.5 = NEP_mean + 1.96*dNEP_mean)

# Calculate KCO2 from metabolism K600, propagate error in Sc and K600
df_KCO2 <- select(df_met_clean, date, contains("K600_mean"), K600_2.5, K600_97.5) |>
  right_join(select(df_sc, date, datetime, Sc_CO2_mean, dSc_CO2_mean), by = "date") |>
  mutate_with_error(KCO2_met_mean ~ K600_mean * (600 / Sc_CO2_mean)^(0.5)) |>
  rename_with(~ str_replace(.x, 
                            pattern = "K600", 
                            replacement = "K600_met"), 
              matches("K600")) 

# Now we want to calculate CO2 fluxes
df_CO2 <- df_carb |>
  df_KCO2 |>
  left_join(select(df_kray, date, depth = `Depth (m)`)) |>
  left_join(df_co2sys) |>
  tidytable::mutate(CO2_flux = depth * (CO2_mmolm3 - CO2_atm) * KCO2_met_mean,
                    co2_dist = tidytable::map2(CO2_mmolm3, dCO2_mmolm3, ~rnorm(10000, .x, .y)),
                    k_dist = tidytable::map2(K600_met_2.5, K600_met_97.5, 
                                             ~sample(c(runif(100, .x, .y),
                                                       runif(100, .x, .y)),
                                                     10000, replace = TRUE)),
                    f_dist = tidytable::map2(co2_dist, k_dist, ~depth * (.x - CO2_atm) * .y),
                    dCO2_flux = tidytable::map_dbl(f_dist, sd, na.rm = T) / sqrt(10000),
                    CO2_flux_2.5 = tidytable::map_dbl(f_dist, quantile, 0.025, na.rm = T),
                    CO2_flux_97.5 = tidytable::map_dbl(f_dist, quantile, 0.975, na.rm = T)) |>
  select(-KCO2_ray_mean, -dKCO2_ray_mean, -co2_dist, -k_dist, -f_dist)
# Save all data for future reference --------------------------------------
df_save <- df_CO2 |>
  select(-contains("ray")) |>
  left_join(select(df_met_err, -contains("K600"))) |>
  mutate(across(contains(c("GPP", "ER", "NEP")), ~.*1000/32)) #from g O2 to mmmol
dt <- as.Date(Sys.time())
saveRDS(df_save, file = file.path("data", "03_CO2",
                                  paste0("CO2_with_uncertainty_dampierre_up_",
                                         dt,
                                         ".RDS")))

write_excel_csv(df_save, file = file.path("data", "03_CO2",
                                          paste0("CO2_with_uncertainty_dampierre_up_",
                                                 dt,
                                                 ".csv")))

