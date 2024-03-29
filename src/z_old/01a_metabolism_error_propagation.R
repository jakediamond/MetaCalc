# 
# Purpose: Propagate error from metabolism to NEP, CO2 flux
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(tidyverse)

# Load error propagation function
source(file.path("src", "000_error_propagation_function.R"))

# Load data and clean---------------------------------------------------------------
df_met <- readRDS(file.path("data", "02_metabolism", "update_metabolism_results_07dec2022.RDS"))
# # Load metabolism data
# list.files(path = file.path("data", "02_metabolism", "newest_metabolism"),
#            pattern = "GetFitDaily",
#            full.names = TRUE) |>
#   map_dfr(read_csv, show_col_types = FALSE, id = "filename")
# 
# # Clean just a bit for site and position
# df_met <- df_met |>
#   mutate(name = gsub(".*GetFitDaily_","", filename),
#          name = gsub(".csv", "", name),
#          site = strsplit(name, "_")[[1]][1],
#          pos = strsplit(name, "_")[[1]][2],
#          period = strsplit(name, "_")[[1]][3])

# get site and position info
df_met <- df_met |>
  dplyr::select(-`...1`) |>
  dplyr::mutate(source = str_extract(source, "([^/]+$)"),
                source = str_remove(source, ".csv")) |>
  separate(col = source, into = c("type", "site", "pos", "period"), sep = "_") |>
  select(-type)


# Load discharge, and K600 from Raymond equations data
df_kray <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                       "DAM_K600_Flux_all_Eq.xlsx"))

# Load pCO2 data from CO2SYS with uncertainty
df_co2sys <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                       "pCO2_uncertainty.xlsx"))

# Clean up column names
df_co2sys <- df_co2sys |>
  select(date = datetime,
         pCO2_ppmv = `pCO2 (ppmv)`, #partial pressure of CO2 in water
         dpCO2_ppmv = `u_pCO2 (ppmv)`, #standard error uncertainty
         CO2_mmolm3 = `CO2 (mmol/m3)`, #same but in mmol/m3
         dCO2_mmolm3 = `u_CO2 (mmol/m3)`,
         CO2_atm = `Atmospheric CO2 (mmol/m3)`)

# Quickly calculate Schmidt number CO2 (2 eqns in Raymond et al. 2012)
df_kray <- df_kray |>
  rename(temp = `Temp (C)`,
         Sc_CO2_1 = Schmidt_CO2) |>
  mutate(#Sc_O2_2 = 1801 - 120.1 * temp + 3.782 * temp^2 - 0.0476 * temp^3,
         #Sc_O2_1 = 1568 - 86.04 * temp + 2.142 * temp^2 - 0.0216 * temp^3,
         Sc_CO2_2 = 1742 - 91.24 * temp + 2.208 * temp^2 - 0.0219 * temp^3) |>
  mutate(#Sc_O2_mean = (Sc_O2_1 + Sc_O2_2) / 2,
         Sc_CO2_mean = (Sc_CO2_1 + Sc_CO2_2) / 2,
         #dSc_O2_mean = abs((Sc_O2_1 - Sc_O2_2)) / sqrt(2),
         dSc_CO2_mean = abs((Sc_CO2_1 - Sc_CO2_2)) / sqrt(2))

# only Dampierre
df_dam <- df_met |>
  filter(site == "dampierre", pos == "up")
  
# Count number of days with GPP and positive ER
sum(df_dam$ER_mean > 0, na.rm = T) #456, = 3.9%
sum(df_dam$ER_mean > 0 & df_dam$ER_2.5pct < 0, na.rm = T) #319, 70.0%
sum(df_dam$GPP_mean < 0, na.rm = T) # 1406 = 11.9%
sum(df_dam$GPP_mean < 0 & df_dam$GPP_97.5pct > 0, na.rm = T) #1251, 89.0%

# Correlation of GPP and ER to account for in error propagation
cor_ge <- cor(df_dam$GPP_mean, y = df_dam$ER_mean, use = "pairwise.complete.obs")

# Data cleaning -----------------------------------------------------------
# Step 1 Remove negative GPP and positive ER
# Set to 0 if biologically impossible, but 95% CI contains 0, otherwise NA
df_met_clean <- df_dam |>
  mutate(GPP_mean = if_else(GPP_mean < 0 & df_dam$GPP_97.5pct > 0, 0, GPP_mean),
         GPP_mean = if_else(GPP_mean < 0 & df_dam$GPP_97.5pct < 0, NA_real_, GPP_mean),
         ER_mean = if_else(ER_mean > 0 & df_dam$ER_2.5pct < 0, 0, ER_mean),
         ER_mean = if_else(ER_mean > 0 & df_dam$ER_2.5pct > NA_real_, 0, ER_mean)) |>
  mutate(across(contains("GPP"), ~ifelse(is.na(GPP_mean), NA_real_, .)),
         across(contains("ER"), ~ifelse(is.na(ER_mean), NA_real_, .)))

# Step 2 is to remove really low K600 values and interpolate them
# Based on observations, anything less than 0.5 1/d is questionable
df_k <- select(df_met_clean, date, contains("K")) |>
  mutate(K600_daily_mean = if_else(K600_daily_mean < 0.5, NA_real_, K600_daily_mean)) |>
  mutate(across(contains("K"), ~ifelse(is.na(K600_daily_mean), NA_real_, .))) |>
  imputeTS::na_kalman()

# Replace these new K values and recalculate flux
df_met_clean <- select(df_met_clean, -contains("K")) |>
  left_join(df_k)
  
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

# Get K600 error from Raymond et al. (2012), propagate error in Sc and K600 to KCO2
df_KCO2_ray <- df_kray |>
  select(date, contains("K600_Raymond")) |>
  pivot_longer(-date) |>
  group_by(date) |>
  summarize(K600_ray_mean = mean(value),
            dK600_ray_mean = sd(value)) |>
  left_join(select(df_kray, date, Sc_CO2_mean, dSc_CO2_mean), by = "date") |>
  mutate_with_error(KCO2_ray_mean ~ K600_ray_mean * (600 / Sc_CO2_mean)^(0.5)) |>
  ungroup()

# Calculate KCO2 from metabolism K600, propagate error in Sc and K600
df_KCO2_met <- select(df_met_clean, date, contains("K600_mean"), K600_2.5, K600_97.5) |>
  left_join(select(df_kray, date, Sc_CO2_mean, dSc_CO2_mean), by = "date") |>
  mutate_with_error(KCO2_met_mean ~ K600_mean * (600 / Sc_CO2_mean)^(0.5)) |>
  rename_with(~ str_replace(.x, 
                            pattern = "K600", 
                            replacement = "K600_met"), 
              matches("K600")) 

# Now we want to calculate CO2 fluxes with the two different K's
df_CO2 <- left_join(df_KCO2_met, df_KCO2_ray) |>
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

