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
source(file.path("src", "000_co2_chem_enhance_functions.R"))

# Load data and clean---------------------------------------------------------------
# Load all the metabolism data for all sites, clean for site info, then get Dampierre
df_met <- readRDS(file.path("data", "dampierre_metabolism_clean.RDS"))

# Load K600 from metabolism
df_o2 <- readRDS(file.path("data", "hourly_o2_system.RDS"))

# Load carbonate system data with uncertainty for CO2
df_carb <- readRDS(file.path("data", "hourly_carbonate_system_v2.RDS"))

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

# Correlation of GPP and ER to account for in error propagation
cor_ge <- cor(df_met$GPP_mean, y = df_met$ER_mean, use = "pairwise.complete.obs")

# Data cleaning -----------------------------------------------------------
# Need to use a sampling of K600 for its standard deviation
# Otherwise we can have -K600, Which is impossible
df_met <- df_met |>
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
df_met_err <- df_met |>
  mutate(NEP_mean = GPP_mean + ER_mean,
         dNEP_mean = sqrt(dGPP_mean^2 + dER_mean^2 + 2 * cor_ge * dGPP_mean * dER_mean),
         NEP_2.5 = NEP_mean - 1.96*dNEP_mean, # estimate 95% credible interval for NEP
         NEP_97.5 = NEP_mean + 1.96*dNEP_mean)

# Calculate KCO2 from metabolism K600, propagate error in Sc and K600
df_KCO2 <- select(df_met, date, contains("K600_mean"), K600_2.5, K600_97.5) |>
  right_join(select(df_sc, date, datetime, Sc_CO2_mean, dSc_CO2_mean), by = "date") |>
  mutate_with_error(KCO2_mean ~ K600_mean * (600 / Sc_CO2_mean)^(0.5)) |>
  mutate(KCO2_2.5 = K600_2.5 * (600 / Sc_CO2_mean)^0.5,
         KCO2_97.5 = K600_97.5 * (600 / Sc_CO2_mean)^0.5)

# Now we calculate CO2 fluxes
# First get all the relevant data together
df_hr <- df_carb |>
  left_join(df_KCO2) |>
  left_join(select(df_o2, date, datetime, depth)) |>
  mutate(enh = chem_enh_fun(temp, pH, KCO2_mean, depth))

# Function for getting flux and uncertainty with monte carlo
fco2_fun <- function (depth, kco2, dkco2, ak, bk, co2, dco2, co2eq, n = 10000) {
  # distribution of CO2
  co2_dist <- rnorm(n, co2, dco2)
  # Distribution of KCO2 (1/d)
  k_dist <- truncnorm::rtruncnorm(n, a = ak, b = bk,
                                  mean = kco2, sd = dkco2)
  f_dist <-  k_dist * depth * (co2_dist - co2eq)
  f_mean <- mean(f_dist)
  f_se <- sd(f_dist) / sqrt(n)
  tibble(FCO2 = f_mean, dFCO2 = f_se)
}

# Do the calculations for mean FCO2 and it's distributions with Monte Carlo
# This takes a while and uses lots of memory (probably a better way)
df_FCO2 <- df_hr |>
  # tail(10) |>
  select(CO2_uM, dCO2_uM, CO2eq_uM, KCO2_mean, KCO2_2.5, KCO2_97.5, dKCO2_mean, depth) |>
  tidytable::mutate(
    fc = tidytable::pmap(list(depth, KCO2_mean, dKCO2_mean, KCO2_2.5, KCO2_97.5, CO2_uM, 
                  dCO2_uM, CO2eq_uM),
              fco2_fun))
  
# Save all data for future reference --------------------------------------
df_save <- df_FCO2 |>
  unnest(fc) |>
  # select(FCO2, dFCO2) |>
  left_join(df_hr) |>
  left_join(select(df_met_err, -contains("K600"))) |>
  mutate(across(contains(c("GPP", "ER", "NEP")), ~.*1000/32)) #from g O2 to mmmol

saveRDS(df_save, file = file.path("data", "nep_fco2_uncertainty_v2.RDS"))

