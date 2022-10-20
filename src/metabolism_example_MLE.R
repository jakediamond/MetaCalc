# 
# Purpose: To estimate metabolism with MLE and raymond equations for K at Dampierre in 2018–2020
# Author: Jake Diamond and An Nguyen
# Date: October 14, 2022
# 

# Set working directory
setwd("C:/Users/diamo/Dropbox/Projects/Loire_DO")

# Load libraries
library(lubridate)
library(streamMetabolizer)
library(tidyverse)

# Look at what data inputs are needed for maximum likelihood estimation (MLE) model
metab_inputs("mle", "data")

# DO data load and clean --------------------------------------------------
# Load datetime, DO, light, depth, and temperature data and get right time zone
# DO should be called 'DO.obs' and temperature should be called 'temp.water'
df <- readRDS("Data/K600_estimates_nighttime_regression_Dampierre")

x = df@fit %>%
  filter(K600.daily > 0, ER.daily <0 , p.value < 0.05,
         year(date) == 2017)

# Force the correct time zone
df$datetime <- force_tz(df$datetime, "Etc/GMT+1")

# Prepare data for stream Metabolizer -------------------------------------
# Estimate DO saturation DO.sat
df$DO.sat <- ifelse(df$temp.water == 0,
                    0,
                    14.652 - 0.41022 * df$temp.water + 0.007991 * 
                      df$temp.water^2 - 0.000077774 * df$temp.water^3)

# Convert to solar time at Gien station
df$solar.time <- calc_solar_time(df$datetime, longitude = 2.5)

# Check that time zone it UTC
tz(df$solar.time)

# Caclculate light for periods without data (and to compare)
df$light_est <- calc_light(solar.time = df$solar.time,
                           latitude = 47.7,
                           longitude = 2.5)

# Replace light if it is missing with estimated light
df <- df %>%
  mutate(light = ifelse(is.na(light), light_est, light)) %>%
  select(-light_est)

# df should have 6 columns:
# solar.time, depth, light, DO.obs, DO.sat, temp.water
colnames(df)

# Load or calculate K600 values -------------------------------------------
# Daily K600 values with one column as "date" and the other as "K600.daily"
df_K600 <- read_csv()

# Configure the model -----------------------------------------------------
# Choose a model structure
# We choose a MLE model with K600 from empirical equations
mod <- mm_name(type = 'mle')

# Metabolism function ---------------------------------------
met_fun <- function(data, data_daily, name = mod){
  
  # Set the specifications
  specs <- specs(model_name = name)
  
  # Do the metabolism
  metab(specs = specs, 
        data = as.data.frame(data), 
        data_daily = as.data.frame(data_daily))
}

# Run the metabolism model  ---------------------------------
df_met_mle <- met_fun(df, df_K600)


