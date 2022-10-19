# 
# Purpose: To estimate metabolism at Dampierre in 2018–2020
# Author: Jake Diamond and An Nguyen
# Date: October 10, 2022
# 

# Set working directory
setwd("")

# Load libraries
library(lubridate)
library(streamMetabolizer)
library(tidyverse)

# Look at what data inputs are needed for bayes model
metab_inputs("bayes", "data")

# Discharge data load and clean -----------------------------------------------------
# Load daily discharge data, make sure discharge is called 'discharge.daily'
df_q <- readxl::read_xlsx("Data/Loire_DO/Dampierre_1994_prepared.xlsx") %>%
  rename(discharge.daily = discharge) %>%
  mutate(date = date(datetime)) %>%
  select(discharge.daily, date) %>%
  group_by(date) %>%
  summarize(discharge.daily = mean(discharge.daily, na.rm = T))

# Get rid of negative values
df_q$discharge.daily <- ifelse(df_q$discharge.daily < 0,
                               NA,
                               df_q$discharge.daily)

# DO data load and clean --------------------------------------------------
# Load datetime, DO, light, and temperature data and get right time zone
# DO should be called 'DO.obs' and temperature should be called 'temp.water'
df <- readxl::read_xlsx("Data/Loire_DO/Dampierre_1994_prepared.xlsx") %>%
  select(-discharge)

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

# Caclculate light for periods without data (and to compare)
df$light_est <- calc_light(solar.time = df$solar.time,
                           latitude = 47.7,
                           longitude = 2.5)
df <- df %>%
  mutate(light = ifelse(is.na(light), light_est, light)) %>%
  select(-light_est)

# Calculate depth
depth <- df_q %>%
  mutate(depth = 0.134 * discharge.daily^0.4125)

# Get rid of datetime
df$datetime <- NULL

# Combine depth (m) with streamMetabolizer data
df <- depth %>%
  right_join(df %>%
               mutate(date = date(solar.time))) %>%
  select(-date)

# Now we should have two dataframes
# df should have 6 columns:
# solar.time, depth, light, DO.obs, DO.sat, temp.water
colnames(df)
# df_q should have 2 columns:
# date, discharge.daily
colnames(df_q)

# Configure the model -----------------------------------------------------
# Choose a model structure
# We choose a Bayesian model with both observation error and process error
# We will pool K600
bayes_mod <- mm_name(type = 'bayes', 
                     pool_K600 = 'binned', 
                     err_obs_iid = TRUE, 
                     err_proc_iid = TRUE)
bayes_mod

# Metabolism function ---------------------------------------
met_fun <- function(data, data_q, bayes_name = bayes_mod){
  # Calculate the natural-log-space centers of the discharge bins
  # These are the bins for the time frame of 
  # Use the width method as in the help file with with = 0.8 log units
  brks <- calc_bins(vec = log(df_q$discharge.daily),
                    method = "width",
                    width = 0.8)$bounds
  # 
  # # Estimate the mean ln(k600) value for the river from O'Connor and direct 
  # # measurements with floating dome
  # # These are the hyperprior mean for k600 in log space 
  k6 <- 4
  # 
  # # Same for standard deviation, super tight prior
  k6_sd <- 2
  # 
  # Set the specifications
  bayes_specs <- specs(model_name = bayes_name,
                       burnin_steps = 1000,
                       saved_steps = 500
                       , K600_lnQ_nodes_centers = brks
                       , K600_lnQ_nodes_meanlog = rep(k6,
                                                      length(brks))
                       , K600_lnQ_nodes_sdlog = rep(k6_sd,
                                                    length(brks))
                       # , GPP_daily_mu = 8 
                       # , GPP_daily_sigma = 6
  )
  
  # Do the metabolism
  metab(specs = bayes_specs, 
        data = as.data.frame(data), 
        data_daily = as.data.frame(data_q))
}

# Run the metabolism model  ---------------------------------
df_met2 <- met_fun(df, df_q)
x = get_params(df_met)

y = get_fit(df_met)

streamMetabolizer::plot_distribs(df_met, "GPP_daily")


x2 = get_params(df_met2)

y2 = get_fit(df_met2)

streamMetabolizer::plot_distribs(df_met, "K600_lnQ_nodes")

plot(x2$K600.daily, x$K600.daily)

plot(x2$K600.daily, x2$ER.daily)

plot(x$K600.daily, x$ER.daily)

write_csv2(x2, "Z:/MetaCalc/dampierre_1994_new_metabolism_biggerK600_french.csv")
