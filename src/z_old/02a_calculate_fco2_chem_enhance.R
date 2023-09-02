# -------------------------------------
# Author: Jake Diamond
# Purpose: Compare NEP and CO2 efflux over time
# Date: 24 October 2022
# -------------------------------------

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(patchwork)
library(tidyverse)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_co2_chem_enhance_functions.R"))
source(file.path("src", "000_lowpass_function.R"))
source(file.path("src", "000_error_propagation_function.R"))
# Load data---------------------------------------------------------------
# Daily C fluxes and all variables
df <- readRDS(file.path("data", "03_CO2", 
                        "CO2_with_uncertainty_dampierre_up_2023-02-21.RDS"))

# Load discharge, pH and temperature, which are needed for the chemical enhancement
df_carb <- read_csv(file.path("data", "03_CO2", 
                              "DAM_full_correct_daily2.csv")) %>%
  select(date = datetime, Discharge, temp = `Temp (C)`, pH) %>%
  mutate(date = mdy(date))

# Calculate chemical enhancement ------------------------------------------
df_enh <- select(df, date, KCO2 = KCO2_met_mean, depth) %>%
  left_join(df_carb) %>%
  mutate(
    # the chemical enhancement factor
    enh = chem_enh_fun(temp, pH, KCO2, depth)
  )

# Quick look
ggplot(data = df_enh,
       aes(x = date,
           y = enh)) +
  geom_point() +
  scale_y_log10()

ggplot(data = df_enh,
       aes(x = enh))+
  geom_histogram() +
  scale_x_log10()
