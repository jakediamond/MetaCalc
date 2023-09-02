# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate daily trophlux states, accounting for FCO2 chemical enhancement
# Date: 24 October 2022
# -------------------------------------

# Load libraries
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

# Load daily mean discharge, pH and temperature, which are needed for the chemical enhancement
df_carb <- read_csv(file.path("data", "03_CO2", 
                                    "DAM_full_correct_daily2.csv")) %>%
  select(date = datetime, Discharge, temp = `Temp (C)`, pH) %>%
  mutate(date = mdy(date))

# Calculate daily mean chemical enhancement -----------------------------------
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

# Clean the data ----------------------------------------------------------
df_clean <- df %>%
  select(date, GPP_mean, GPP_2.5, GPP_97.5,
                   ER_mean, ER_2.5, ER_97.5, NEP_mean, NEP_2.5, NEP_97.5,
                   CO2_mean = CO2_flux, CO2_2.5 = CO2_flux_2.5, 
                   CO2_97.5 = CO2_flux_97.5) %>%
  left_join(select(df_enh, date, enh_mean = enh)) %>%
  # Calculate chemical enhancement
  mutate(CO2_meanenh = CO2_mean * enh_mean, 
         CO2_2.5enh = CO2_2.5 * enh_mean, 
         CO2_97.5enh = CO2_97.5 * enh_mean) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  group_by(type, val_type) %>%
  arrange(type, val_type, date) %>%
  nest() %>%
  mutate(
    # Apply a lowpass function to get rid of excess noise
    filt = map(data, lowpass_fun)
    ) %>%
  unnest(cols = filt) %>%
  select(-data) %>%
  pivot_wider(names_from = c(type, val_type), values_from = c(value, filtered)) %>%
  ungroup() %>%
  drop_na() %>%
  left_join(select(df_carb, date, Q = Discharge)) %>%
  arrange(date)

# Calculate trophlux states -------------------------------------------------------
# Get only the necessary data
df_nepco2 <- df_clean %>%
  select(date, Q, filtered_enh_mean, contains(c("NEP", "CO2")))

# Define based on archetype
df_nepco2 <- df_nepco2 %>%
  mutate(troph = if_else(filtered_NEP_mean > 0, "autotrophic", "heterotrophic"),
         sourcesink = if_else(filtered_CO2_meanenh > 0, "source", "sink"),
         archetype = str_c(troph, sourcesink, sep = "_")) %>%
  mutate(year = year(date),
         month = month(date),
         across(contains("NEP"), ~.*-1)) # get NEP from atmosphere perspective)

# Save this data
saveRDS(df_nepco2, file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

df_nepco2 %>%
  group_by(archetype) %>%
  summarize(n = n(),
            per = n()/nrow(df_nepco2))

# What occurs just before the river becomes a heterotrophic sink?
# First get some group information about each archetype series
df_archgroup <- df_nepco2 %>%
  select(date, archetype) %>%
  group_by(archetype) %>%
  mutate(group = cumsum(c(1, diff(date) > 1))) %>% # group by archeytpe events
  mutate(archgroup = paste0(archetype, group)) %>%
  mutate(arch_l = sequence(rle(archgroup)$lengths)) # length of events

# What is the mean length of events
df_archgroup |>
  group_by(archetype, archgroup) |>
  summarize(d = max(arch_l)) |>
  ungroup() |>
  group_by(archetype) |>
  summarize(mean = mean(d),
            sd = sd(d))

r <- rep(0,nrow(df_nepco2))
df_hetsink <- within(df_archgroup, 
                     identifier <- replace(r,
                                           sapply(which(archetype=="heterotrophic_sink" &
                                                          arch_l == 1),
                                                  `+`, -3:3),
                                           -3:3))
# 71/92 (77%) events are prefaced by autotrophic sinks
filter(ungroup(df_hetsink), identifier == -1) %>%
  group_by(archetype) %>%
  summarize(n = n())
