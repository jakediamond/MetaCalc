# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate daily trophlux states, accounting for FCO2 chemical enhancement
# Date: 24 October 2022
# -------------------------------------

# Load libraries
library(patchwork)
library(tidyverse)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_lowpass_function.R"))
source(file.path("src", "000_error_propagation_function.R"))
# Load data---------------------------------------------------------------
# Hourly data
df <- readRDS(file.path("data", "hourly_data.RDS"))

# Clean the data ----------------------------------------------------------
# Get daily data first
df_d <- df |>
  select(date, datetime, FCO2_enh, dFCO2_enh) |>
  mutate(FCO2_hr = FCO2_enh /24,
         dFCO2_hr = dFCO2_enh / 24 * sqrt(10000)) |>
  group_by(date) |>
  summarize(FCO2_mean = sum(FCO2_hr, na.rm = T),
            dFCO2 = sqrt(sum(dFCO2_hr^2))) |>
  mutate(FCO2_2.5 = FCO2_mean - 1.96 * dFCO2,
         FCO2_97.5 = FCO2_mean + 1.96 * dFCO2) |>
  select(-dFCO2) |>
  left_join(distinct(df, date, NEP_mean, NEP_2.5, NEP_97.5))

df_clean <- df_d |>
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) |>
  group_by(type, val_type) |>
  arrange(type, val_type, date) |>
  nest() |>
  mutate(
    # Apply a lowpass function to get rid of excess noise
    filt = map(data, lowpass_fun)
    ) |>
  unnest(cols = filt) |>
  select(-data) |>
  pivot_wider(names_from = c(type, val_type), values_from = c(value, filtered)) |>
  ungroup() |>
  drop_na() |>
  left_join(distinct(df, date, Q_m3s)) |>
  arrange(date)

# Calculate trophlux states -------------------------------------------------------
# Define based on trophic state and co2 source/sink
df_nepco2 <- df_clean |>
  mutate(troph = if_else(filtered_NEP_mean > 0, "autotrophic", "heterotrophic"),
         sourcesink = if_else(filtered_FCO2_mean > 0, "source", "sink"),
         trophlux = paste(troph, sourcesink)) |>
  mutate(year = year(date),
         month = month(date))

# Save this data
saveRDS(df_nepco2, file.path("data", "daily_trophlux.RDS"))

df_nepco2 |>
  group_by(trophlux) |>
  summarize(n = n(),
            per = n()/nrow(df_nepco2))

# What occurs just before the river becomes a heterotrophic sink?
# First get some group information about each trophlux series
df_group <- df_nepco2 |>
  select(date, trophlux) |>
  group_by(trophlux) |>
  mutate(group = cumsum(c(1, diff(date) > 1))) |> # group by archeytpe events
  mutate(archgroup = paste0(trophlux, group)) |>
  mutate(arch_l = sequence(rle(archgroup)$lengths)) # length of events

# What is the mean length of events
df_group |>
  group_by(trophlux, archgroup) |>
  summarize(d = max(arch_l)) |>
  ungroup() |>
  group_by(trophlux) |>
  summarize(mean = mean(d),
            sd = sd(d))

r <- rep(0,nrow(df_nepco2))
df_hetsink <- within(df_group, 
                     identifier <- replace(r,
                                           sapply(which(trophlux=="heterotrophic sink" &
                                                          arch_l == 1),
                                                  `+`, -3:3),
                                           -3:3))
# 53/70 (76%) events are prefaced by autotrophic sinks
filter(ungroup(df_hetsink), identifier == -1) |>
  group_by(trophlux) |>
  summarize(n = n())


df_hetsink2 <- ungroup(df_group) %>%
  group_by(archgroup) %>%
  mutate(endgroup = if_else(arch_l == max(arch_l), 1, 0))


df_hetsinkend <- within(df_hetsink2, 
                     identifier <- replace(r,
                                           sapply(which(trophlux=="heterotrophic sink" &
                                                          endgroup == 1),
                                                  `+`, -3:3),
                                           -3:3))
# 34/70 (50%) events are followed by autotrophic sinks and heterotrophic sources
filter(ungroup(df_hetsinkend), identifier == 1) |>
  group_by(trophlux) |>
  summarize(n = n())
