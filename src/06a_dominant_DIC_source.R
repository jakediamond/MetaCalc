# -------------------------------------
# Author: Jake Diamond
# Purpose: to determine dominant DIC sources for autotrophs
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
# Hourly data
df <- readRDS(file.path("data", "hourly_data.RDS"))

df |>
  filter(NEP > 0, FCO2_enh < 0) |>
  summarize(p=median(pCO2_uatm, na.rm = T))

quantile(df$FCO2_enh, na.rm = T) / 1000 * 12
mean(df$FCO2_enh, na.rm = T) / 1000 * 12

quantile(df$NEP_mean, na.rm = T) / 1000 * 12
mean(df$NEP_mean, na.rm = T) / 1000 * 12

mean(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
sd(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12 / sqrt(283777)
median(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
quantile(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
mean(df$pCO2_uatm, na.rm = T)
sd(df$pCO2_uatm, na.rm = T)
median(df$pCO2_uatm, na.rm = T)

# daily trophlux data
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS")) |>
  drop_na(troph) |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic"))

df_tf |>
  group_by(troph) |>
  summarize(n = n())

df_y <- df_tf |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  filter(n() > 360) |>
  summarize(fco2_y = sum(filtered_FCO2_mean, na.rm = T),
            nep_y = sum(filtered_NEP_mean, na.rm = T))

mean(df_y$fco2_y) / 1000 * 12
sd(df_y$fco2_y) / 1000 * 12
quantile(df_y$fco2_y) / 1000 * 12

mean(df_y$nep_y) / 1000 * 12
sd(df_y$nep_y) / 1000 * 12
quantile(df_y$nep_y) / 1000 * 12

# Read in all the daily models
df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))
  
# Was CO2 depleted on a day? ----------------------------------------------
# determine if CO2 was depleted during the day
df_dep <- df |>
  mutate(exCO2 = CO2_uM - CO2eq_uM) |>
  mutate(dep = if_else(between(hr, 10, 20) & 
                         NEP > 0 &
                         exCO2 <= 0,
                       "depleted",
                       "not")) |>
  group_by(year, date) |>
  summarize(depleted = sum(dep == "depleted"),
            SI = max(SI))

# Quick summary of depleted days
df_dep |>
  ungroup() |>
  summarize(count = sum(depleted >6, na.rm = T))

# Calculate statistics on slopes ------------------------------------------
# Check for uptake 
df_test <- df_mods |>
  mutate(
    # Z score and p-value for if DIC slope is different than HCO3 (p > 0 is no difference)
    Z_hco3dic = abs((slo_hco3 - slo_dic)) / sqrt(se_hco3^2 + se_dic^2),
    p_hco3dic = pt(Z_hco3dic, 46, lower.tail = FALSE),
    # Z score and p-value for if DIC slope is = 0.5 (p> 0 is no difference)
    Z_dic0.5 = abs((slo_dic - -0.5)) / sqrt(se_dic^2),
    p_dic0.5 = pt(Z_dic0.5, 23, lower.tail = FALSE),
    # Z score and p-value for if HCO3 slope is = 0.5 (p> 0 is = 0.5)
    Z_hco30.5 = abs((slo_hco3 - -0.5)) / sqrt(se_hco3^2),
    p_hco30.5 = pt(Z_hco30.5, 23, lower.tail = FALSE),
    # Z score and p-value for if NEP slope is = 0.13
    Z_spc0.13 = abs((slo_spc - -0.13)) / sqrt(se_spc^2),
    p_spc0.13 = pt(Z_spc0.13, 23, lower.tail = FALSE),
    # Z score and p-value for if NEP slope is 0.17
    Z_spc0.17 = abs((slo_spc - -0.17)) / sqrt(se_spc^2),
    p_spc0.17 = pt(Z_spc0.17, 23, lower.tail = FALSE),
    # Z score and p-value for if NEP slope is 0.04 (p > 0 is = -0.04)
    Z_spc0.04 = abs((slo_spc - -0.04)) / sqrt(se_spc^2),
    p_spc0.04 = pt(Z_spc0.04, 23, lower.tail = FALSE))
  
df_up <- left_join(df_dep, df_test) |>
  left_join(df_tf) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
      depleted > 6 & p_hco3dic > 0.01  & p_dic0.5 < 0.01 ~ "HCO3",
      # DIC slope is different than HCO3 slope, but NEP slope = -0.04
      depleted > 6 & p_hco3dic < 0.01  & p_dic0.5 < 0.01 & p_spc0.04 > 0.01 ~ "HCO3",
      # depleted > 6 & p_hco3dic < 0.01  & p_dic0.5 < 0.01 & p_spc0.04 < 0.01 &
      #   between(slo_hco3, -1.2, -0.8) & between(slo_spc, -0.06, -0.03)~ "HCO3",
      depleted > 6 & SI > 0.4 & p_dic0.5 > 0.01 & r2_dic > 0.9 ~ "CaCO3", # need to add SI here
      depleted > 6 & SI > 0.4 & p_dic0.5 > 0.01 & r2_dic <= 0.9 & between(slo_spc, -0.2, -0.04) ~ "CaCO3",
      depleted > 6 & SI > 0.4 & between(slo_spc, -0.2, -0.06) & between(slo_hco3, -0.7, -0.4) ~ "CaCO3",
      TRUE ~ "CO2"
    )) 


# Summary of that data
df_up |>
  # data before 1994 is not reliable for SpC...no hourly data
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  group_by(source) |>
  # group_by(troph , source) |>
  # group_by(regime, source) |>
  # group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) 


# by years
df_y <- df_up |>
  filter(between(year, 1993, 2022), between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(year, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count))

# By year, then summarize
df_up |>
  filter(between(year, 1993, 2022), between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(year, troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  ungroup() %>%
  group_by(regime, troph, sourcesink, source) |> #
  summarize(mean = mean(x, na.rm = T) * 100,
           sd = sd(x, na.rm = T) * 100) |>
  mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))

# By year, then summarize but only regime
df_y |>
  select(-count) |>
  ungroup() %>%
  group_by(regime, source) |> #
  summarize(mean = mean(x, na.rm = T) * 100,
            sd = sd(x, na.rm = T) * 100) |>
  mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime))

# Overall
df_up |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
  group_by(regime, source) |>
  # group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) #|>
  # arrange(desc(regime), troph, desc(sourcesink))

# Overall for flux states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(sourcesink, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(sourcesink))

# Overall for trophic states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(troph)

# Overall for trophlux states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, troph) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(troph, desc(sourcesink))

# Overall for trophlux states and regime
df_up |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # ungroup() %>%
  # group_by(regime, troph, sourcesink, source) |> #
  # summarize(mean = mean(x, na.rm = T) * 100,
  #           sd = sd(x, na.rm = T) * 100) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))


# Overall for trophlux states
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # ungroup() %>%
  # group_by(regime, troph, sourcesink, source) |> #
  # summarize(mean = mean(x, na.rm = T) * 100,
  #           sd = sd(x, na.rm = T) * 100) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))

# Get mean pCO2 and Q for each year
df_qc <- df |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # mutate(wy = if_else(month > 8, year + 1, year)) |>
  group_by(year) |>
  summarize(qQ = median(Q_m3s),
            qC = median(pCO2_uatm))


mean(df$CO2_uM, na.rm = T)
mean(df$pCO2_cor_uatm, na.rm = T)

abs((-0.507 - -0.5)) / sqrt(0.023^2)

t.test(scale(1:23)*0.023 * sqrt(23) + -0.507, mu = -0.5)

t.test(scale(1:23)*0.01 + -0.05, scale(1:n2)*sd2 + mean2, ...)

