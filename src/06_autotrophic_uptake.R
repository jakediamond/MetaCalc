# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot exO2 vs exCO2 and exDIC
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
# library(plotly)
# library(htmltools)
library(patchwork)
library(tidyverse)
# library(tidytable)
source(file.path("src", "000_carbonate_functions.R"))

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data_final.RDS")) |>
  mutate(trophlux = str_replace(trophlux, "_", " ")) |>
  drop_na(trophlux)
colnames(df)
df |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
  group_by(regime) |>
  summarize(co2 = median(CO2_uM, na.rm = T))

df |>
  filter(date == ymd(20170616)) |>
  ggplot(aes(x = NEP / depth_smooth,
             y = SpC - lag(SpC),
             color = exCO2)) +
  geom_point() +
  scale_color_viridis_c()

# determine if CO2 was depleted during the day
df_dep <- df |>
  mutate(dep = if_else(between(hr, 10, 20) & 
                         NEP > 0 &
                         exCO2 <= 0,
                       "depleted",
                       "not")) |>
  group_by(year, date, trophlux, troph, sourcesink) |>
  summarize(depleted = sum(dep == "depleted"),
            SI = max(SI))

df_dep |>
  ungroup() |>
  summarize(count = sum(depleted >=6, na.rm = T))

# Read in all the models
df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))

ggplot(filter(df_mods, between(slo_spc, -0.2, 0)), 
       aes(x = slo_spc,
           y = slo_hco3)) +
  stat_summary_bin()

# Check for uptake 
df_test <- df_mods |>
  left_join(df_dep) |>
  # Just look at high quality data
  # filter(r2_dic > 0.66) |>
  # Ignore winter days when GPP is too low
  # filter(!(trophlux == "heterotrophic source" & slo_dic > -0.6)) |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
  mutate(Z_hco3dic = abs((slo_hco3 - slo_dic)) / sqrt(se_hco3^2 + se_dic^2),
         p_hco3dic = pt(Z_hco3dic, 46, lower.tail = FALSE),
         Z_dic0.5 = abs((slo_dic - -0.5)) / sqrt(se_dic^2),
         p_dic0.5 = pt(Z_dic0.5, 23, lower.tail = FALSE),
         Z_hco30.5 = abs((slo_hco3 - -0.5)) / sqrt(se_hco3^2),
         p_hco30.5 = pt(Z_hco30.5, 23, lower.tail = FALSE),
         Z_spc0.13 = abs((slo_spc - -0.13)) / sqrt(se_spc^2),
         p_spc0.13 = pt(Z_spc0.13, 23, lower.tail = FALSE),
         Z_spc0.17 = abs((slo_spc - -0.17)) / sqrt(se_spc^2),
         p_spc0.17 = pt(Z_spc0.17, 23, lower.tail = FALSE),
         Z_spc0.04 = abs((slo_spc - -0.04)) / sqrt(se_spc^2),
         p_spc0.04 = pt(Z_spc0.04, 23, lower.tail = FALSE)) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
      depleted > 0 & p_hco3dic > 0.05  & p_dic0.5 < 0.05 ~ "HCO3",
      depleted > 0 & p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 > 0.05 ~ "HCO3",
      depleted > 0 & p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 < 0.05 & 
        between(slo_hco3, -1.2, -0.8) & between(slo_spc, -0.06, -0.03)~ "HCO3",
      depleted > 0 & p_dic0.5 > 0.05 & r2_dic > 0.7 ~ "CaCO3",
      depleted > 0 & p_dic0.5 > 0.05 & r2_dic <= 0.7 & between(slo_spc, -0.2, -0.04) ~ "CaCO3",
      depleted > 0 & between(slo_spc, -0.2, -0.06) & between(slo_hco3, -0.7, -0.4) ~ "CaCO3",
      TRUE ~ "CO2"
    )) 

x = df_test |>
  filter(year > 1994) |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
  filter(source == "HCO3")
# Summary of that data
df_test |>
  # filter(year > 1994) |>
  # filter(year > 1994, between(yday(date), 90 ,270)) |>
  group_by(source) |>
  # group_by(troph , source) |>
  # group_by(regime, source) |>
  # group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) 
  # arrange(desc(regime), troph, desc(sourcesink))


df_test |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
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


mean(df$CO2_uM, na.rm = T)
mean(df$pCO2_cor_uatm, na.rm = T)

abs((-0.507 - -0.5)) / sqrt(0.023^2)

t.test(scale(1:23)*0.023 * sqrt(23) + -0.507, mu = -0.5)

t.test(scale(1:23)*0.01 + -0.05, scale(1:n2)*sd2 + mean2, ...)

