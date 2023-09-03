# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot exO2 vs exCO2 and exDIC
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS"))
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS"))
df <- left_join(df,
                select(df_tf, date, troph, sourcesink, trophlux)) |>
  drop_na(troph)

# hourly exO2 vs del DIC each day
df_mod_O2DIC <- df |>
  drop_na(DIC_uM, exO2) |>
  group_by(year, month, date) |>
  filter(n() > 20) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$exO2 ~ .$DIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Regression results
df_O2DIC <- ungroup(df_mod_O2DIC) |>
  hoist(gl, "r.squared") |>
  select(-gl, -mod, -data) |>
  tidyr::unnest(cols = tid) |>
  filter(term == ".$DIC_uM") |>
  select(date, slo_dic = estimate, se_dic = std.error, p_dic = p.value, r2_dic =r.squared)

# hourly exO2 vs HCO3 each day
df_mod_O2HCO3 <- df |>
  drop_na(HCO3_uM, exO2) |>
  group_by(year, month, date) |>
  filter(n() > 20) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$exO2 ~ .$HCO3_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Regression results
df_O2HCO3 <- ungroup(df_mod_O2HCO3) |>
  hoist(gl, "r.squared") |>
  select(-gl, -mod, -data) |>
  tidyr::unnest(cols = tid) |>
  filter(term == ".$HCO3_uM") |>
  select(date, slo_hco = estimate, se_hco = std.error, p_hco = p.value, r2_hco =r.squared)

# hourly exO2 vs HCO3 each day
df_mod_CNEP <- df |>
  mutate(delC = SpC - lag(SpC)) |>
  drop_na(delC, NEP) |>
  group_by(year, month, date) |>
  filter(n() > 20) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$delC ~ .$NEP)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Regression results
df_CNEP <- ungroup(df_mod_CNEP) |>
  hoist(gl, "r.squared") |>
  select(-gl, -mod, -data) |>
  tidyr::unnest(cols = tid) |>
  filter(term == ".$NEP") |>
  select(date, slo_nep = estimate, se_nep = std.error, p_nep = p.value, r2_nep =r.squared)

# all results
df_all <- left_join(df_O2DIC, df_O2HCO3) |>
  left_join(df_CNEP)

saveRDS(df_all, file.path("data", "daily_regressions.RDS"))
