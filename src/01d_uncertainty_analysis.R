# 
# Purpose: Analyze errors associated with mean estimates for CO2 and metabolism
# Author: Jake Diamond
# Date: 27 January 2023
# 

# Load libraries
library(plotly)
library(lubridate)
library(tidyverse)

# Load error propagation function
source(file.path("src", "000_error_propagation_function.R"))
source(file.path("src", "000_lowpass_function.R"))

# Load data ---------------------------------------------------------------
df <- readRDS(file = file.path("data", "03_CO2",
                                    "CO2_with_uncertainty_dampierre_up_2023-02-21.RDS")) 

# Load discharge, and K600 from Raymond equations data
df_q <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                       "DAM_K600_Flux_all_Eq.xlsx")) |>
  select(date, Discharge)

# Clean data --------------------------------------------------------------
# df_clean <- df |>
#   mutate(GPP_2.5 = if_else(GPP_mean == 0, 0, GPP_2.5), #no uncertainty if GPP < 0
#          GPP_97.5 = if_else(GPP_mean == 0, 0, GPP_97.5),
#          ER_2.5 = if_else(ER_mean == 0, 0, ER_2.5), #or if ER > 0
#          ER_97.5 = if_else(ER_mean == 0, 0, ER_97.5)) |>


# Look at uncertainty distributions ----------------------------------------
# First get into nice format for only variables of concern
df_cv <- df |>
  mutate(ER_mean = abs(ER_mean),
         dCO2_flux = dCO2_flux*sqrt(1000)) |>
  # mutate_with_error(PR_mean ~ GPP_mean/ER_mean) |>
  select(date, contains(c("mean", "d")), pCO2_ppmv, CO2_mmolm3, CO2_flux_mean = CO2_flux,
         -depth, -dSc_CO2_mean, -Sc_CO2_mean) |>
  pivot_longer(cols = -date) |>
  mutate(type = if_else(str_detect(name, "d"), "sd", "mean"),
         name2 = str_remove(name, "_mean"),
         name3 = str_extract(name2, "(?<=d).*")) |>
  mutate(name4 = if_else(is.na(name3), name2, name3)) |>
  select(date, name = name4, type, value) |>
  pivot_wider(names_from = type, values_from = value) |>
  mutate(name = case_match(name,
                           "K600_met" ~ "K600",
                           "CO2_flux" ~ "fCO2",
                           "pCO2_ppmv" ~ "pCO2",
                           .default = name))

# calculate cv for each day (sd/mean)
df_cv <- df_cv |>
  mutate(cv = abs(sd / mean))

# get summary of that data
cv_summary <- df_cv |> 
  filter(!is.infinite(cv),
         !(name %in% c("CO2_mmolm3", "KCO2_ray", "KCO2_met", "PR"))) |>
  drop_na() |>
  group_by(name) |>
  summarize(median = round(median(cv), 2),
            mean = round(mean(cv), 2),
            sd = round(sd(cv), 2)) |>
  mutate(lab = paste("median = ", median, "\nmean = ", mean, "\nsd = ", sd))

# Plot distributions of cv for each
df_cv |>
  filter(!is.infinite(cv),
         !(name %in% c("CO2_mmolm3", "KCO2_ray", "KCO2_met", "PR"))) |>
  drop_na() |>
  ggplot(aes(x = cv, group = name)) +
  geom_histogram() +
  theme_classic() + 
  geom_vline(data = cv_summary,
             aes(xintercept = median), color = "red", size = 1.25) +
  geom_text(data = cv_summary, aes(label = lab), x = 1.2, y = 7500) +
  facet_wrap(~name) +
  scale_x_continuous(limits = c(0, 2)) +
  labs(x = "daily coefficient of variation (SD/mean)")
ggsave(filename = file.path("results", "uncertainty", "cv_distributions.png"),
       dpi = 300,
       units = "cm",
       height = 12,
       width = 18.4)
# Calculate central tendency with different measures ----------------------
df_cent <- df_cv |>
  mutate(wt = 1 / sd^2) |>
  group_by(name) |>
  summarize(mn = mean(mean, na.rm = T),
            med = median(mean, na.rm = T),
            wtdmn = weighted.mean(mean, wt, na.rm = T),
            wtdmed = Hmisc::wtd.quantile(mean, wt, probs = 0.5, normwt = T))

# Estimate how often NEP or fCO2 were not different than 0 ----------------
# distributions of the variables based on sampling 1000 times
df_dist <- df_cv |>
  select(-cv) |>
  filter(name %in% c("NEP", "fCO2")) |>
  pivot_wider(names_from = name, values_from = c(mean, sd)) |>
  mutate(fco2_dist = map2(mean_fCO2, sd_fCO2, ~rnorm(1000, .x, .y)),
         nep_dist = map2(mean_NEP, sd_NEP, ~rnorm(1000, .x, .y)))

# T tests
df_test <- drop_na(df_dist) |>
  mutate(co2_t = map(fco2_dist, t.test),
         nep_t = map(nep_dist, t.test))

# Tidy output
df_tidy <- df_test |>
  mutate(co2_p = map(co2_t, ~.$p.value),
         nep_p = map(nep_t, ~.$p.value)) |>
  select(date, co2_p, nep_p)

sum(df_tidy$co2_p > 0.01)
sum(df_tidy$nep_p > 0.01)

# Plot time series with uncertainty ---------------------------------------
# Just plot GPP and ER first
# Get data in a nice format
df_p_met <- select(df, date, GPP_mean, GPP_2.5, GPP_97.5,
                   ER_mean, ER_2.5, ER_97.5) |>
  # mutate(GPP_2.5 = if_else(GPP_mean == 0, 0, GPP_2.5),
  #        GPP_97.5 = if_else(GPP_mean == 0, 0, GPP_97.5),
  #        ER_2.5 = if_else(ER_mean == 0, 0, ER_2.5),
  #        ER_97.5 = if_else(ER_mean == 0, 0, ER_97.5)) |>
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) |>
  pivot_wider(names_from = val_type, values_from = value) |>
  arrange(type, date) |>
  # mutate(mean = if_else(type == "GPP" & `2.5` < 0, NA_real_, mean),
  #        mean = if_else(type == "ER" & `97.5` > 0, NA_real_, mean)) |>
  imputeTS::na_kalman()

# Smooth the data for a nicer plot
df_p_met_smooth <- df_p_met |> 
  pivot_longer(cols = mean:`97.5`) |>
  group_by(type, name) |>
  arrange(type, name, date) |>
  nest() |>
  mutate(filt = map(data, lowpass_fun)) |>
  unnest(cols = filt) |>
  select(-data) |>
  pivot_wider(names_from = name, values_from = c(value, filtered)) |>
  ungroup() |>
  drop_na() |>
  left_join(select(df_q, date, Q = Discharge)) |>
  arrange(date)

# Plotly plot
p_met_q <- plot_ly(data = df_p_met_smooth, x=~date,
                   color = ~type, 
                   colors = c("dark green", "black")) |>
  add_trace(y= ~ filtered_mean, type = "scatter", mode='lines',
            showlegend = FALSE) |>
  add_ribbons(ymin = ~filtered_2.5,
              ymax = ~ filtered_97.5,
              showlegend = FALSE) |>
  add_trace(y= ~ Q, type = "scatter", mode='lines',
            color = I("dark blue"),
            yaxis="y2",
            showlegend = FALSE) |>
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{discharge (} m^{3}  s^{-1})"),
                       tickfont = list(color = "darkblue"),
                       titlefont = list(color = "darkblue"),
                       range = list(0,10000)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\text{flux (mmol } CO_{2} m^{-2} d^{-1})")),
         title = "GPP and ER (river's perspective)",
         margin = list(r = 50)) |>
  config(mathjax = "cdn")

htmltools::browsable(p_met_q)

# Now plot NEP and CO2 flux -----------------------------------------------
# Get data together and smooth a bit
df_p_NEP_CO2 <- df |> 
  mutate(across(contains("NEP"), ~.*-1)) |>
  select(date, CO2_mean = CO2_flux, CO2_2.5 = CO2_flux_2.5, 
         CO2_97.5 = CO2_flux_97.5, NEP_mean, NEP_2.5, NEP_97.5) |>
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) |>
  group_by(type, val_type) |>
  arrange(type, val_type, date) |>
  nest() |>
  mutate(filt = map(data, lowpass_fun)) |>
  unnest(cols = filt) |>
  select(-data) |>
  pivot_wider(names_from = c(val_type), values_from = c(value, filtered)) |>
  ungroup()

# Plot NEP and CO2
p_NEP_CO2 <- plot_ly(data = df_p_NEP_CO2, 
                     x=~date,
                     color = ~type, 
                     colors = c("#1E88E5", "#FFC107")) |>
  add_trace(y= ~ filtered_mean, type = "scatter", mode='lines') |>
  add_ribbons(ymin = ~value_2.5,
              ymax = ~value_97.5,
              showlegend = FALSE) |>
  layout(xaxis = list(title = ""),
         yaxis = list(title = TeX("\\text{flux (mmol } m^{-2} d^{-1})")),
         title = TeX("\\text{NEP [} O_{2} \\text{] and } CO_{2,K_{met}} \\text{flux (atm. perspective)}")) |>
  config(mathjax = "cdn")

htmltools::browsable(p_NEP_CO2)


