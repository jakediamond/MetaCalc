# 
# Purpose: Analyze errors associated with mean estimates for FCO2 and metabolism
# Author: Jake Diamond
# Date: 27 January 2023
# 

# Load libraries
library(tidyverse)

# Load uncertainty propagation code
source(file.path("src", "000_error_propagation_function.R"))

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS")) 
x <- select(df, GPP_mean, ER_mean, dGPP_mean, dER_mean)
y <- select(df, date, datetime, CO2_uM, dCO2_uM, KCO2_mean, dKCO2_mean, depth_m, 
            CO2eq_uM, K600_mean, dK600_mean, dSc_CO2_mean,
            FCO2, dFCO2) |>
  mutate(ddepth_m = 0, dCO2eq_uM = 0) |>
  mutate_with_error(fco2 ~ depth_m * KCO2_mean * (CO2_uM - CO2eq_uM))

# Clean data --------------------------------------------------------------
# Don't want to inflate uncertainty for GPP/ER for days without data
df <- df |>
  mutate(dGPP_mean = if_else(GPP_mean == 0, 0, dGPP_mean), #no uncertainty if GPP < 0
         dER_mean = if_else(ER_mean == 0, 0, dER_mean)) #or if ER > 0

# Get daily data for simplicity
df_d <- df |>
  select(date, datetime, FCO2, dFCO2, CO2_uM, dCO2_uM) |>
  mutate(FCO2_hr = FCO2 /24,
         dFCO2_hr = dFCO2 / 24 * sqrt(10000)) |>
  group_by(date) |>
  summarize(fCO2 = sum(FCO2_hr, na.rm = T),
            # was se not sd from monte carlo, multiply by sqrt(n)
            dfCO2 = sqrt(sum(dFCO2_hr^2)),
            # dfCO2 = se_magwt(FCO2, FCO2),
            CO2 = mean(CO2_uM),
            dCO2 = sqrt(sum(dCO2_uM^2))) |>
  left_join(distinct(df, date, K600 = K600_mean, dK600 = dK600_mean,
                     GPP = GPP_mean, dGPP = dGPP_mean,
                     ER = ER_mean, dER = dER_mean,
                     NEP = NEP_mean, dNEP = dNEP_mean))


z <- filter(y, date == ymd(19900101)) |>
  mutate(FCO2_hr = FCO2 /24,
         dFCO2_hr = dFCO2 / 24 * sqrt(10000))
sqrt(sum(z$dFCO2_hr^2))
# Look at uncertainty distributions ----------------------------------------
# First get into nice format for only variables of concern
df_cv <- df_d |>
  mutate(ER = abs(ER),
         fCO2 = abs(fCO2)) |> # was se not sd from monte carlo
  pivot_longer(cols = -date) |>
  mutate(type = if_else(str_detect(name, "d"), "sd", "mean"),
         name2 = str_remove(name, "_mean"),
         name3 = str_extract(name2, "(?<=d).*")) |>
  mutate(name4 = if_else(is.na(name3), name2, name3)) |>
  select(date, name = name4, type, value) |>
  pivot_wider(names_from = type, values_from = value)

# calculate cv for each day (sd/mean)
df_cv <- df_cv |>
  mutate(cv = abs(sd / mean))

# get summary of that data
cv_summary <- df_cv |> 
  filter(!is.infinite(cv)) |>
  drop_na() |>
  group_by(name) |>
  summarize(median = round(median(cv), 2),
            mean = round(mean(cv), 2),
            sd = round(sd(cv), 2)) |>
  mutate(lab = paste("median = ", median, "\nmean = ", mean, "\nsd = ", sd))

# Plot distributions of cv for each
df_cv |>
  filter(!is.infinite(cv)) |>
  drop_na() |>
  ggplot(aes(x = cv, group = name)) +
  geom_histogram() +
  theme_classic() + 
  geom_vline(data = cv_summary,
             aes(xintercept = median), color = "red", linewidth = 1.1) +
  geom_text(data = cv_summary, aes(label = lab), x = 1.2, y = 3500) +
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
  drop_na() |>
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


