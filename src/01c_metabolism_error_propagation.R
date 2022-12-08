# 
# Purpose: Propagate error from metabolism to NEP, CO2 flux
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(plotly)
library(lubridate)
library(tidyverse)

# Load error propagation function
source(file.path("src", "00_error_propagation_function.R"))

# Load data and clean---------------------------------------------------------------
# Load metabolism data
df_met <- list.files(path = file.path("data", "02_metabolism", "new_metabolism"), 
                 pattern = "GetFitDaily",
                 full.names = TRUE) %>% 
  map_dfr(read_csv, show_col_types = FALSE)

# Load discharge, K600, and CO2 data
df_cq <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                      "DAM_K600_Flux_all_Eq.xlsx"))

# Quickly calculate Schmidt number for O2 and CO2 (2 eqns each in Raymond et al. 2012)
df_cq <- df_cq %>%
  rename(temp = `Temp (C)`,
         Sc_CO2 = Schmidt_CO2) %>%
  mutate(#Sc_O2_2 = 1801 - 120.1 * temp + 3.782 * temp^2 - 0.0476 * temp^3,
         Sc_O2 = 1568 - 86.04 * temp + 2.142 * temp^2 - 0.0216 * temp^3,
         dSc_CO2 = 0,
         dSc_O2 = 0) #%>%
         #Sc_CO2_2 = 1742 - 91.24 * temp + 2.208 * temp^2 - 0.0219 * temp^3) %>%
  # mutate(Sc_O2_mean = (Sc_O2_1 + Sc_O2_2) / 2,
  #        Sc_CO2_mean = (Sc_CO2_1 + Sc_CO2_2) / 2,
  #        dSc_O2_mean = abs((Sc_O2_1 - Sc_O2_2)) / sqrt(2),
  #        dSc_CO2_mean = abs((Sc_CO2_1 - Sc_CO2_2)) / sqrt(2))

# Remove negative GPP and positive ER, only select columns we want
df_met_clean <- df_met %>%
  filter(year(date) > 1992) %>% #data before 1992 unreliable
  mutate(GPP_mean = if_else(GPP_mean < 0, 0, GPP_mean),
         ER_mean = if_else(ER_mean > 0, 0, ER_mean)) %>%
  select(date, 
         GPP_mean, dGPP_mean = GPP_daily_sd, 
         GPP_2.5 = GPP_2.5pct, GPP_97.5 = GPP_97.5pct,
         ER_mean, dER_mean = ER_daily_sd, 
         ER_2.5 = ER_2.5pct, ER_97.5 = ER_97.5pct,
         K600_mean = K600_daily_mean, dK600_mean = K600_daily_sd,
         K600_2.5 = K600_daily_2.5pct, K600_97.5 = K600_daily_97.5pct)

# Propagate errors ----------------------------------------------
# Propagate error from GPP and ER to NEP
# Assume symmetric error around the mean (mostly true) and minimal covariance 
# between GPP and ER (not true, but the added effect is small)
df_met_err <- df_met_clean %>%
  mutate_with_error(NEP_mean ~ GPP_mean + ER_mean) %>%
  mutate(NEP_2.5 = NEP_mean - 1.96*dNEP_mean, # estimate 95% credible interval for NEP
         NEP_97.5 = NEP_mean + 1.96*dNEP_mean)

# a <- df_met_err %>%
#   mutate(gpp_dist = map2(GPP_mean, dGPP_mean, ~rnorm(1000, .x, .y)),
#          er_dist = map2(ER_mean, dER_mean, ~rnorm(1000, .x, .y)),
#          nep_dist = map2(gpp_dist, er_dist, ~.x +  .y),
#          nep_mean_mcmc = map_dbl(nep_dist, mean, na.rm = T),
#          nep_sd_mcmc = map_dbl(nep_dist, sd, na.rm = T))
# 
# mean(a$nep_sd_mcmc - a$dNEP_mean, na.rm = T)
# Get K600 error from Raymond et al. (2012), propagate error in Sc and K600 to KCO2
df_KCO2_ray <- df_cq %>%
  filter(year(date) > 1992) %>%
  select(date, contains("K600_Raymond")) %>%
  pivot_longer(-date) %>%
  group_by(date) %>%
  summarize(K600_ray_mean = mean(value),
            dK600_ray_mean = sd(value)) %>%
  left_join(select(df_cq, date, Sc_CO2, dSc_CO2), by = "date") %>%
  mutate_with_error(KCO2_ray_mean ~ K600_ray_mean/((600/Sc_CO2)^(-0.5))) %>%
  ungroup()

# Calculate KCO2 from metabolism K600, propagate error in Sc and K600
df_KCO2_met <- select(df_met_clean, date, contains("K600_mean")) %>%
  left_join(select(df_cq, date, Sc_CO2, dSc_CO2), by = "date") %>%
  mutate_with_error(KCO2_met_mean ~ K600_mean/((600/Sc_CO2)^(-0.5))) %>%
  rename_with(~ str_replace(.x, 
                            pattern = "K600", 
                            replacement = "K600_met"), 
              matches("K600")) 

# Now we want to calculate CO2 fluxes with the two different K's
df_CO2 <- left_join(df_KCO2_met, df_KCO2_ray) %>%
  left_join(select(df_cq, CO2_w = `CO2 (mmol/m3)`, CO2_a = `CO2_atm (mmol/m3)`,
                   depth = `Depth (m)`, date)) %>%
  mutate(dCO2_w = 0, dCO2_a = 0, ddepth = 0) %>% #no uncertainty in [CO2] and depth
  mutate_with_error(CO2_ray_mean ~ depth * (CO2_w - CO2_a) * KCO2_ray_mean) %>%
  mutate_with_error(CO2_met_mean ~ depth * (CO2_w - CO2_a) * KCO2_met_mean) %>%
  mutate(CO2_ray_2.5 = CO2_ray_mean - 1.96*dCO2_ray_mean, # estimate 95% credible interval for CO2 fluxes
         CO2_ray_97.5 = CO2_ray_mean + 1.96*dCO2_ray_mean,
         CO2_met_2.5 = CO2_met_mean - 1.96*dCO2_met_mean,
         CO2_met_97.5 = CO2_met_mean + 1.96*dCO2_met_mean)

# Plot time series with uncertainty ---------------------------------------
# Just plot GPP and ER first
# Get data in a nice format
df_p_met <- select(df_met_err, date, GPP_mean, GPP_2.5, GPP_97.5,
                   ER_mean, ER_2.5, ER_97.5) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  arrange(type, date) %>%
  imputeTS::na_interpolation()

# Smooth the data for a nicer plot
df_p_met_smooth <- df_p_met %>% 
  group_by(type) %>%
  mutate(across(where(is.numeric),
                ~stats::filter(., rep(1/5, 9), sides = 2))) %>%
  ungroup() %>%
  drop_na() %>%
  left_join(select(df_cq, date, Q = Discharge))

# Plotly plot
p_met_q <- plot_ly(data = df_p_met_smooth, x=~date,
                   color = ~type, 
                   colors = c("dark green", "black")) %>%
  add_trace(y= ~ mean, type = "scatter", mode='lines',
            showlegend = FALSE) %>%
  add_ribbons(ymin = ~`2.5`,
              ymax = ~`97.5`,
              showlegend = FALSE) %>%
  add_trace(y= ~ Q, type = "scatter", mode='lines',
            color = I("dark blue"),
            yaxis="y2",
            showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{discharge (} m^{3}  s^{-1})"),
                       tickfont = list(color = "darkblue"),
                       titlefont = list(color = "darkblue"),
                       range = list(0,10000)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\text{flux (g } O_{2} m^{-2} d^{-1})")),
         title = "GPP and ER (river's perspective)",
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_met_q)

# Now plot NEP and CO2 flux -----------------------------------------------
# Get data together and smooth a bit
df_p_NEP_CO2 <- df_CO2 %>% 
  select(date, starts_with("CO2")) %>%
  select(-CO2_w, -CO2_a) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "method", "val_type")) %>%
  bind_rows(select(df_met_err, date, starts_with("NEP")) %>%
              mutate(across(where(is.numeric), ~.*-1000/32)) %>% #get NEP from atmosphere perspective in mmol
              pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
              mutate(method = "met")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  group_by(type, method) %>%
  arrange(type, method, date) %>%
  mutate(across(where(is.numeric),
                ~stats::filter(., rep(1/5, 9), sides = 2))) %>% #smoothing for visualization
  drop_na(mean, `2.5`, `97.5`) %>%
  ungroup()

# Plot NEP and CO2
p_NEP_CO2 <- plot_ly(data = filter(df_p_NEP_CO2, method == "met"), x=~date,
                     color = ~type, 
                     colors = c("#1E88E5", "#FFC107")) %>%
  add_trace(y= ~ mean, type = "scatter", mode='lines') %>%
  add_ribbons(ymin = ~`2.5`,
              ymax = ~`97.5`,
              showlegend = FALSE) %>%
  # add_trace(y= ~ Q, type = "scatter", mode='lines',
  #           color = I("dark blue"), linetype = I("solid"),
  #           yaxis="y2",
  #           showlegend = FALSE) %>%
  layout(xaxis = list(title = ""),
         yaxis = list(title = TeX("\\text{flux (mmol } m^{-2} d^{-1})")),
         title = TeX("\\text{NEP [} O_{2} \\text{] and } CO_{2,K_{met}} \\text{flux (atm. perspective)}")) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_NEP_CO2)

# Plot CO2 from metabolism K and raymond K
p_CO2_ray_met <- plot_ly(data = filter(df_p_NEP_CO2, type == "CO2"), x=~date,
                     color = ~method, 
                     colors = c("#1E88E5", "#FFC107")) %>%
  add_trace(y= ~ mean, type = "scatter", mode='lines') %>%
  add_ribbons(ymin = ~`2.5`,
              ymax = ~`97.5`,
              showlegend = FALSE) %>%
  # add_trace(y= ~ Q, type = "scatter", mode='lines',
  #           color = I("dark blue"), linetype = I("solid"),
  #           yaxis="y2",
  #           showlegend = FALSE) %>%
  layout(xaxis = list(title = ""),
         yaxis = list(title = TeX("\\text{flux (mmol } m^{-2} d^{-1})")),
         title = TeX("CO_{2,K_{Ray}} \\text{ and } CO_{2,K_{met}} \\text{flux (atm. perspective)}")) %>%
  config(mathjax = "cdn")
htmltools::browsable(p_CO2_ray_met)

# Quick ratio of NEP/CO2 efflux
df_ratio <- df_CO2 %>% 
  select(date, CO2_ray_mean, dCO2_ray_mean, CO2_met_mean, dCO2_met_mean) %>%
  left_join(select(df_met_err, date, NEP_mean, dNEP_mean) %>%
              mutate(across(where(is.numeric), ~.*-1000/32))) %>% #get NEP from atmosphere perspective in mmol
  mutate_with_error(ray_mean ~ NEP_mean / CO2_ray_mean) %>%
  mutate_with_error(met_mean ~ NEP_mean / CO2_met_mean)

# for plotting, clean outliers and smooth
df_p_ratio <- df_ratio %>% 
  select(date, ray_mean, dray_mean, met_mean, dmet_mean) %>%
  mutate(across(where(is.numeric), ~if_else(abs(.) > 100, NA_real_, .))) %>%
  mutate(ray_2.5 = ray_mean - 1.96*dray_mean, # estimate 95% credible interval for NEP
         ray_97.5 = ray_mean + 1.96*dray_mean,
         met_2.5 = met_mean - 1.96*dmet_mean,
         met_97.5 = met_mean + 1.96*dmet_mean) %>%
  select(-dray_mean, -dmet_mean) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  group_by(type) %>%
  arrange(type, date) %>%
  mutate(across(where(is.numeric),
                ~imputeTS::na_kalman(.))) %>%
  mutate(across(where(is.numeric),
                ~stats::filter(., rep(1/5, 9), sides = 2))) %>%
  drop_na(mean, `2.5`, `97.5`) %>%
  ungroup()

# Plot ratios from different methods over time
p_ratios <- plot_ly(data = df_p_ratio, x=~date,
                     color = ~type, 
                     colors = c("#1E88E5", "#FFC107")) %>%
  add_trace(y= ~ mean, type = "scatter", mode='lines') %>%
  # add_ribbons(ymin = ~`2.5`,
  #             ymax = ~`97.5`,
  #             showlegend = FALSE) %>%
  # add_trace(y= ~ Q, type = "scatter", mode='lines',
  #           color = I("dark blue"), linetype = I("solid"),
  #           yaxis="y2",
  #           showlegend = FALSE) %>%
  layout(xaxis = list(title = ""),
         yaxis = list(title = "ratio NEP:CO2 efflux"))

htmltools::browsable(p_ratios)


# Plot CO2 over tiime
p_CO2 <- plot_ly(data = df_cq, x=~date) %>%
  add_trace(y= ~ `CO2 (mmol/m3)`, type = "scatter", mode='lines') %>%
  # add_ribbons(ymin = ~`2.5`,
  #             ymax = ~`97.5`,
  #             showlegend = FALSE) %>%
  add_trace(y= ~ Discharge, type = "scatter", mode='lines',
            color = I("dark blue"),
            yaxis="y2",
            showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{discharge (} m^{3}  s^{-1})"),
                       tickfont = list(color = "darkblue"),
                       titlefont = list(color = "darkblue"),
                       range = list(0,10000)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("[CO_{2}]\\text{(mmol} m^{-3})")),
         title = "CO2",
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_CO2)

# Plot all the plots together ---------------------------------------------

htmltools::browsable(htmltools::tagList(p_CO2, p_met_q, p_NEP_CO2, 
                                        p_CO2_ray_met, p_ratios))
# Old ---------------------------------------------------------------------
# p <- ggplot(data = df_p_met, aes(x = date)) +
#   geom_line(aes(y = mean,
#                 color = type)) +
#   geom_hline(yintercept = 0) +
#   geom_ribbon(aes(ymin = `2.5`, ymax = `97.5`, fill = type), alpha = 0.7) +
#   theme_bw() +
#   scale_color_manual(name = "flux", values = c("dark green", "black")) +
#   scale_fill_manual(name = "flux", values = c("dark green", "black")) +
#   labs(x = "", y = "metabolism (g O2/m2/d)") +
#   theme(axis.title.x = element_blank(),
#         legend.position = "none")
# 
# ggplotly(p)
# 
# p_smooth <- ggplot(data = df_p_met_smooth, aes(x = date)) +
#   geom_line(aes(y = mean,
#                 color = type)) +
#   geom_hline(yintercept = 0) +
#   geom_ribbon(aes(ymin = `2.5`, ymax = `97.5`, fill = type), alpha = 0.7) +
#   theme_bw() +
#   scale_color_manual(name = "flux", values = c("dark green", "black")) +
#   scale_fill_manual(name = "flux", values = c("dark green", "black")) +
#   labs(x = "", y = "metabolism (g O2/m2/d)") +
#   theme(axis.title.x = element_blank(),
#         legend.position = "none")
# 
# ggplotly(p_smooth)
# 
# htmltools::browsable(htmltools::tagList(ggplotly(p_smooth)))



df_use <- df_cq %>% 
  left_join(select(df_met_err, date, NEP_mean, dNEP_mean, GPP_mean, dGPP_mean, 
                   ER_mean, dER_mean, K600_mean, dK600_mean) %>%
              mutate(across(where(is.numeric), ~.*-1000/32))) %>% #get NEP from atmosphere perspective in mmol
  mutate_with_error(ray_mean ~ NEP_mean / CO2_ray_mean) %>%
  mutate_with_error(met_mean ~ NEP_mean / CO2_met_mean)

