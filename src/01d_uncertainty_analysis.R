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

# Load data ---------------------------------------------------------------
df <- readRDS(file = file.path("data", "03_CO2",
                                    "CO2_with_uncertainty_dampierre_up_updated_again.RDS")) 

# Load discharge, and K600 from Raymond equations data
df_q <- readxl::read_xlsx(file.path("data", "03_CO2", 
                                       "DAM_K600_Flux_all_Eq.xlsx")) %>%
  select(date, Discharge)


# Clean data --------------------------------------------------------------
# df_clean <- df %>%
#   mutate(GPP_2.5 = if_else(GPP_mean == 0, 0, GPP_2.5), #no uncertainty if GPP < 0
#          GPP_97.5 = if_else(GPP_mean == 0, 0, GPP_97.5),
#          ER_2.5 = if_else(ER_mean == 0, 0, ER_2.5), #or if ER > 0
#          ER_97.5 = if_else(ER_mean == 0, 0, ER_97.5)) %>%


# Look at uncertainty distributions ----------------------------------------
# First get into nice format for only variables of concern
df_cv <- df %>%
  mutate(ER_mean = abs(ER_mean),
         dCO2_flux = dCO2_flux*sqrt(1000)) %>%
  mutate_with_error(PR_mean ~ GPP_mean/ER_mean) %>%
  select(date, contains(c("mean", "d")), pCO2_ppmv, CO2_mmolm3, CO2_flux_mean = CO2_flux,
         -depth, -dSc_CO2_mean, -Sc_CO2_mean) %>%
  pivot_longer(cols = -date) %>%
  mutate(type = if_else(str_detect(name, "d"), "sd", "mean"),
         name2 = str_remove(name, "_mean"),
         name3 = str_extract(name2, "(?<=d).*")) %>%
  mutate(name4 = if_else(is.na(name3), name2, name3)) %>%
  select(date, name = name4, type, value) %>%
  pivot_wider(names_from = type, values_from = value)


# calculate cv for each day (sd/mean)
df_cv <- df_cv %>%
  mutate(cv = abs(sd / mean))

# get summary of that data
cv_summary <- df_cv %>% 
  filter(!is.infinite(cv),
         !(name %in% c("CO2_mmolm3", "KCO2_ray", "KCO2_met", "PR"))) %>%
  drop_na() %>%
  group_by(name) %>%
  summarize(median = round(median(cv), 2),
            mean = round(mean(cv), 2),
            sd = round(sd(cv), 2)) %>%
  mutate(lab = paste("median = ", median, "\nmean = ", mean, "\nsd = ", sd))

# Plot distributions of cv for each
df_cv %>%
  filter(!is.infinite(cv),
         !(name %in% c("CO2_mmolm3", "KCO2_ray", "KCO2_met", "PR"))) %>%
  drop_na() %>%
  ggplot(aes(x = cv, group = name)) +
  geom_histogram() +
  theme_classic() + 
  geom_vline(data = cv_summary,
             aes(xintercept = median), color = "red", size = 1.25) +
  geom_text(data = cv_summary, aes(label = lab), x = 1.2, y = 7500) +
  facet_wrap(~name) +
  scale_x_continuous(limits = c(0, 2)) +
  labs(x = "daily coefficient of variation (SD/mean)")
ggsave(filename = file.path("results", "uncertainty", "cv_distributions_updated_again.png"),
       dpi = 600,
       units = "cm",
       height = 16,
       width = 18.4)

# Plot time series with uncertainty ---------------------------------------
# Just plot GPP and ER first
# Get data in a nice format
df_p_met <- select(df, date, GPP_mean, GPP_2.5, GPP_97.5,
                   ER_mean, ER_2.5, ER_97.5) %>%
  # mutate(GPP_2.5 = if_else(GPP_mean == 0, 0, GPP_2.5),
  #        GPP_97.5 = if_else(GPP_mean == 0, 0, GPP_97.5),
  #        ER_2.5 = if_else(ER_mean == 0, 0, ER_2.5),
  #        ER_97.5 = if_else(ER_mean == 0, 0, ER_97.5)) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  arrange(type, date) %>%
  # mutate(mean = if_else(type == "GPP" & `2.5` < 0, NA_real_, mean),
  #        mean = if_else(type == "ER" & `97.5` > 0, NA_real_, mean)) %>%
  imputeTS::na_kalman()

# Smooth the data for a nicer plot
df_p_met_smooth <- df_p_met %>% 
  pivot_longer(cols = mean:`97.5`) %>%
  group_by(type, name) %>%
  arrange(type, name, date) %>%
  nest() %>%
  mutate(filt = map(data, lowpass_fun)) %>%
  unnest(cols = filt) %>%
  select(-data) %>%
  pivot_wider(names_from = name, values_from = c(value, filtered)) %>%
  ungroup() %>%
  drop_na() %>%
  left_join(select(df_q, date, Q = Discharge)) %>%
  arrange(date)

# Plotly plot
p_met_q <- plot_ly(data = df_p_met_smooth, x=~date,
                   color = ~type, 
                   colors = c("dark green", "black")) %>%
  add_trace(y= ~ filtered_mean, type = "scatter", mode='lines',
            showlegend = FALSE) %>%
  add_ribbons(ymin = ~filtered_2.5,
              ymax = ~ filtered_97.5,
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
         yaxis = list(title = TeX("\\text{flux (mmol } CO_{2} m^{-2} d^{-1})")),
         title = "GPP and ER (river's perspective)",
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_met_q)

# Now plot NEP and CO2 flux -----------------------------------------------
# Get data together and smooth a bit
df_p_NEP_CO2 <- df %>% 
  # mutate(ER_mean = abs(ER_mean)) %>%
  # mutate_with_error(PR_mean ~ GPP_mean/ER_mean) %>%
  # mutate(PR_2.5 = PR_mean - 1.96 * dPR_mean,
  #        PR_97.5 = PR_mean + 1.96 * dPR_mean) %>%
  mutate(across(contains("NEP"), ~.*-1)) %>%
  select(date, CO2_mean = CO2_flux, CO2_2.5 = CO2_flux_2.5, 
         CO2_97.5 = CO2_flux_97.5, NEP_mean, NEP_2.5, NEP_97.5) %>%
         # PR_mean, PR_2.5, PR_97.5) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  group_by(type, val_type) %>%
  arrange(type, val_type, date) %>%
  nest() %>%
  mutate(filt = map(data, lowpass_fun)) %>%
  unnest(cols = filt) %>%
  select(-data) %>%
  pivot_wider(names_from = c(type, val_type), values_from = c(value, filtered)) %>%
  ungroup() #%>%
  # pivot_wider(names_from = val_type, values_from = value) %>%
  # group_by(type) %>%
  # arrange(type, date) %>%
  # select(date, starts_with("CO2")) %>%
  # select(-CO2_mmolm3, -CO2_atm) %>%
  # pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "method", "val_type")) %>%
  # bind_rows(select(df_met_err, date, starts_with("NEP")) %>%
  #             mutate(across(where(is.numeric), ~.*-1000/32)) %>% #get NEP from atmosphere perspective in mmol
  #             pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  #             mutate(method = "met")) %>%
  # pivot_wider(names_from = val_type, values_from = value) %>%
  # group_by(type, method) %>%
  # arrange(type, method, date) %>%
  # imputeTS::na_kalman() %>%
  # mutate(across(where(is.numeric),
  #               ~stats::filter(., rep(1/5, 9), sides = 2))) %>% #smoothing for visualization
  # drop_na(mean, `2.5`, `97.5`)# %>%
  # ungroup()

w = df_p_NEP_CO2 %>% 
  filter(filtered_CO2_mean < filtered_NEP_mean)
w2 = filter(w, filtered_CO2_97.5 < filtered_NEP_97.5)
w3 = filter(w2, filtered_CO2_mean < 0, filtered_NEP_mean > 0)
1865/11822
# Plot NEP and CO2
p_NEP_CO2 <- plot_ly(data = df_p_NEP_CO2, 
                     x=~date,
                     color = ~type, 
                     colors = c("#1E88E5", "#FFC107")) %>%
  add_trace(y= ~ filtered_mean, type = "scatter", mode='lines') %>%
  add_ribbons(ymin = ~value_2.5,
              ymax = ~value_97.5,
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

x = df_p_NEP_CO2 %>%
  select(date, type, filtered_mean) %>%
  pivot_wider(names_from = type, values_from = filtered_mean) %>%
  filter(CO2 < 0 & NEP > 0)
y = filter(df, NEP_mean < 0 & CO2_flux < 0)
z = anti_join(y,x)

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



# Lowpass function --------------------------------------------------------

lowpass_fun <- function(data,
                        cutoff_frequency = 7) {
  # Impute all NAs so that there are none with kalman filtering
  data$value <- imputeTS::na_kalman(data$value)
  # Order the data, just in case
  data <- data[with(data, order(date)),]
  # Sampling rate [s^-1]
  sr <- 1 / 86400
  # Nyquist frequency = half the sampling rate
  nyq <- sr / 2
  # Cutoff frequency (days^-1)
  cutoff <- 1 / (cutoff_frequency * 60 * 60 * 24)
  # Normalized cutoff frequency for Butterworth filter
  W <- cutoff / nyq
  # Butterworth low-pass filter, digital, 2nd order
  myfilter <- signal::butter(2, W, type = 'low', plane = 'z')
  # Forward-reverse filter to remove phase-shift 
  # associated with Butterworth filter (must be in vector-form)
  vec <- as.vector(data$value)
  filtered <- signal::filtfilt(myfilter, vec)
  # Filtered data
  data$filtered <- filtered
  data <- data[with(data, order(date)), ]
  # rem <- sr / cutoff
  # data <- data[-c(1:rem, (nrow(data) - rem):nrow(data)),]
}
