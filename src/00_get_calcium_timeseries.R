# 
# Purpose: To estimate calcium continuously
# Author: Jake Diamond
# Date: 27 February 2023
# 

# Load libraries
library(plotly)
library(lubridate)
library(broom)
library(tidyverse)

# Load functions ----------------------------------------------------------
source(file.path("src", "000_calcite_saturation_index_function.R"))

# Load data---------------------------------------------------------------
# Load water chemistry data
df_wq <- readxl::read_xlsx(file.path("data", "00_water chemistry", 
                                     "WQ_Loire_Moyenne.xlsx"))

# Hourly data at Dampierre
df_hr <- read_csv(file.path("data", "03_CO2", 
                            "DAM_full_correct_hourly.csv"))


# Estimate Calcium from Conductivity --------------------------------------
# Get daily info of important things
df_day <- df_hr %>%
  mutate(date = date(datetime)) %>%
  distinct(date, Discharge, Conductivity, Alkalinity) %>%
  left_join(mutate(df_hr, date = date(datetime)) %>%
              group_by(date) %>%
              summarize(HCO3_dam = mean(`HCO3 (mmol/m3)`, na.rm = T)) %>%
              ungroup())

# Quick clean
df_wq <- dplyr::mutate(df_wq,
                dplyr::across(everything(), ~ifelse(. == -1, NA_real_, .))) %>%
  filter(site_no %in% c("4046000",
                        "4046800",
                        "4048000"))

# Add discharge, temp. alk, and conductivity to daily chemistry data
df_wq_q <- df_day %>%
  right_join(mutate(df_wq,
                    date = ymd(paste(year, month, day))))

# Look at discharge effects on HCO3
p_qH <- ggplot(data = df_wq_q,
                aes(x = log(Discharge),
                    y = HCO3*1000 /(61.0168))) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  scale_y_continuous(limits = c(500, 3000)) +
  theme_classic(base_size = 10) +
  labs(x = expression("ln(Q) ("*m^3~s^{-1}*")"),
       y = expression(HCO[3]^{`-`}~"("*mu*M*")"))
p_qH

# Look at discharge effects on Ca
p_qCa <- ggplot(data = df_wq_q,
                aes(x = log(Discharge),
                    y = Ca*1000 /(40.078))) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  scale_y_continuous(limits = c(0, 1500)) +
  theme_classic(base_size = 10) +
  labs(x = expression("ln(Q) ("*m^3~s^{-1}*")"),
       y = expression(Ca^{`2+`}~"("*mu*M*")"))
p_qCa

# Look at discharge effects on SC
p_qSC <- ggplot(data = df_wq_q,
                aes(x = log(Discharge),
                    y = SC)) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  scale_y_continuous(limits = c(0, 1500)) +
  theme_classic(base_size = 10) +
  labs(x = expression("ln(Q) ("*m^3~s^{-1}*")"),
       y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_qSC

# Look at Ca vs SC
p_CaSC <- ggplot(data = df_wq_q,
                aes(x = SC,
                    y = Ca*1000 /(40.078))) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  scale_y_continuous(limits = c(0, 1800)) +
  scale_x_continuous(limits = c(100, 400)) +
  theme_classic(base_size = 10) +
  labs(x = expression(C[25]~"("*mu*S~cm^{-1}*")"),
       y = expression(Ca^{`2+`}~"("*mu*M*")"))
p_CaSC

# Look at HCO3 vs SC
p_HSC <- ggplot(data = filter(df_wq_q, between(month,5,9)),
                 aes(x = Conductivity,
                     y = HCO3*1000 /(61.0168))) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  scale_y_continuous(limits = c(500, 3000)) +
  scale_x_continuous(limits = c(100, 400)) +
  theme_classic(base_size = 10) +
  labs(x = expression(C[25]~"("*mu*S~cm^{-1}*")"),
       y = expression(HCO[3]^{`-`}~"("*mu*M*")"))
p_HSC

p <- (p_CaSC + p_HSC) / (p_qCa + p_qH) +plot_annotation(tag_levels = "a")

ggsave(plot = p,
       filename = file.path("results", "ca_hco3_q_sc.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12)


# Quick look at specific conductivity from NAIDES vs EDF
ggplot(data = df_wq_q,
       aes(x = log(Discharge),
           y = log(HCO3),
           color = HCO3)) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., 
                                                   ..adj.rr.label.., 
                                                   sep = "~~~~"))) +
  theme_bw()

# Robust regression indicates high agreement, suggesting interchangeability
summary(MASS::rlm(SC ~ Conductivity, data = df_wq_q))

# Want to estimate time series of Ca from data
summary(lm(Ca ~ HCO3 + log(Discharge) + SC, 
           data = df_wq_q %>%
             mutate(HCO3 = HCO3 *1000 /(61.0168),
                    Ca = Ca * 1000 / (40.078))))

# Get calcium from conductivity and bicarbonate measured at Dampierre
df_hr <- df_hr %>%
  mutate(Ca = -12.5861 + 0.0515 * Conductivity + 0.0115 * `HCO3 (mmol/m3)` +
           1.9636 * log(Discharge))

# See how that lines up with Naiades measurements
df_comp_ca <- df_hr %>%
  mutate(date = date(datetime)) %>%
  group_by(date) %>%
  summarize(Ca_est = mean(Ca)) %>%
  right_join(df_wq_q)

ggplot(data = df_comp_ca,
       aes(x = Ca_est,
           y = Ca)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# Underestimates by 35%
summary(MASS::rlm(Ca ~ Ca_est, data = filter(df_comp_ca, site_no %in% c("4046000",
                                                              "4046800",
                                                              "4048000"))))

# Get calcite preciptation estimates --------------------------------------
# Get SI and CCPP
df_ccpp <- df_hr %>%
  mutate(LSI = SI_cal_fun(temp = `Temp (C)`,
                          pH = pH,
                          cond = Conductivity,
                          Ca = Ca,
                          HCO3 = `HCO3 (mmol/m3)` / 1E6),
         CCPP = CCPP_fun(temp = `Temp (C)`,
                         pH = pH,
                         cond = Conductivity,
                         alk = Alkalinity,
                         Ca = Ca))

# Plot LSI and CCPP
p_LSI_CCPP <- plot_ly(data = df_ccpp,
                    x = ~datetime) %>%
  add_trace(y = ~ LSI, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y = ~ CCPP, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE, yaxis="y2") %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\color{#FFC107}{CCPP~(mg~CaCO_{3}~L^{-1})}")),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{LSI~(-)}")),
         title = TeX("\\text{Langelier Saturation Index and Calcium Carbonate Precipitation Potential}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_LSI_CCPP)


