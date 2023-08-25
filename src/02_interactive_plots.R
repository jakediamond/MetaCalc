# 
# Purpose: Plot time series of CO2, NEP, carbonate system, etc.
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(tidyverse)

# Load data---------------------------------------------------------------
# NEED TO CLARIFY IF USED CONDUCTIVITY OR SPECIFIC CONDUCTANCE!!!!

# Hourly carbonate system
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS"))
df2 <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS")) %>%
  filter(site == "dampierre", pos == "down")

df_raw <- readRDS(file.path("data", "01_EDF", "raw", 
                            "Raw_cond_oxy_pH_temp_1992_2022.rds"))

dam <- filter(df_raw,
             site == "Dampierre",
             position == "upstream")
dam_l <- filter(df_raw, site == "Dampierre") %>%
  mutate(alpha = 0.0192+8E-5 * Temperature, # alpha constant for spc based on logger specs
         spc = Conductivity / (1 + alpha * (Temperature - 25))) %>%
  pivot_longer(cols = c(Conductivity:spc))

# Plot upstream downstream conductivity
p_cond <- plot_ly(data = filter(dam_l, name == "spc"),
                    x=~datetime,
                    y = ~value,
                    color = ~position) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}Conductivity~(mu*S~cm^{-1})"),
                      range = list(0, 460)),
         title = TeX("\\text{Conductivity upstream downstream}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_cond)


# Plot conductivity and discharge
p_cond_q <- plot_ly(data = left_join(df, select(dam, datetime, Conductivity)) %>%
                      mutate(alpha = 0.0192+8E-5 * Temperature, # alpha constant for spc based on logger specs
                             spc = Conductivity / (1 + alpha * (Temperature - 25))), # calculate spc based on logger specs), 
                      x=~datetime) %>%
  add_trace(y= ~ Conductivity, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ spc, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE, yaxis="y2") %>%
  # add_trace(y= ~ pH, type = "scatter", mode='lines', color = I("black"),
  #           yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{discharge}"),
                       range = list(0, 3600)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}Conductivity~(mu*S~cm^{-1})"),
                      range = list(0, 460)),
         title = TeX("\\text{Conductivity and discharge}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_cond_q)


# Plot carbonate system and pH
p_carb_sys <- plot_ly(data = df, 
                      x=~datetime) %>%
  add_trace(y= ~ `CO2 (mmol/m3)`, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ `HCO3 (mmol/m3)`/6, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE) %>%
  add_trace(y= ~ pH, type = "scatter", mode='lines', color = I("black"),
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{pH}"),
                       range = list(6, 11)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{[CO_{2}]}~or~\\color{#FFC107}{[HCO_3^{-}]}~(mmol~m^{-3})"),
                      range = list(0, 460)),
         title = TeX("\\text{carbonate system}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_carb_sys)


# Plot alkalinity and discharge
p_alk_q <- plot_ly(data = df, 
                      x=~datetime) %>%
  add_trace(y= ~ alk, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE) %>%
  add_trace(y= ~ discharge, type = "scatter", mode='lines', color = I("dark blue"),
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{discharge}~(m^{3}~s^{-1})"),
                       range = list(40, 3200)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#FFC107}{alkalinity}~(mmol~m^{-3})"),
                      range = list(1, 2.5)),
         title = TeX("\\text{alkalinity and discharge}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_alk_q)

# Plot O2 and CO2
p_C_O <- plot_ly(data = df2, 
                      x=~solar.time) %>%
  # add_trace(y= ~ `CO2 std (mmol/m3)`, type = "scatter", mode='lines', 
  #           color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ DO.obs / DO.sat * 100, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE,
            yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\color{#FFC107}{[O_2]}~(mmol~m^{-3})"),
                       range = list(50, 1000)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{[CO_{2}]}~(mmol~m^{-3})"),
                      range = list(0, 360)),
         title = TeX("O_2~and~CO_2"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_C_O)

# Plot O2 and CO2 excess
p_C_Oex <- plot_ly(data = df, 
                 x=~datetime) %>%
  add_trace(y= ~ CO2ex, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ O2ex, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE,
            yaxis="y2", showlegend = FALSE) %>%
  # add_trace(y= ~ pH, type = "scatter", mode='lines', color = I("black"),
  #           yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\color{#FFC107}{[exO_2]}~(\\mu M)"),
                       range = list(-270, 400)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{[exCO_{2}]}~(\\mu M)"),
                      range = list(-25, 450)),
         title = TeX("excess~O_2~and~CO_2"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_C_Oex)


# Plot alkalinity and pH
p_alk_carb <- plot_ly(data = df, 
                    x=~datetime) %>%
  add_trace(y= ~ Alkalinity, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ pH, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE,
            yaxis="y2", showlegend = FALSE) %>%
  # add_trace(y= ~ pH, type = "scatter", mode='lines', color = I("black"),
  #           yaxis="y2", showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\color{#FFC107}{pH}"),
                       range = list(7, 10.5)),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{Alkalinity}"),
                      range = list(1, 3)),
         title = TeX("Alkalinity and pH"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

browsable(p_alk_ph)

# Plot all the plots together ---------------------------------------------

browsable(tagList(p_alk_carb, p_alk_q, p_carb_sys, p_C_O, p_C_Oex))


# Other plots -------------------------------------------------------------



ggplot(data = filter(df, month(datetime) == 8,
                     between(CO2ex, -40, 60),
                     between(O2ex, -100, 150)),
       aes(x = CO2ex,
           y = O2ex,
           color = `HCO3 (mmol/m3)`),
       alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = -1.2, intercept = 0) +
  # stat_smooth(method = "lm", se = FALSE) +
  # ggpubr::stat_regline_equation(label.x = 3, label.y = 32) +
  scale_color_viridis_c() +
  facet_wrap(~year(datetime))

ggsave(file.path("results", "exo2_vs_exco2_august_by_year.png"),
       dpi = 1200,
       width = 18,
       height = 16,
       units = "cm")

mod = lm(O2ex ~ CO2ex, data = filter(df_carb, month(datetime) == 7))
summary(mod)
