# -------------------------------------
# Author: Jake Diamond
# Purpose: Plot all time series for SM
# Date: 8 March 2023
# -------------------------------------
# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(patchwork)
library(tidyverse)
library(tidytable)

# Load all data
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete_updateSpC_AT_Ca.RDS"))

# Load just daily NEP and co2 data
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

# Get data as monthly time series to show the data better
# Nice plotting format first
df_nepco2 <- df_nepco2 %>%
  distinct(date, year, month, NEP_mean = filtered_NEP_mean, 
         NEP_2.5 = filtered_NEP_2.5, NEP_97.5 = filtered_NEP_97.5,
         CO2_mean = filtered_CO2_meanenh,
         CO2_2.5 = filtered_CO2_2.5enh, CO2_97.5 = filtered_CO2_97.5enh) %>%
  mutate(year = year(date),
         month = month(date),
         across(contains("NEP"), ~.*-1)) %>%
  pivot_longer(cols = -c(month, year, date), names_sep = "_", 
               names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  mutate(wts = 1/abs(`2.5` - `97.5`)^2) %>%
  group_by(type) %>%
  arrange(type, date) %>%
  mutate(rollmean = slider::slide_dbl(
    .x = cur_data(),
    .f = ~ weighted.mean(
      x = .x$mean,
      w = .x$wts
    ), .before = 90
  ))

# Weighted average time series
p_ts_mo <- ggplot(data = df_nepco2,
                  aes(x = date,
                      color = type)) +
  geom_line(aes(y = rollmean),  linewidth = 1.2) + 
  scale_color_manual(name = "", values = c("#1E88E5", "#FFC107")) +
  theme_classic(base_size = 10) + 
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(0.05, 0.85),
        legend.key.height = unit(0.25, "cm"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = "",
       y = expression(FCO[2]~"("*mmol~m^{-2}~d^{-1}*")"))

p_ts_mo

# Excess O2 plot
p_do <- ggplot(data = df,
               aes(x = datetime)) +
  geom_line(aes(y = O2ex),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  labs(x = "",
       y = expression(exO[2]~"("*mmol~m^{-3}*")"))
p_do  

# Conductivity plot
p_c <- ggplot(data = df,
               aes(x = datetime)) +
  geom_line(aes(y = SpC),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "",
       y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_c 


# pH plot
p_pH <- ggplot(data = df,
              aes(x = datetime)) +
  geom_line(aes(y = pH),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "",
       y = "pH")
p_pH 

# pH plot
p_alk <- ggplot(data = mutate(df, alk = if_else(is.na(AT), Alk_molkg * 1e6, AT*1000)) %>%
                  filter(between(alk, 1200, 2300)),
               aes(x = datetime)) +
  geom_line(aes(y = alk),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "",
       y = expression(A[T]~"("*mmol~m^{-3}*")"))
p_alk 

# temp plot
p_t <- ggplot(data = df,
                aes(x = datetime)) +
  geom_line(aes(y = temp),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "",
       y = expression("temp. ("*degree~C*")"))
p_t 

# discharge plot
p_q <- ggplot(data = distinct(df, date, discharge),
              aes(x = date)) +
  geom_line(aes(y = discharge),  linewidth = 0.8) + 
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "",
       y = expression("discharge ("*m^{3}~s^{-1}*")"))
p_q 

p <- p_do / p_t / p_pH / p_c / p_alk / p_q / p_ts_mo 
# p

ggsave(plot = p,
       filename = file.path("results", "FigS1_timeseries.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 24)
