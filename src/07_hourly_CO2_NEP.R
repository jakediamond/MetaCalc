# -------------------------------------
# Author: Jake Diamond
# Purpose: Look at hourly NEP and CO2 fluxes
# Date: 9 fevrier 2023
# -------------------------------------

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
# library(seacarb)
library(tidyverse)

# Load clean DO data
df_do <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS")) %>%
  filter(site == "dampierre", pos == "up")

# Quick clean
df_do <- df_do %>%
  mutate(datetime = floor_date(solar.time, "hours"))

# Load hourly CO2 data
df_co2 <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS")) %>%
  janitor::clean_names()

# Get data together
df <- select(df_do, do_clean = DO.obs, dosat_clean = DO.sat, 
             solartime = solar.time, datetime, light) %>%
  right_join(df_co2)

# Calculate hourly NEP 
df <- df %>%
  mutate(Sc_O2 = 1801 - 120.1 * temperature + 3.782 * temperature^2 - 0.0476 * temperature^3,
         Sc_CO2 = 1742 - 91.24 * temperature + 2.208 * temperature^2 - 0.0219 * temperature^3,
         KO2_md = k_co2_meta_1_d * (Sc_O2 / Sc_CO2)^-0.5,
         NEP = ((do_clean - lag(do_clean)) - (KO2_md / 24) * (dosat_clean - do_clean)) * depth)

# Add some additional information
df <- df %>%
  mutate(hour = hour(datetime),
         month = month(datetime),
         year= year(datetime))

# Diel plot of CO2 and NEP
ggplot(data = filter(df, date(datetime+hours(4)) == ymd(20020401))) + 
  geom_path(aes(x = datetime,
                y = -NEP), size= 1.5) +
  geom_path(aes(x = datetime,
                y = co2_flux_mmol_m2_d), size= 1.5, color = "red") +
  scale_color_viridis_c() +
  theme_bw() +
    geom_hline(yintercept = 0) +
  labs(y = expression("flux (mmol "~m^-2~h^-1*")"),
       x = "")
ggsave(file.path("results", "hourly_nep_co2_20200501.png"),
       dpi = 600,
       width = 16,
       height = 16,
       units = "cm")

# lag between min and max
df_dif <- df %>%
  mutate(date = date(datetime+ hours(4)),
         hour2 = hour + 24) %>%
  group_by(date, month, year) %>%
  select(hour2, nep = NEP, co2 = co2_flux_mmol_m2_d) %>%
  mutate(nep = -nep) %>%
  pivot_longer(cols = nep:co2) %>%
  group_by(date, month, year, name) %>%
  filter(value == max(value) | value == min(value)) %>%
  mutate(type = if_else(value == max(value), "max", "min")) %>%
  select(-value) %>%
  pivot_wider(names_from = "name", values_from = hour2, values_fn = mean) %>%
  mutate(dif = co2 - nep)

ggplot(data = filter(df_dif, type == "max"),
       aes(x = month, y = -dif)) +
  stat_summary() +
  theme_classic() +
  geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = seq(1,12,1))+
  labs(y = "hourly difference between co2 and nep",
       x = "month")
ggsave(file.path("results", "hourly_nep_co2_hourly_dif_max_stat_summary.png"),
       dpi = 600,
       width = 18,
       height = 14,
       units = "cm")

# Hysteresis between NEP and CO2
ggplot(data = filter(df, date(datetime) == ymd(20121011)),
       aes(x = NEP *1000 / 32,
           y = co2_flux_mmol_m2_d / 24,
           color = hco3_mmol_m3)) + 
  geom_point() +
  geom_path(size= 1.5) +
  scale_color_viridis_c() +
  theme_bw() +
  geom_abline(slope = -1, intercept = 0, color = "blue") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  labs(x = expression("NEP (mmol "*O[2]~m^-2~h^-1*")"),
       y = expression(CO[2]~"flux (mmol "*m^-2~h^-1*")"))

# ggsave(file.path("results", "hourly_nep_co2_hysteresis_hetsink_20050806.png"),
#        dpi = 600,
#        width = 18,
#        height = 14,
#        units = "cm")

ggplot(data = filter(df, month == 1, year == 2019),
       aes(x = NEP *1000 / 32,
           y = co2_flux_mmol_m2_d / 24,
           color = alkalinity)) + 
  geom_point() +
  geom_path(size= 1.5) +
  scale_color_viridis_c() +
  theme_bw() +
  geom_abline(slope = -1, intercept = 0, color = "blue") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  labs(x = expression("NEP (mmol "*O[2]~m^-2~h^-1*")"),
       y = expression(CO[2]~"flux (mmol "*m^-2~h^-1*")"))

ggsave(file.path("results", "hourly_nep_co2_hysteresis_example_201901.png"),
       dpi = 600,
       width = 18,
       height = 14,
       units = "cm")


# x = filter(df, month == 2, year == 2015)
x = filter(df, date(datetime) == ymd(20050806)) %>%
  summarize(nep = sum(NEP*1000/32),
            co2 = sum(co2_flux_mmol_m2_d/24))


ggplot(data = df, #filter(df, date(datetime) == ymd(20200812)),
       aes(x = o2ex,
           y = co2ex,
           color = hco3_mmol_m3)) + 
  geom_point(alpha = 0.5) +
  # geom_path(size= 1.5) +
  scale_color_viridis_c() +
  theme_bw() +
  geom_abline(slope = -1, intercept = 0, color = "blue") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  labs(x = expression(exO[2]~"(mmol "*O[2]~m^-3*")"),
       y = expression(exCO[2]~"(mmol "*m^-3*")"))
