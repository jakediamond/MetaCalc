# -------------------------------------
# Author: Jake Diamond
# Purpose: To calculate the excess DIC compared to atmosphere
# Date: 25 February 2022
# -------------------------------------

# Load libraries
library(seacarb)
library(tidyverse)

# Load data---------------------------------------------------------------
# Hourly carbonate system
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS"))

df_carb <- carb(flag = 8, # input pH and ALK
                df$pH,
                df$Alkalinity / 1000, #mmol/L to mol/kg 
                S = 0,
                T = df$Temperature)
df <- tail(df, 200)
# Potential carbonate system
df_pot <-  seacarb::carb(flag = 24, # input eq pCO2 with atm. and ALK
                df$`Atmosphere CO2 (ppmv)`,
                df$Alkalinity / 1000, #mmol/L to mol/kg
                S = 0,
                pHscale = "F",
                T = df$Temperature)

df_pot2 <-  seacarb::carb(flag = 21, # input eq pCO2 with atm. and ALK
                         df$`Atmosphere CO2 (ppmv)`,
                         df$pH, #mmol/L to mol/kg
                         S = 0,
                         pHscale = "F",
                         T = df$Temperature)

df_test <- data.frame(temp = df$temp, Alk = x, atm_co2 = df_o2co2$atmosphere_co2_ppmv)
x = ifelse(is.na(df$Alk), df_o2co2$alkalinity/1000, is.na(df$Alk))

exDIC = (ifelse(is.na(df$DIC), NA_real_, df$DIC/1E9) - df_pot$DIC) * 1000

# 
# # Test sd of +-pH 0.1
# df_carb1 <- seacarb::carb(flag = 8, # input pH and ALK
#                 7,
#                 0.002, #mmol/L to mol/kg 
#                 S = 0,
#                 T = 15)
# df_carb2 <- seacarb::carb(flag = 8, # input pH and ALK
#                           6.9,
#                           0.002, #mmol/L to mol/kg 
#                           S = 0,
#                           T = 15)
# df_carb3 <- seacarb::carb(flag = 8, # input pH and ALK
#                           7.1,
#                           0.002, #mmol/L to mol/kg 
#                           S = 0,
#                           T = 15)
# 
# sd(c(df_carb1$CO2, df_carb2$CO2, df_carb3$CO2))


# add exDIC
df <- df %>%
  mutate(exDIC = (df_carb$DIC - df_pot$DIC) * 1000,
         delDIC = (df_carb$DIC - lag(df_carb$DIC)) * 1000,
         DIC = df_carb$DIC,
         delDO = O2_mmolm3 - lag(O2_mmolm3)) #mmol/kg

saveRDS(df, file.path("data", "03_CO2", "all_hourly_data_exDIC_calc.RDS"))

# The whole dang thing
filter(df, year(datetime) > 1992) %>%
  mutate(Q_breaks = cut(log(discharge, base = 10), 6)) %>%
  ggplot(aes(x = exDIC,
           y = -O2ex / 1000, #mmol/m3 to  mmol/kg, turn to AOU
           color = `HCO3 (mmol/m3)`),
       alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  facet_wrap(~Q_breaks) + 
  # geom_abline(slope = -1.2, intercept = 0) +
  stat_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c() +
  labs(x = expression("exDIC ("*mmol~kg^{-1}*")"), 
       y = expression("AOU ("*mmol~kg^{-1}*")"))

ggsave(file.path("results", "AOU_vs_exDIC_all_data_by_discharge_bin.png"),
       dpi = 1200,
       width = 28,
       height = 20,
       units = "cm")

# August data
ggplot(data = filter(df, month(datetime) == 8,
                     year(datetime) > 1992),
       aes(x = exDIC,
           y = -O2ex / 1000, #mmol/m3 to  mmol/kg, turn to AOU
           color = `HCO3 (mmol/m3)`),
       alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  # geom_abline(slope = -1.2, intercept = 0) +
  stat_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c() +
  facet_wrap(~year(datetime)) +
  labs(x = expression("exDIC ("*mmol~kg^{-1}*")"), 
       y = expression("AOU ("*mmol~kg^{-1}*")"))

ggsave(file.path("results", "AOU_vs_exDIC_august_by_year.png"),
       dpi = 1200,
       width = 28,
       height = 20,
       units = "cm")

#  April
ggplot(data = filter(df, month(datetime) == 4,
                     # date(datetime) == ymd(20190724),
                     year(datetime) > 1992),
       aes(x = exDIC,
           y = -O2ex / 1000, #mmol/m3 to  mmol/kg, turn to AOU
           color = `HCO3 (mmol/m3)`),
       alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  # geom_abline(slope = -1.2, intercept = 0) +
  stat_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c() +
  facet_wrap(~year(datetime)) +
  labs(x = expression("exDIC ("*mmol~kg^{-1}*")"), 
       y = expression("AOU ("*mmol~kg^{-1}*")"))

ggsave(file.path("results", "AOU_vs_exDIC_april_by_year.png"),
       dpi = 1200,
       width = 28,
       height = 20,
       units = "cm")

# Look at daily cycles
ggplot(data = filter(df, year(datetime) == 2019) %>%
                     group_by(month(datetime)) %>%
                       slice_head(n = 24),
       aes(x = exDIC,
           y = -O2ex / 1000, #mmol/m3 to  mmol/kg
           color = `HCO3 (mmol/m3)`),
       alpha = 0.5) +
  geom_point() +
  geom_text(aes(label = hour(datetime)), hjust = "outward") +
  geom_path() +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) + 
  # geom_abline(slope = -1.2, intercept = 0) +
  stat_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c() +
  facet_wrap(~date(datetime)) +
  labs(x = expression("exDIC ("*mmol~kg^{-1}*")"), 
       y = expression("AOU ("*mmol~kg^{-1}*")"))

ggsave(file.path("results", "AOU_vs_exDIC_hourly_2019_by_month.png"),
       dpi = 1200,
       width = 28,
       height = 20,
       units = "cm")

# AOU vs del DIC
df_mod <- df %>%
  mutate(date = date(datetime),
         year = year(datetime),
         month = month(datetime),
         AOU = -O2ex / 1000) %>%
  filter(year > 1992) %>%
  group_by(year, month, date) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$AOU ~ .$exDIC)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

df_res <- select(df_mod, -data,- mod) %>%
  unnest()

# Plot of daily slopes
df_res %>%
  filter(term == ".$exDIC") %>%
  ggplot(aes(x = date,
             y = estimate)) + 
  geom_point(aes(alpha = r.squared)) + 
  scale_y_continuous(limits = c(-3,3)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%y",
               name = "") +
  stat_smooth() +
  theme_bw() +
  labs(y = "AOU vs exDIC slope")


ggsave(file.path("results", "AOU_vs_exDIC_slope_daily.png"),
       dpi = 1200,
       width = 28,
       height = 16,
       units = "cm")

# Changes in monthly slope over time
df_res %>%
  filter(term == ".$exDIC",
         between(estimate, -10, 10)) %>%
  ggplot(aes(x = year,
             y = estimate)) + 
  stat_summary() + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  # scale_y_continuous(limits = c(-3,3)) +
  # scale_x_date(date_breaks = "1 year", date_labels = "%y",
  #              name = "") +
  stat_smooth() +
  facet_wrap(~month) + 
  theme_bw() +
  labs(y = "AOU vs exDIC slope")

ggsave(file.path("results", "AOU_vs_exDIC_slope_yearly_by_month.png"),
       dpi = 1200,
       width = 28,
       height = 16,
       units = "cm")
