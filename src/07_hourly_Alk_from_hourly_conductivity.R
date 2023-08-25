# -------------------------------------
# Author: Jake Diamond
# Purpose: Get hourly alkalinity data instead of daily from conductivity
# Date: 2023 March 02
# -------------------------------------

# Load libraries
library(plotly)
library(broom)
library(tidyverse)

# Load functions
source(file.path("src", "000_lowpass_function.R"))

# Load raw hourly data --------------------------------------------------------
# All time series from EDF, just want conductivity
df_cond <- readRDS(file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))

# daily data for discharge and temperature and pH
df_d <- readRDS(file.path("data", "03_CO2", "dampierre_all_daily_data.RDS"))

# Load best hourly measurements for Dampierre (not conductivity)
df_dam_hr <- read_csv(file.path("data", "03_CO2", "DAM_full_correct_hourly.csv"))

# Get needed data together
df <- select(df_dam_hr, year = Year, month = Month, datetime, Q_m3s = Discharge, 
             depth_m = Depth, temp = `Temp (C)`, pH, Alk_molm3 = Alkalinity, O2 = Oxy) %>%
  right_join(df_damup) %>%
  mutate(date = date(datetime))

# Look at regression
ggplot(data = df,
       aes(x = filtered,
           y = Alk_molm3)) + 
  geom_point()

df_mod <- df %>%
  group_by(date) %>%
  summarize(across(where(is.numeric), mean))

ggplot(data = df_mod,
       aes(x = filtered,
           y = Alk_molm3)) + 
  geom_point()

mod <- lm(Alk_molm3 ~ filtered*Q_m3s*O2+temp+pH+month, data = df_mod)
summary(mod)
plot(mod)
summary(MASS::rlm(Alk_molm3 ~ filtered*Q_m3s*O2+temp+pH+month, data = df_mod))
# Alk (mol/m3) = 0.0061 * Spc (uS/cm) + 0.1804

# Estimate hourly Alkalinity (mol/m3) = (mmol/L), need to divide by 1000 for mol/L or mol/kg
df <- mutate(df,
             Alk_molkg = (0.0061 * filtered + 0.1804) / 1000,
             Alk2 = (-0.069 + 0.00696 * filtered + 0.00083 * Q_m3s + 0.0186 * O2 +
                       0.0015 * temp + 0.0118 * pH - 0.0017 * month - 
                       4.1E-6 * filtered * Q_m3s - 8.56E-5 * filtered * O2 -
                       4.1E-5 * Q_m3s * O2 + 1.93E-7 * filtered * Q_m3s * O2)
                      / 1000)


# Estimate CO2 system
df_clean <- drop_na(df)
carb <- seacarb::carb(flag = 8,
                        var1 = df_clean$pH,
                        var2 = df_clean$Alk_molkg,
                        S = 0,
                        T = df_clean$temp)
df_CO2 <- df_clean %>%
  mutate(CO2 = carb$CO2 * 1E6,
         HCO3 = carb$HCO3 * 1E6,
         CO3 = carb$CO3 * 1E6,
         DIC = carb$DIC * 1E6)

# Get calcium from conductivity and bicarbonate measured at Dampierre
df_hr <- df_CO2 %>%
  mutate(Ca = -12.5861 + 0.0515 * filtered + 0.0115 * HCO3 +
           1.9636 * log(Q_m3s))


# Get SI and CCPP
df_ccpp <- df_hr %>%
  mutate(LSI = SI_cal_fun(temp = temp,
                          pH = pH,
                          cond = SpC,
                          Ca = Ca,
                          HCO3 = HCO3 / 1E6),
         CCPP = CCPP_fun(temp = temp,
                         pH = pH,
                         cond = SpC,
                         alk = Alk_molkg / 1000,
                         Ca = Ca))


# Estimate calcification
df_calc <- df_ccpp %>%
  mutate(alkdif_uM = Alk_molkg * 1E6 - lag(Alk_molkg * 1E6)) %>%
  mutate(calc = if_else(LSI > 0.4 & CCPP > 0.4,
                        -0.5 * alkdif_uM,
                        0))

# Look at daily calcification
df_calc_d <- df_calc %>%
  group_by(year, month, date) %>%
  summarize(calc = sum(calc, na.rm = T))

ggplot(data = df_calc_d,
       aes(x=date,
           y = calc)) +
  geom_point()
saveRDS(df_calc, file.path("data", "03_CO2", "hourly_calcite_precip.RDS"))

# Daily C fluxes and all variables
df_c <- readRDS(file.path("data", "03_CO2", 
                        "all_hourly_data_exDIC_calc.RDS")) %>%
  left_join(df_calc_d)
