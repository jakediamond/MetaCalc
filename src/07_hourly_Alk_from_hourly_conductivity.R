# -------------------------------------
# Author: Jake Diamond
# Purpose: Get hourly alkalinity data instead of daily from conductivity
# Date: 2023 March 02
# -------------------------------------

# Load libraries
library(plotly)
library(lubridate)
library(broom)
library(tidyverse)

# Load functions
source(file.path("src", "000_lowpass_function.R"))

# Load raw hourly data --------------------------------------------------------
# All time series from EDF, just want conductivity
df_raw <- readRDS(file.path("data", "01_EDF", "raw", 
                            "Raw_cond_oxy_pH_temp_1992_2022.rds"))

# Only for Dampierre upstream
df_dam <- filter(df_raw,
                 site == "Dampierre",
                 position == "upstream") %>%
  select(datetime, SpC = Conductivity)

# Apply a lowpass filter to Conductivity
df_dam <- lowpass_fun(df_dam, cutoff_frequency = 8)

p <- plot_ly(data = df_dam,
            x  = ~datetime) %>%
  add_trace(y = ~ SpC, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y = ~ filtered, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE, yaxis="y2") %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\color{#FFC107}{CCPP~(mg~CaCO_{3}~L^{-1})}")),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{LSI~(-)}")),
         title = TeX("\\text{Langelier Saturation Index and Calcium Carbonate Precipitation Potential}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p)

# Load best hourly measurements for Dampierre (not conductivity)
df_dam_hr <- read_csv(file.path("data", "03_CO2", "DAM_full_correct_hourly.csv"))

# Get needed data together
df <- select(df_dam_hr, year = Year, month = Month, datetime, Q_m3s = Discharge, 
             depth_m = Depth, temp = `Temp (C)`, pH, Alk_molm3 = Alkalinity) %>%
  right_join(df_dam) %>%
  mutate(date = date(datetime))

# Look at regression
ggplot(data = df,
       aes(x = filtered,
           y = Alk_molm3)) + 
  geom_point()

summary(MASS::rlm(Alk_molm3 ~ filtered, data = df))
# Alk (mol/m3) = 0.0061 * Spc (uS/cm) + 0.1804

# Estimate hourly Alkalinity (mol/m3) = (mmol/L), need to divide by 1000 for mol/L or mol/kg
df <- mutate(df,
             Alk_molkg = (0.0061 * filtered + 0.1804) / 1000)


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
