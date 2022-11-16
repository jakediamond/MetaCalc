# -------------------------------------
# Author: Jake Diamond
# Purpose: Get temperature, discharge, and radiation data together
# Date: 27 October 2022
# -------------------------------------
# load libraries
library(lubridate)
library(plotly)
library(tidyverse)
library(tidytable)

# Load all discharge data -------------------------------------------------
# Data directory
data_dir <- file.path("data", "01a_discharge")

# files
files <- fs::dir_ls(data_dir, regexp = "\\.txt$")

# Load discharge data
df_q <- files %>% 
  map_dfr(read_delim, skip = 41, col_names = TRUE, delim = ";",
          trim_ws = TRUE,
          .id = "filename")

# get data in good format
df_q <- df_q %>%
  mutate(date = ymd(Date),
         Q_m3s = Qls / 1000,
         site = case_when(str_detect(filename, "L140") ~ "civaux",
                          str_detect(filename, "K418") ~ "dampierre",
                          str_detect(filename, "K683") ~ "chinon")
  ) %>%
  select(site, date, Q_m3s) %>%
  mutate(Q_m3s = if_else(Q_m3s < 0, NA_real_, Q_m3s))

# Get data in wide format for plotting and modeling
df_qw <- pivot_wider(df_q, names_from = site, values_from = Q_m3s) 

# Only need data from 1990, add a month and year column 
df_qw <- filter.(df_qw, year(date) > 1989) %>%
  mutate(month = month(date),
         year = year(date))

# 2804 NAs for Chinon, 787 for Dampierre, 1 for Civaux
summary(df_qw)

# For some reason Dampierre is missing data, but I have some of it from before
df_qdam <- readRDS(file.path("data", "01a_discharge", 
                             "dampierre_discharge_daily")) %>%
  bind_rows(readRDS(file.path("data", "01a_discharge", 
                              "dampierre_prepared_1995_discharge"))) %>%
  bind_rows(readRDS(file.path("data", "01a_discharge", 
                              "dampierre_prepared_1996_discharge"))) %>%
  drop_na() %>%
  distinct() %>%
  mutate(discharge.daily = if_else(discharge.daily < 0, NA_real_, discharge.daily)) 

# Get that data in there
df_qw <- left_join(df_qw, df_qdam) %>%
  mutate(dampierre = if_else(is.na(dampierre), discharge.daily, dampierre))

# There's also some additional data from An for 2021/2022, get daily
df_q_2122 <- read_csv(file.path("data", "01_EDF", "hourly EDF data", 
                                "daily_discharge.csv")) %>%
  mutate(date = date(datetime)) %>%
  group_by(date) %>%
  summarize(chinon = mean(CHB, na.rm = T),
            dampierre = mean(DAM, na.rm = T)) %>%
  mutate(year = year(date),
         month = month(date))

# Add this to the end of the dataframe
df_qw <- filter.(df_q_2122, date > max(df_qw$date)) %>%
  bind_rows(df_qw) %>%
  arrange(date)

# Now dampierre only has 158 NAs
summary(df_qw)

# Look at relationships for gapfilling when data is missing
# Dampierre-Chinon is good
df_qw %>%
  ggplot(aes(x = dampierre, y = chinon)) +
  geom_point() + stat_smooth(method = "lm") +
  facet_wrap(~month(date)) +
  ggpubr::stat_regline_equation()

# Dampierre Civaux is less good, but usable
df_qw %>%
  ggplot(aes(x = civaux, y = dampierre)) +
  geom_point() + stat_smooth(method = "lm") +
  facet_wrap(~month(date)) +
  ggpubr::stat_regline_equation()

# Fill missing data for chinon
# Get models using robust iterated least squares
mod_dam_chi <- MASS::rlm(chinon ~ dampierre * month * year, data = df_qw)
summary(mod_dam_chi)

mod_dam_civ <- MASS::rlm(civaux ~ dampierre * month * year, data = df_qw)
summary(mod_dam_civ)

# Fill in missing NAs for dampierre using seasonal kalman filtering
dampierre_ts <- ts(df_qw$dampierre, deltat = 1/365)
dampierre_ts_clean <- imputeTS::na_seasplit(dampierre_ts, algorithm = "kalman")
df_qw$dampierre <- as.numeric(dampierre_ts_clean)

# Predictions for Chinon and Civaux from Dampierre to fill missing data
df_qw <- mutate(df_qw,
                pred_chi = as.numeric(predict(mod_dam_chi, newdata = df_qw)),
                pred_civ = as.numeric(predict(mod_dam_civ, newdata = df_qw)))

# Clean up a tad
df_qw <- mutate(df_qw, chinon = if_else(is.na(chinon), pred_chi, chinon),
                civaux = if_else(is.na(civaux), pred_civ, civaux)) %>%
  select(-pred_chi, -pred_civ, -discharge.daily)

# Get in long format
df_qclean <- df_qw %>%
  pivot_longer(cols = c(chinon, civaux, dampierre), names_to = "site",
               values_to = "Q_m3s")

# Read in temperature data ------------------------------------------------
# Hourly temperature for three Loire stations
df_tl <- read_csv(file.path("data", "01_EDF", "hourly EDF data",
                           "hourly_temp.csv"))

# Hourly temperature for Vienne
df_tv <- read_tsv(file.path("data", "01_EDF", "hourly EDF data", "Vienne",
                            "Temperature.txt")) %>%
  rename(datetime = Date, civaux = Value)

# Get all temperature data together
df_t <- df_tl %>%
  mutate(datetime = mdy_hm(datetime)) %>%
  rename(belleville = BEL, dampierre = DAM, chinon = CHB) %>%
  left_join(df_tv)

# 120198 missing hours of temperature for Civaux
summary(df_t)

# Fill missing data for vienne (civaux)
# Get models to fill in from other sites
mod_t<- MASS::rlm(civaux ~ dampierre + belleville + chinon, data = df_t)
summary(mod_t)

# Predictions for Civaux
df_t <- mutate(df_t,
                pred_civ = as.numeric(predict(mod_t, newdata = df_t)))

# Clean up a tad
df_t <- mutate(df_t, civaux = if_else(is.na(civaux), pred_civ, civaux)) %>%
  select(-pred_civ)

# Data all clean
summary(df_t)

# Long format
df_tclean <- df_t %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temp_C")
            
# Load radiation data -----------------------------------------------------
# Read Vienne radation from SAFRAN
df_rv <- read_csv(file.path("data", "01b_radiation", "radiation_vienne.csv"))

# Read in radiation from Orleans
df_ro <- readRDS(file.path("data", "01b_radiation", "light_dampierre_for_metab"))

# Clean the Orleans data (in umol/m2/s, want W/m2...weird past conversion)
df_ro <- df_ro %>%
  mutate(rad_Wm2 = if_else(is.na(light) & (hour(datetime) < 6 | hour(datetime) >= 21), 
                           0, light / 2.1)) %>%
  select(datetime, orleans = rad_Wm2)

# Fill in missing NAs for using seasonal kalman filtering
ro_ts <- ts(df_ro$orleans, deltat = 1/(365*24))
ro_ts_clean <- imputeTS::na_seasplit(ro_ts, algorithm = "kalman")
df_ro$rad_clean <- as.numeric(ro_ts_clean)

# Join this to vienne data
df_rad <- select(df_ro, datetime, orleans = rad_clean) %>%
  left_join(df_rv %>%
              mutate(datetime = mdy_hm(Date)) %>%
              select(datetime, vienne = SW))

# 80352 NAs for vienne
summary(df_rad)

# Fill missing data for vienne
# Get models to fill in from other orleans, 0 intercept
mod_r <- lm(vienne ~ 0 + orleans, data = df_rad)
summary(mod_r)

# Predictions for vienen
df_rad <- mutate(df_rad,
                 pred_civ = as.numeric(predict(mod_r, newdata = df_rad)))

# Clean up a tad
df_rad <- mutate(df_rad, vienne = if_else(is.na(vienne), pred_civ, vienne)) %>%
  select(-pred_civ)

# Get it in long format
df_radclean <- df_rad %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "rad_Wm2") %>%
  mutate(rad_Wm2 = if_else(rad_Wm2 < 0, 0, rad_Wm2))

# Combine all data --------------------------------------------------------
df_all <- df_tclean %>%
  mutate(date = date(datetime),
         qmatch = case_when(site == "belleville" ~ "dampierre",
                            TRUE ~ site),
         rmatch = case_when(site == "civaux" ~ "vienne",
                            TRUE ~ "orleans")) %>%
  left_join(select(df_qclean, qmatch = site, date, Q_m3s)) %>%
  left_join(select(df_radclean, rmatch = site, datetime, rad_Wm2))

# Still a lot of NAs for radiation because missing for recent years. fill with 
# streamMetabolizer
summary(df_all)

# For missing data in recent years, use estimated light from lat/long
df_all <- df_all %>%
  mutate(soltime = case_when(rmatch == "orleans" ~ 
                               streamMetabolizer::calc_solar_time(datetime, 2.6),
                             rmatch == "vienne" ~ 
                               streamMetabolizer::calc_solar_time(datetime, 0.2))) %>%
  mutate(rad_Wm2 = case_when(is.na(rad_Wm2) & rmatch == "orleans" ~ 
                               streamMetabolizer::calc_solar_insolation(soltime, 47.7),
                             is.na(rad_Wm2) & rmatch == "vienne" ~
                               streamMetabolizer::calc_solar_insolation(soltime, 47.2),
                             TRUE ~ rad_Wm2))

# Clean data and save
df_all <- df_all %>%
  select(-qmatch, -rmatch)
 
# Add the depth info ------------------------------------------------------
df_all <- df_all %>%
  mutate(depth_m = case_when(site == "belleville" ~ 0.3207 + sqrt(Q_m3s) * 0.052,
                             site == "dampierre" ~ 0.6171 + sqrt(Q_m3s) * 0.0716,
                             site == "chinon" ~ 0.4092 + sqrt(Q_m3s) * 0.0721,
                             site == "civaux" ~ -0.916 + log(Q_m3s) * 0.647)
  )

# save everything
saveRDS(df_all, file.path("data", "05_hourly_data_clean", "temp_discharge_rad_data.RDS"))  
