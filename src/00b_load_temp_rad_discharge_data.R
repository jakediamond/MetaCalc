# -------------------------------------
# Author: Jake Diamond
# Purpose: Get temperature, discharge, and radiation data together
# Date: 27 October 2022
# -------------------------------------
# load libraries
library(lubridate)
library(tidyverse)
library(tidytable)

# Load all discharge data -------------------------------------------------
# Data directory
data_dir <- file.path("data", "01a_discharge")

# files
files <- fs::dir_ls(data_dir, regexp = "\\.txt$")

# Load DO data
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
df_qw <- filter(df_qw, year(date) > 1989) %>%
  mutate(month = month(date),
         year = year(date))

summary(df_qw)

# For some reason Dampierre is missing data, but I have some of it from before
df_qdam <- readRDS(file.path("data", "01a_discharge", "dampierre_discharge_daily")) %>%
  bind_rows(readRDS(file.path("data", "01a_discharge", "dampierre_prepared_1995_discharge"))) %>%
  bind_rows(readRDS(file.path("data", "01a_discharge", "dampierre_prepared_1996_discharge"))) %>%
  drop_na() %>%
  distinct() %>%
  mutate(discharge.daily = if_else(discharge.daily < 0, NA_real_, discharge.daily)) 

# Get that data in there
df_qw <- left_join(df_qw, df_qdam) %>%
  mutate(dampierre = if_else(is.na(dampierre), discharge.daily, dampierre))

# Look at relationships
df_qw %>%
  ggplot(aes(x = dampierre, y = chinon)) +
  geom_point() + stat_smooth(method = "lm") +
  facet_wrap(~month(date)) +
  ggpubr::stat_regline_equation()

df_qw %>%
  ggplot(aes(x = civaux, y = dampierre)) +
  geom_point() + stat_smooth(method = "lm") +
  facet_wrap(~month(date)) +
  ggpubr::stat_regline_equation()

# Fill missing data for chinon
# Get models
mod_dam_chi <- MASS::rlm(chinon ~ dampierre * month * year, data = df_qw)
summary(mod_dam_chi)

# Predictions
df_qw <- mutate(df_qw,
                pred = as.numeric(predict(mod_dam_chi, newdata = df_qw)))

summary(df_qw)

# Clean up a tad
df_qw <- mutate(df_qw, chinon = if_else(is.na(chinon), pred, chinon)) %>%
  select(-discharge.daily, -pred)

# Fill in NAs and get in long format
df_qclean <- df_qw %>%
  # mutate(dampierre = ts(dampierre, frequency = 1))
  imputeTS::na_seasplit(algorithm = "kalman")




# Read in temperature data ------------------------------------------------
# Hourly temperature for three Loire stations
df_tl <- read_csv(file.path("data", "01_EDF", "hourly EDF data",
                           "hourly_temp.csv"))

# Hourly temperature for Vienne
df_tv <- read_tsv(file.path("data", "01_EDF", "hourly EDF data", "Vienne",
                            "Temperature.txt"))

# Get all temperature data together
df_t <- df_tl %>%
  mutate(datetime = mdy_hm(datetime)) %>%
  rename(belleville = BEL, dampierre = DAM, chinon = CHB) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temp_C") %>%
  bind_rows(rename(df_tv, datetime = Date, temp_C = Value) %>%
              mutate(site = "civaux"))

            
# Load radiation data -----------------------------------------------------
# Read Vienne radation from SAFRAN
df_rv <- read_csv(file.path("data", "01_EDF", "hourly EDF data", "Vienne",
                            "Radiation.csv"))

# Read in radiation from Orleans
df_ro <- readRDS(file.path("data", "light_dampierre_for_metab"))

# Combine this data (in umol/m2/s, want W/m2...weird past conversion)
df_rad <- df_ro %>%
  mutate(rad_Wm2 = if_else(is.na(light) & (hour(datetime) < 6 | hour(datetime) > 21), 
                           0, light / 2.1),
         site = "orleans") %>%
  select(site, datetime, rad_Wm2) %>%
  bind_rows(df_rv %>%
              mutate(datetime = mdy_hm(Date),
                     site = "vienne") %>%
              select(site, datetime, rad_Wm2 = SW))


# Combine all data --------------------------------------------------------
df_all <- df_t %>%
  mutate(date = date(datetime),
         qmatch = case_when(site == "belleville" ~ "dampierre",
                            TRUE ~ site),
         rmatch = case_when(site == "civaux" ~ "vienne",
                            TRUE ~ "orleans")) %>%
  left_join(select(df_q, qmatch = site, date, Q_m3s)) %>%
  left_join(select(df_rad, rmatch = site, datetime, rad_Wm2))

# Clean data and save
df_all <- df_all %>%
  select(-qmatch, -rmatch)

saveRDS(df_all, file.path("data", "temp_discharge_rad_data.RDS"))  
