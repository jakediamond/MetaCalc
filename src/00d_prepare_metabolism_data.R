# -------------------------------------
# Author: Jake Diamond
# Purpose: Prepare the data for use in streamMetabolizer
# Date: 05 November 2022
# -------------------------------------

# Load libraries
library(lubridate)
library(tidyverse)

# Load physical data
df_phys <- readRDS(file.path("data", "05_hourly_data_clean", "temp_discharge_rad_data.RDS"))

# Load DO data
df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_final.RDS"))

# Combine data and clear workspace
df <- left_join(df_DO, df_phys) %>%
  ungroup()

rm(df_DO, df_phys)

# Calculate DO saturation
df_DOsat <- select(df, site, datetime, pos, wtr = temp_C) %>% 
  arrange(site, pos, datetime) %>%
  ungroup() %>%
  select(datetime, wtr) %>%
  as.data.frame() %>%
  LakeMetabolizer::o2.at.sat() %>%
  select(-datetime)

# save data for metabolism
savefile <- select(df, site, date, solar.time = soltime,
                   period, pos, DO.obs = DO, 
                   temp.water = temp_C, rad_Wm2, depth = depth_m,
                   discharge = Q_m3s) %>%
  ungroup() %>%
  mutate(light = streamMetabolizer::convert_SW_to_PAR(rad_Wm2),
         DO.sat = df_DOsat$do.sat) %>%
  select(-rad_Wm2, -discharge, -date) %>%
  distinct(site, pos, period, solar.time, .keep_all = TRUE)

# All data
saveRDS(savefile, file.path("data", "02_metabolism", "hourly_inputs.RDS"))
savefile <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS"))

# All discharge data
qfile <- select(df, site, pos, period, Q_m3s, date) %>%
  group_by(site, pos, period, date) %>%
  summarize(discharge.daily = mean(Q_m3s)) %>%
  ungroup()

# All discharge data
saveRDS(qfile, file.path("data", "02_metabolism", "daily_Q_inputs.RDS"))

# Split data
select(savefile, -discharge, -date) %>%
  mutate(filename = paste(site, pos, period, sep = "_")) %>%
  group_split(site, pos, period) -> list_of_dfs

# name of each datafile
list_of_dfs %>%
  map(~pull(., filename)) %>%
  map(~unique(.)) -> names(list_of_dfs)

imap(list_of_dfs, ~saveRDS(.x, file = file.path("data", "02_metabolism", 
                                                paste0(.y, "_for_meta", ".RDS"))))

# daily discharge, too
select(savefile, site, pos, period, discharge, date) %>%
  group_by(site, pos, period, date) %>%
  summarize(discharge.daily = mean(discharge)) %>%
  ungroup() %>%
  mutate(filename = paste(site, pos, period, sep = "_")) %>%
  group_split(site, pos, period) -> list_of_dfs_q

# name of each datafile
list_of_dfs_q %>%
  map(~pull(., filename)) %>%
  map(~unique(.)) -> names(list_of_dfs_q)

imap(list_of_dfs_q, ~saveRDS(.x, file = file.path("data", "02_metabolism", 
                                                  paste0(.y, "_daily_Q",".RDS"))))

# 
# 
# library(tidyverse)
# savefile <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS"))
# qfile <- readRDS(file.path("data", "02_metabolism", "daily_Q_inputs.RDS"))
# x = distinct(df_full, solar.time, DO.obs)
# 
# # Specific site and position
# df_full <- savefile[savefile$site == "dampierre" & savefile$pos == "up", ]
# df_q_full <- qfile[qfile$site == "dampierre" & qfile$pos == "up", ]
# 
# #n_period = unique(df_full$period)
# #n_site = unique(df_full$site)
# #n_pos = unique(df_full$pos)
# 
# #split the data into smaller groups
# #df_split <- split(df_full, ~cut(df_full$datetime, 5) + df_DO$site + df_DO$pos)
# #df_q_split <- split(df_q_full, ~cut(df_q_full$date, 5) + df_DO$site + df_DO$pos)
# 
# #unique periods
# n_periods <- unique(df_full$period)
# 
# df = df_full[df_full$period == n_periods[2],]
# df_q = df_q_full[df_q_full$period == n_periods[2],]
# 
# 
# x = filter(savefile, pos == "up", site == "dampierre", period == 1)
# y = distinct(x, solar.time) %>%
#   semi_join(x)
# x[duplicated(x$solar.time),]
# z = distinct(x, solar.time, .keep_all = TRUE)
# 
# filter(x, date(solar.time) == ymd(19941231))
# 
# savefile <- savefile %>%
#   distinct(site, pos, period, solar.time, .keep_all = TRUE)
