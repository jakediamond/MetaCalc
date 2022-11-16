# 
# Purpose: Get all metabolism data
# Author: Jake Diamond
# Date: 15 November 2022
# 

# Load libraries
library(plotly)
library(lubridate)
library(tidyverse)

# Load data and clean---------------------------------------------------------------
# Data directory
data_dir <- file.path("data", "02_metabolism")

# File paths
files <- fs::dir_ls(data_dir, regexp = "GetFitDaily")

# Load all metabolism data
df_met <- files %>% 
  map_dfr(read_csv, .id = "source")

saveRDS(df_met, file.path("data", "02_metabolism", "update_metabolism_results.RDS"))

# get site and position info
df_met <- df_met %>%
  dplyr::mutate(source = str_extract(source, "([^/]+$)"),
                source = str_remove(source, ".csv")) %>%
  separate(col = source, into = c("type", "site", "pos", "period"), sep = "_")

# clean a bit
df_met_clean <- df_met %>%
  mutate(GPP = if_else(GPP_daily_mean > 0, GPP_daily_mean, NA_real_),
         ER = if_else(ER_daily_mean < 0, ER_daily_mean, NA_real_),
         K600 = K600_daily_mean) %>%
  select(site, pos, period, date, GPP, ER, K600) %>%
  mutate(site_f = factor(site),
         year = year(date)) %>%
  mutate(site_f = fct_relevel(site_f, "belleville", "dampierre", "chinon", 
                             "civaux"))

# Little plot
ggplot(data = df_met_clean,
       aes(x = date,
           y = GPP)) +
  # stat_summary() + 
  # stat_smooth(method = "lm") +
  geom_point() +
  theme_classic() +
  facet_grid(rows = vars(pos), cols = vars(site_f))


a = filter(df_met_clean, site == "civaux", pos == "down")
gpp = imputeTS::na_kalman(a$K600)
gppts <- ts(gpp, deltat = 1/365)
d = decompose(gppts)
plot(d)
