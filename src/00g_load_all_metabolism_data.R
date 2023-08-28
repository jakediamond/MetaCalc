# 
# Purpose: Get all metabolism data
# Author: Jake Diamond
# Date: 15 November 2022
# 

# Load libraries
library(plotly)
library(imputeTS)
library(lubridate)
library(tidyverse)
# Load data and clean---------------------------------------------------------------
# Data directory
data_dir <- file.path("data", "02_metabolism", "newest_metabolism")

# File paths
files <- fs::dir_ls(data_dir, regexp = "GetFitDaily")

# Load all metabolism data
df_met <- files %>% 
  map_dfr(read_csv, .id = "source")

saveRDS(df_met, file.path("data", "02_metabolism", "update_metabolism_results_07dec2022.RDS"))
df_met <- readRDS(file.path("data", "02_metabolism", 
                            "update_metabolism_results_07dec2022.RDS"))

# get site and position info
df_met <- df_met %>%
  dplyr::mutate(source = str_extract(source, "([^/]+$)"),
                source = str_remove(source, ".csv")) %>%
  separate(col = source, into = c("type", "site", "pos", "period"), sep = "_")

# clean a bit
df_clean <- df_met %>%
  mutate(GPP = if_else(GPP_mean > 0, GPP_mean, NA_real_),
         ER = if_else(ER_mean < 0, ER_mean, NA_real_),
         K600 = K600_daily_mean) %>%
  select(site, pos, period, date, GPP, ER, K600) %>%
  mutate(site_f = factor(site),
         year = year(date),
         pos_f = factor(pos)) %>%
  mutate(site_f = fct_relevel(site_f, c("belleville", "dampierre", "civaux", "chinon")),
         pos_f = fct_relevel(pos_f, c("up", "down")))

# Save this
df_clean %>%
  select(-site_f, -pos_f, -year) %>%
  mutate(sitepos = paste(site, pos, sep = "_")) %>%
  group_split(sitepos) -> list_of_dfs

# name of each datafile
list_of_dfs %>%
  map(~pull(., sitepos)) %>%
  map(~unique(.)) -> names(list_of_dfs)

list_of_dfs %>%
  writexl::write_xlsx(file.path("data", "02_metabolism", 
                                "loire_metabolism.xlsx"))
# Little plot
ggplot(data = df_clean,
       aes(x = date,
           y = ER)) +
  # stat_summary() + 
  # stat_smooth(method = "lm") +
  geom_point() +
  theme_classic() +
  facet_grid(rows = vars(pos_f), cols = vars(site_f))

# ER vs K600
ggplot(data = df_clean,
       aes(x = K600,
           y = ER,
           color = year(date))) +
  scale_color_viridis_c() +
  geom_point() +
  theme_classic() +
  facet_grid(rows = vars(pos_f), cols = vars(site_f)) +
  scale_y_continuous(limits = c(-40, 0)) +
  scale_x_continuous(limits = c(0, 15))

