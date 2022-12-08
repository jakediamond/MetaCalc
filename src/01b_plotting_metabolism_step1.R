# -------------------------------------
# Author: Jake Diamond
# Purpose: Plot all metabolism
# Date: 16 November 2022
# -------------------------------------

# Load libraries
library(plotly)
library(lubridate)
library(tidyverse)

# Load all DO data---------------------------------------------------------------
# Load all the new gap data
# Data directory
data_dir <- file.path("data", "02_metabolism")

# Load the new metabolism data
df <- readRDS(file.path(data_dir, "update_metabolism_results.RDS"))

# files
files <- fs::dir_ls(data_dir, regexp = "GetFitDaily")

# Load gap fill metabolism data
df_met <- files %>% 
  map_dfr(read_csv, .id = "filename")

# Get site data
df_met_site <- df_met %>%
  mutate(source = str_extract(filename, "([^/]+$)")) %>%
  mutate(source = str_remove(source, ".csv")) %>%
  separate(source, into = c("type", "site", "pos", "period"), sep = "_")

# all data
df <- bind_rows(df, df_met_site)

# clean a bit
df_clean <- df_clean %>%
  mutate(GPP = if_else(GPP_mean < 0, NA_real_, GPP_mean),
         ER = if_else(ER_mean > 0, NA_real_, ER_mean),
         K600 = K600_daily_mean) %>%
  mutate(site_f = factor(site),
         pos_f = factor(pos)) %>%
  mutate(site_f = fct_relevel(site_f, c("belleville", "dampierre", "civaux", "chinon")),
         pos_f = fct_relevel(pos_f, c("up", "down")))

saveRDS(df_clean, file.path("data", "new_metabolism_results_all.RDS"))

df_clean <- readRDS(file.path("data", "02_metabolism", 
                              "update_metabolism_results.RDS"))

# Quick plot
ggplot(data = df_clean,
       aes(x = date,
           y = ER)) +
  geom_point() +
  theme_classic() +
  facet_grid(rows = vars(pos), cols = vars(site))

# ER vs K600
ggplot(data = df_clean,
       aes(x = K600,
           y = ER)) +
  geom_point() +
  theme_classic() +
  facet_grid(rows = vars(pos_f), cols = vars(site_f)) +
  scale_y_continuous(limits = c(-40, 0)) +
  scale_x_continuous(limits = c(0, 15))

# Get data into clean time series with no NA
df_ts <- df_clean %>%
  select(site, pos, date, GPP, ER, K600) %>%
  pivot_longer(cols = c(GPP, ER, K600)) %>%
  arrange(site, pos, name, date) %>%
  select(-date) %>%
  group_by(site, pos, name) %>%
  nest() %>%
  mutate(ts = map(data, ~ ts(.x, deltat = 1/365)),
         ts_clean = map(ts, imputeTS::na_seasplit))

# Decompose the data
df_ts <- df_ts %>%
  mutate(deco = map(ts_clean, decompose))

# Get it into usable format
df_trend <- df_ts %>%
  mutate(trend = map(deco, pluck, "trend"),
         # season = map(deco, pluck, "seasonal"),
         trend = map(trend, as.numeric)) %>%
         # season = map(season, as.numeric)) %>%
  select(site, pos, name, trend) %>%
  unnest(trend) %>%
  bind_cols(df_clean %>%
              select(site_f, site, pos, pos_f, date, GPP, ER, K600) %>%
              pivot_longer(cols = c(GPP, ER, K600)) %>%
              arrange(site, pos, name, date) %>%
              select(date, site_f, pos_f))

p_ts <- ggplot(data = filter(df_trend, name != "K600"),
       aes(x = date, color = site_f, group = interaction(site_f, name))) +
  # geom_line(aes(y = season)) +
  geom_line(aes(y = trend)) +
  theme_classic(base_size = 12) + 
  facet_grid(rows = vars(pos_f)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("orange", "blue", "black", "red"))
p_ts

ggsave(plot = p_trends,
       filename = file.path("results", "metabolism_trends.png"),
       dpi = 1200,
       units = "cm",
       height = 24,
       width = 28)


a=pluck(df_ts, 7, 8)
plot(a)
pluck(df_ts, 2, )
c = pluck(df_trend, 4, 1)
b =pluck(a, "seasonal")
a$trend
str(a)
plot(c)
