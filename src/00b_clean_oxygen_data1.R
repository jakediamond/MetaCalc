# 
# Purpose: To clean DO time-series for middle Loire and Vienne, reduce outliers
# Author: Jake Diamond
# Date: 25 October 2022
# 

# Load libraries
library(zoo)
library(signal)
library(imputeTS)
library(plotly)
library(lubridate)
library(tidyverse)
library(tidytable)

# Load functions
source(file.path("src", "000_DO_cleaning_functions.R"))

# Load all DO data---------------------------------------------------------------
# Data directory
df <- readRDS(file.path("data", "01_EDF", "do_timeseries.RDS"))
  
# Plot all sites upstream and downstream
p_DO_raw <- plot_ly(data = dplyr::filter(df, site != "stlaurent"),
                    x = ~datetime,
                    y = ~DO_mgL,
                    color = ~site,
                    linetype = ~pos,
                    colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
          #             title = TeX("\\text{pH}"),
           #            range = list(6, 11)),
         xaxis = list(title = ""),
         yaxis = list(title = "DO (mg/L)",
                      range = list(0, 30)),
         title = "raw DO data")

browsable(p_DO_raw)

# Nest data by site and period
# Split into 5 year periods
df_n <- df %>%
  distinct() %>%
  select(-QC) %>%
  mutate(yr = year(datetime) - min(year(datetime))) %>%
  group_by(site, pos) %>%
  nest() %>%
  mutate(period = map(data, ~as.numeric(cut_interval(.$yr, length = 5)))) %>%
  unnest(c(data, period)) %>%
  group_by(site, pos, period) %>%
  nest()

# Remove original data to keep workspace clean
rm(df)

# Apply functions ---------------------------------------------------------
# Clean the data and use the lowpass filter
df_n <- df_n %>%
  filter.(site != "stlaurent") %>% #not enough data at st laurent to care
  mutate(clean = map(data, clean_fun))

# Get all in one dataframe
df_DO <- unnest(df_n, clean)

# do a kalman imputation, max 12 hours
df_DO <- df_DO %>%
  group_by(site, pos) %>%
  arrange(datetime) %>%
  mutate(DO_use = na_kalman(DO, maxgap = 12))

# Plot all sites upstream and downstream
p_DO_clean <- plot_ly(data = dplyr::filter(df_DO, site == "belleville"), 
                    x = ~datetime,
                    y = ~DO_use,
                    color = ~pos) %>%
                    # colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)),
    title = "raw DO data")

htmltools::browsable(p_DO_clean)

# Save data
saveRDS(df_DO, file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))
