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

# Function to clean data --------------------------------------------------
lowpass_fun <- function(data, 
                        cutoff_frequency = 30) {
  # Re-interpolate all NAs so that there are none with kalman filtering
  data$value_an_int <- na_kalman(data$value)
  # Order the data, just in case
  data <- data[with(data, order(date)),]
  # Sampling rate [s^-1]
  sr <- 1 / (as.numeric(data$date[2] - data$date[1]) * 60 * 60 * 24)
  # Nyquist frequency = half the sampling rate
  nyq <- sr / 2
  # Cutoff frequency (days^-1)
  cutoff <- 1 / (cutoff_frequency * (60*60*24))
  # Normalized cutoff frequency for Butterworth filter
  W <- cutoff / nyq
  # Butterworth low-pass filter, digital, 2nd order
  myfilter <- signal::butter(2, W, type = 'low', plane = 'z')
  # Forward-reverse filter to remove phase-shift 
  # associated with Butterworth filter (must be in vector-form)
  vec <- as.vector(data$value_an_int)
  filtered <- signal::filtfilt(myfilter, vec)
  # Filtered data
  data$filtered <- filtered
  data <- data[with(data, order(date)), ]
  rem <- sr / cutoff
  data <- data[-c(1:rem, (nrow(data) - rem):nrow(data)),]
}

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

# name of each sheet will be the site
savefile %>%
  group_split(site, pos) -> names(list_of_dfs)

list_of_dfs %>%
  writexl::write_xlsx(path = "Headwaters/Data/DO_time_series.xlsx")

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

# Estimate of K600 amont et aval-----------------------------------------------
df_k <- select(df_clean, site_f, pos_f, date, K600) %>%
  mutate(pos_fr = fct_recode(pos_f, amont = "up", aval = "down"))

p_k <- ggplot(data = df_k,
              aes(x = site_f,
                  y = K600,
                  fill = pos_fr)) +
  stat_summary(fun=mean,
               position=position_dodge2(),geom="bar")+
  stat_summary(fun.data=mean_cl_normal,
               position=position_dodge2(),geom="errorbar") + 
  scale_fill_manual(values = c("#FFC107", "#1E88E5")) +
  theme_classic(base_size = 18) +
  labs(y = expression(K[600]~"("*d^{-1}*")"),
       x = "") +
  theme(axis.title.x = element_blank(), legend.position = "none")
p_k

ggsave(plot = p_k,
       filename = file.path("results", "K600_amont_aval.png"),
       dpi = 1200,
       units = "cm",
       height = 18,
       width = 20)

# Plot metabolism ---------------------------------------------------------
# Long format time series
df_l <- df_clean %>%
  select(site, pos, site_f, pos_f, date, GPP, ER, K600) %>%
  pivot_longer(cols = c(GPP, ER, K600)) %>%
  arrange(site, pos, name, date)

# Nest and lowpass filter
df_filt <- df_l %>%
  group_by(site, pos, site_f, pos_f, name) %>%
  nest() %>%
  mutate(filt = map(data, lowpass_fun))

# unnest
df_filt <- df_filt %>%
  select(-data) %>%
  tidyr::unnest(cols = filt) %>%
  mutate(pos_fr = fct_recode(pos_f, amont = "up", aval = "down"))

# metabolism plots
p_met <- ggplot(data = filter(df_filt, name != "K600"),
                   aes(x = date,
                       y = filtered,
                       color = name)) +
  geom_line() +
  theme_classic(base_size = 22) + 
  facet_grid(rows = vars(pos_fr), cols = vars(site_f)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("dark green", "black")) +
  labs(y = expression("metabolic flux (g "*O[2]~m^{-2}~d^{-1}*")"),
       x = "") +
  theme(axis.title.x = element_blank(), legend.position = "none")
p_met

ggsave(plot = p_met,
       filename = file.path("results", "metabolism_filtered.png"),
       dpi = 1200,
       units = "cm",
       height = 36,
       width = 56)


# Look at NEP -------------------------------------------------------------
# metabolism plots
df_nep <- filter(df_filt, name != "K600") %>%
  select(site_f, pos_fr, date, name, filtered) %>%
  tidyr::pivot_wider(names_from = name, values_from = filtered, values_fn = mean) %>%
  mutate(NEP = GPP + ER,
         NEPc = NEP * 32/12 * 0.85) %>%
  select(-GPP, -ER, -NEP) %>%
  pivot_longer(cols = NEPc) %>%
  group_by(site_f, pos_fr, name) %>%
  nest() %>%
  mutate(filtnep = map(data, lowpass_fun)) %>%
  select(-data) %>%
  tidyr::unnest(filtnep)

p_nep <- ggplot(data = df_nep,
                aes(x = date,
                    y = -filtered,
                    color = pos_fr)) +
  # geom_line() +
  stat_smooth(method = "loess", span = 1/10) +
  theme_classic(base_size = 22) + 
  facet_grid(cols = vars(site_f)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("#FFC107", "#1E88E5")) +
  labs(y = expression("NEP flux (g "*C~m^{-2}~d^{-1}*")"),
       x = "") +
  theme(axis.title.x = element_blank(), legend.position = "none")
p_nep

ggsave(plot = p_nep,
       filename = file.path("results", "nep_filtered.png"),
       dpi = 1200,
       units = "cm",
       height = 36,
       width = 56)

# Trend analysis ----------------------------------------------------------
# Get data into clean time series with no NA
df_ts <- ungroup(df_filt) %>%
  select(-date, -value, -value_an_int) %>%
  group_by(site, pos, name, pos_f, site_f, pos_fr) %>%
  nest() %>%
  mutate(ts = map(data, ~ ts(.x, deltat = 1/365)),
         ts_clean = map(ts, imputeTS::na_seasplit))

# x = pluck(df_ts, 6, 1)

# Decompose the data
df_ts <- df_ts %>%
  mutate(deco = map(ts_clean, decompose, "multiplicative")) #, filter = rep(1/5000, 5000))

# plot(pluck(df_ts, 10, 11))
# Get it into usable format
df_trend <- df_ts %>%
  mutate(trend = map(deco, pluck, "trend"),
         season = map(deco, pluck, "seasonal"),
         trend = map(trend, as.numeric),
         season = map(season, as.numeric)) %>%
  select(site, pos, name, trend, season) %>%
  unnest(c(trend, season)) %>%
  bind_cols(ungroup(df_filt) %>% select(date))

p_trends <- ggplot(data = filter(df_trend, name != "K600"),
                   aes(x = date, color = site_f, 
                       group = interaction(site_f, pos_f, name))) +
  # geom_line(aes(y = season)) +
  geom_line(aes(y = trend, linetype = name), size = 1.2) +
  theme_classic(base_size = 12) + 
  facet_grid(rows = vars(pos_fr)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("black", "#E69F00", "#009E73", "#0072B2")) +
  labs(y = "trend in metabolism") +
  theme(axis.title.x = element_blank())
p_trends

ggsave(plot = p_trends,
       filename = file.path("results", "metabolism_trends_v2.png"),
       dpi = 1200,
       units = "cm",
       height = 24,
       width = 28)
