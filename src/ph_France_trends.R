# -------------------------------------
# Author: Jacob Diamond
# Purpose: Look at pH trends for France
# Date: 30 November 2022
# -------------------------------------

# Load libraries
library(lubridate)
library(sf)
library(tmap)
library(tidyverse)
library(tidytable)

# Load the data
load(file.path("data", "00_water chemistry", "pH_data.RData"))

# Filter the data for at least 30 years of data
df <- thisdata %>%
  rename(site = cd_site) %>%
  mutate(year = year(date_opecont)) %>%
  group_by(site) %>%
  filter(n_distinct(year) > 30)

# Get mean winter values
df_win <- df %>%
  filter(month(date_opecont) %in% c(11,12,1,2,3)) %>%
  group_by(site, year) %>%
  summarize(pH = mean(resultat, na.rm = T))

# Get trends
df_sen <- ungroup(df_win) %>%
  dplyr::nest_by(site) %>%
  dplyr::mutate(model = list(trend::sens.slope(data$pH))) %>%
  dplyr::summarise(senslope = model["estimates"] %>% unlist(),
                   broom::tidy(model))

# Distribution of trends
p_hist_trend <- ggplot(data = df_sen,
                       aes(x = senslope)) +
  geom_histogram() +
  geom_vline(aes(xintercept = median(senslope)), color="red") +
  geom_vline(aes(xintercept = median(senslope)), color="red") +
  theme_classic() +
  labs(x = "winter pH change per year (Sen's slope)",
       title = "histogram of 1017 sites with at least 30 years of data")
p_hist_trend


# Plot of all sites in space
df_sen_sf <- df_sen %>%
  left_join(select(thisdata4v2, site = cd_site, x, y)) %>%
  st_as_sf(coords = c("x","y"), crs = "EPSG:2154") %>%
  dplyr::mutate(sens30 = 30 * senslope) %>%
  dplyr::filter(p.value < 0.05)

p_map <- tm_shape(df_sen_sf) +
  tm_bubbles(col = "sens30", style = "jenks", size = 0.5, legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE, title = "30-yr Sen's slope") 
p_map
