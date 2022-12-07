# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
library(lubridate)
library(sf)
library(tidyverse)


load("Data/Loire_DO/Loire_data_all.RData")

df <- tabPCAll %>%
  mutate(month = month(date_opecont),
         year = year(date_opecont),
         season = case_when(month %in% c(12,1,2) ~ "winter",
                            month %in% c(3,4,5) ~ "spring",
                            month %in% c(6,7,8) ~ "summer",
                            month %in% c(9,10,11) ~ "fall"))

ggplot(data = df,
       aes(x = year,
           y = pH,
           color = cd_site,
           group = cd_site)) +
  # stat_summary(alpha = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_grid(cols = vars(season)) +
  stat_smooth(method = "lm")
distinct(tabPC, cd_site)


df_spat <- df %>%
  left_join(sLoire %>%
              select(cd_site,
                     site,
                     lat = y,
                     long = x))

df_wint_sf <- filter(df_spat, season == "winter") %>%
  group_by(cd_site, year, lat, long) %>%
  summarize(meanpH = mean(pH, na.rm = T)) %>%
  ungroup() %>%
  group_by(cd_site, lat, long) %>%
  filter(n() > 10) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(meanpH ~ year, data = .)),
         tid = map(mod, broom::tidy)) %>%
  select(-data, -mod) %>%
  unnest(tid) %>%
  filter(term == "year", p.value < 0.05) %>%
  st_as_sf(coords = c("long", "lat"))
library(tmap)
tm_shape(df_wint_sf) +
  tm_symbols(col = "estimate")
