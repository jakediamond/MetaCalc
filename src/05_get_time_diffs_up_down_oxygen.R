# -------------------------------------
# Author: Jake Diamond
# Purpose: determine differences in hourly max DO between up and downstream
# Date: 15 November 2022
# -------------------------------------
# load libraries
library(lubridate)
library(plotly)
library(tidyverse)
# library(tidytable)

# Load all oxygen data -------------------------------------------------
df <- readRDS(file.path("data", "hourly_inputs_smaller_pds.RDS"))

# plot fev
plot_ly(data = filter.(df, site == "dampierre"), 
        x = ~solar.time,
        y = ~temp.water,
        color = ~pos,
        colors = c("blue", "orange")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "température (°C)",
                 range = list(0, 30)))

# hour of daily max
df_max <- df %>%
  mutate(date = date(solar.time),
         hour = hour(solar.time)) %>%
  group_by(site, pos, date) %>%
  filter.(DO.obs == max(DO.obs))
df_max %>%
  dplyr::group_by(site, date, pos) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 


# compare up and down
df_com <- ungroup(df_max) %>%
  select(site, pos, date, hour) %>%
  filter.(!(site == "chinon" & date == ymd(20101231))) %>%
  tidyr::pivot_wider(names_from = pos, values_from = hour) %>%
  mutate(dif = up - down)

# Summarize
df_sum <- df_com %>%
  group_by(site) %>%
  filter.(abs(dif) < 5) %>%
  drop_na()

df_sum %>%
  mutate(month = month(date),
         year = year(date)) %>%
  group_by(site, year, month) %>%
  summarize(difmean = mean(dif)) %>%
  mutate(time = ymd(paste(year, month, 1))) %>%
  ggplot(aes(x = time,
             y = difmean)) +
  geom_point() +
  facet_wrap(~site) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  stat_smooth() +
  labs(y = "difference in max DO timing (hrs)",
       x = "")
