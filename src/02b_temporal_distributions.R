# -------------------------------------
# Author: Jacob Diamond
# Purpose: Look at temporal distrubtion of FCO2 as it relates to discharge
# Date: 19-03-2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
# Load daily NEP:CO2 archetypes
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS")) |>
  mutate(wyear = ifelse(month > 9, year + 1, year))

# flux by discharge ---------------------------------------------------
# Percent of fco2 occuring greater than mean
df_nepco2 |>
  # column for if the discharge is above or below mean
  mutate(more_less = if_else(Q > mean(Q, na.rm = T), "more", "less")) |>
  group_by(wyear, more_less) |>
  arrange(wyear, more_less) |>
  # Cumulative annual sum of FCO2 if Q is less/more than mean
  mutate(cusum_lessmean = cumsum(filtered_CO2_meanenh)) |>
  ungroup() |>
  group_by(wyear) |>
  mutate(totflux = sum(filtered_CO2_meanenh, na.rm = T)) |>
  group_by(wyear, more_less) |>
  # Get the total FCO2 if less or more than Q mean
  filter(cusum_lessmean == max(cusum_lessmean)) |>
  # Get the proportion of total annual flow for that break by year
  mutate(prop = cusum_lessmean/totflux) |>
  ungroup() |>
  group_by(more_less) |>
  # Get the summary stats overall
  summarize(mn = mean(prop),
            med = median(prop),
            sd = sd(prop))

# Proportion of total annual +FCO2 and Q (by water year)
df_prop <- mutate(df_nepco2, 
                  CO2 = if_else(filtered_CO2_meanenh < 0,
                                0,
                                filtered_CO2_meanenh)) |>
  group_by(wyear) |>
  arrange(wyear, date) |>
  drop_na(Q) |>
  mutate(Qsum = cumsum(Q),
         CO2sum = cumsum(CO2)) |>
  mutate(Qfrac = Qsum / max(Qsum),
         CO2frac = CO2sum / max(CO2sum))

# Take a look
ggplot(data = df_prop,
       aes(x = Qfrac,
           y = CO2frac,
           color = wyear,
           group = wyear)) +
  geom_line() +
  geom_abline() +
  theme_classic() +
  scale_color_viridis_c() +
  labs(x = "cumulative % of flow",
       y = "cumulative % of CO2") 

# Determine the date at which 80% of +FCO2 is emitted
df_80 <- df_prop |>
  group_by(wyear) |>
  mutate(jday = julian(date, origin = min(date))) |>
  slice(which.min(abs(CO2frac - 0.8)))
  
mean(df_80$jday) + ymd(20200901)
mean(df_80$Qfrac)
sd(df_80$Qfrac)

# Cumulative fluxes -------------------------------------------------------
# Quick look at cumulative NEP and CO2 over time
df_cumNEP <- df_nepco2 |>
  mutate(wyear = ifelse(month > 9, year + 1, year)) |>
  group_by(wyear) |>
  drop_na() |>
  mutate(dCO2 = (filtered_CO2_meanenh - filtered_CO2_2.5enh) / 1.96,
         dNEP = (filtered_NEP_mean - filtered_NEP_97.5) / 1.96) |>
  mutate(cumNEP = cumsum(filtered_NEP_mean),
         dcumNEP = sqrt(sum(dNEP^2)),
         cumCO2 = cumsum(filtered_CO2_meanenh),
         dcumCO2 = sqrt(sum(dCO2^2)),
         jday = julian(date, origin = min(date)))

# What day of year on average does 80% of the FCO2 come
df_cumNEP |>
  filter(cumCO2 / max(cumCO2) <= 0.8) |>
  slice_max(cumCO2, n = 1) |>
  ungroup() |>
  summarize(jd = mean(jday) + ymd(20200901))

# Data ends for plot text
data_ends <- df_cumNEP |>
  group_by(wyear) |>
  top_n(1, jday)
data_ends

# estimates of cumnep and cumco2
df_tot <- df_cumNEP |>
  summarize(maxNEP = max(cumNEP, na.rm = T),
            dNEP = mean(dNEP),
            maxCO2 = max(cumCO2, na.rm = T),
            dCO2 = mean(dCO2)) |>
  mutate(NEP2.5 = maxNEP - 1.96 * dNEP,
         NEP97.5 = maxNEP + 1.96 * dNEP,
         CO22.5 = maxCO2 - 1.96 * dCO2,
         CO297.5 = maxCO2 + 1.96 *dCO2)

ungroup(df_tot) |>
  # mutate(state = if_else(wyear < 2005, "planktonic", "benthic")) |>
  # group_by(state) |>
  summarize(nep = mean(maxNEP),
            nepsd = sd(maxNEP),
            co2 = mean(maxCO2),
            co2sd = sd(maxCO2))

# Plot of cumulative NEP
p_cumnep <- ggplot(data = df_cumNEP,
                   aes(color = wyear,
                       group = wyear)) +
  geom_line(aes(x = jday,
                y = cumNEP),
            size = 1.5) +
  geom_hline(yintercept = 0) + 
  ggrepel::geom_text_repel(aes(label = wyear, x = jday + 15, 
                               y = cumNEP), data = data_ends, size = 4) +
  # geom_ribbon(aes(x = jday, ymin = cumNEP - 1.96 * dcumNEP, 
  #                 ymax = cumNEP + 1.96 * dcumNEP,
  #                 fill = wyear), alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_c(name = 'water year') +
  scale_fill_viridis_c(name = 'water year') +
  theme(legend.position = "none") + #c(0.85, 0.2)) +
  labs(x=  "day in water year",
       y = expression("cumulative NEP flux to atmosphere (g"~C~m^{-2}~d^{-1}*")"))
p_cumnep
# Plot of cumulative NEP
p_cumco2 <- ggplot(data = df_cumNEP,
                   aes(color = wyear,
                       group = wyear)) +
  geom_line(aes(x = jday,
                y = cumCO2),
            size = 1.5) +
  geom_hline(yintercept = 0) + 
  ggrepel::geom_text_repel(aes(label = wyear, x = jday + 15, 
                               y = cumCO2), data = data_ends, size = 4) +
  # geom_ribbon(aes(x = jday, ymin = cumCO2 - 1.96 * dcumCO2, 
  #                 ymax = cumCO2 + 1.96 * dcumCO2,
  #                 fill = wyear), alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_c(name = 'water year') +
  scale_fill_viridis_c(name = 'water year') +
  theme_classic() +
  scale_color_viridis_c(name = 'water year') +
  theme(legend.position = c(0.1, 0.85)) +
  labs(x=  "day in water year",
       y = expression("cumulative"~CO[2]~"(g"~O[2]~m^{-2}~d^{-1}*")"))

p_cumco2


# Gini coefficient and lorenz curves --------------------------------------
# Lorenz curves 
(ggplot(data = mutate(df_nepco2, 
                      CO2 = if_else(filtered_CO2_meanenh <0,
                                    0,
                                    filtered_CO2_meanenh),
                      wyear = ifelse(month > 9, year + 1, year)) |>
          arrange(CO2),
        aes(x = CO2,
            color = wyear,
            group = wyear)) +
    gglorenz::stat_lorenz(desc = FALSE, size = 1.5) +
    geom_abline() +
    theme_classic() +
    scale_color_viridis_c() +
    labs(x = "cumulative % of time",
         y = "cumulative % of CO2")) |>
  ggplotly()

df_gini_Q <- mutate(df_nepco2, 
                    CO2 = if_else(filtered_CO2_meanenh <0,
                                  0,
                                  filtered_CO2_meanenh),
                    wyear = ifelse(month > 9, year + 1, year)) |>
  group_by(wyear) |>
  arrange(wyear, date) |>
  nest() |>
  mutate(gini = map(data, ~ineq::Gini(.x$Q))) |>
  select(-data) |>
  unnest(gini)

ineq::Lc()

mean(df_gini$gini)
mean(df_gini_Q$gini)
df_gini |>
  ungroup() |>
  mutate(state = if_else(wyear < 2006, "planktonic", "benthic")) |>
  group_by(state) |>
  summarize(mean = mean(gini),
            sd = sd(gini))
