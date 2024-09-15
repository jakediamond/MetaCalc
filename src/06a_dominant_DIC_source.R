# -------------------------------------
# Author: Jake Diamond
# Purpose: to determine dominant DIC sources for autotrophs
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
# Hourly data
df <- readRDS(file.path("data", "hourly_data.RDS"))

# daily trophlux data
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS")) |>
  drop_na(troph) |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic"))

# Read in all the daily models
# df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))
df_mods <- readRDS(file.path("data", "daily_regressions.RDS"))


# Look at GPP vs diel SC --------------------------------------------------
df_sc <- df %>%
  group_by(date) %>%
  summarize(scamp = max(SpC) - min(SpC),
            SI = max(SI)) %>%
  left_join(distinct(df, date, GPP_mean, NEP_mean, Q_m3s))

left_join(df_sc, df_tf) %>%
  group_by(trophlux) %>%
  filter(Q_m3s < 250) %>%
  summarize(sc = quantile(scamp, na.rm = T))

ggplot(data = filter(df_sc, 
                     between(month(date), 5, 9),
                     scamp >0,
                     Q_m3s < 250),
       aes(x = GPP_mean,
           y = scamp,
           color = SI)) +
  geom_point() +
  ggpubr::stat_cor() +
  scale_color_viridis_c() +
  facet_wrap(~year(date))

# Get photoperiod
df_pp <- df %>%
  distinct(date) %>%
  mutate(jday = yday(date),
         photoperiod = meteor::photoperiod(x = jday, latitude = 47.6))

# Was CO2 depleted on a day? ----------------------------------------------
# determine if CO2 was depleted during the day
df_hr <- df |>
  # select(date, datetime, hr, NEP, CO2_uM, CO2eq_uM, depth_m, FCO2) %>%
  mutate(exCO2 = CO2_uM - CO2eq_uM) |>
  # estimate hourly GPP from NEP + ER / 24, and available CO2
  mutate(GPP_hr = if_else(NEP - ER_mean / 24 < 0, 0, NEP - ER_mean / 24),
         CO2_avail = CO2_uM + (-(ER_mean + FCO2_enh) / 24 * depth_m)) %>%
  group_by(date) %>%
  mutate(GPP_hr_mean = mean(GPP_hr * 24)) %>%
  ungroup() %>%
  mutate(ratio = if_else(is.na(GPP_mean / GPP_hr_mean), 1, GPP_mean / GPP_hr_mean))

df_dep <- df_hr %>%
  mutate(dep = if_else(
    # if the sun is shining
    between(hr, 6, 22) & 
      # and if at least half of GPP is greater than the CO2 available
      (GPP_hr / ratio) * depth_m  * 0.5 > CO2_avail,
      # and if at least half of NEP is greater than the CO2 available
      # NEP * depth_m  * 0.5 > CO2_uM,
                         # NEP * depth_m > CO2_uM,# - FCO2*depth_m,
                         # NEP > 0 &
                         # exCO2 <= 0,
                       "depleted",
                       "not")) |>
  group_by(year, date) |>
  summarize(depleted = sum(dep == "depleted"),
            SI = max(SI))
df_dep %>%
  ungroup() %>%
  summarize(count = sum(depleted, na.rm = T))

# Quick summary of depleted days
df_dep |>
  ungroup() |>
  left_join(df_pp) %>%
  # summarize(count = sum(depleted > photoperiod * 0.5, na.rm = T))
  summarize(count = sum(depleted > 6, na.rm = T))

# Plot some examples ------------------------------------------------------
# Quick look at depleted days
# Just 6 random examples
dep_samp <- ungroup(df) %>%
  semi_join(ungroup(df_dep) %>%
              group_by(date) %>%
              filter(sum(depleted) >6) %>%
              ungroup() %>%
              slice_sample(n = 6),
            by = "date")

p_co2 <- ggplot(data = dep_samp,
                aes(x = hr,
                    y = CO2_uM)) +
  geom_point(aes(color = NEP*depth_m)) +
  scale_color_viridis_c(name = "NEP") +
  theme_classic() +
  geom_line(aes(y = CO2eq_uM)) +
  facet_wrap(~date) #+
  # scale_y_continuous(limits = c(0, 500))
p_co2
ggsave(plot =p_co2,
       filename = file.path("results", "depleted_Examples.png"),
       dpi = 300,
       width = 18,
       height = 12,
       units = "cm")
p_dic <- ggplot(data = dep_samp,
                aes(x = hr,
                    y = DIC_uM - CO2_uM)) +
  geom_point(color = "red") +
  theme_classic() +
  facet_wrap(~date) +
  scale_y_continuous(position = "right") +#, limits = c(1500, 2000)) +
  theme(axis.text.y = element_text(color = "red"),
        axis.title.y = element_text(color = "red"),
        axis.ticks.y = element_line(color= "red"),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) +
  labs(x = "",
       y = expression(DIC~'('*mu*M*")"))

p <- p_co2 + p_dic + plot_layout(design = c(area(t = 1, l = 1, b = 5, r = 5),
                                            area(t = 1, l = 1, b = 5, r = 5)))
p
# Calculate statistics on slopes ------------------------------------------
# Check for uptake 
df_test <- df_mods |>
  mutate(
    # # Z score and p-value for if DIC slope is different than HCO3 (p > 0 is no difference)
    # Z_hco3dic = abs((slo_hco - slo_dic)) / sqrt(se_hco^2 + se_dic^2),
    # p_hco3dic = pt(Z_hco3dic, 46, lower.tail = FALSE),
    # Z score and p-value for if DIC slope is = 0.5 (p> 0 is no difference)
    Z_dic0.5 = abs((slo_dic - -0.5)) / sqrt(se_dic^2),
    p_dic0.5 = pt(Z_dic0.5, 23, lower.tail = FALSE),
    # # Z score and p-value for if HCO3 slope is = 0.5 (p> 0 is = 0.5)
    # Z_hco30.5 = abs((slo_hco - -0.5)) / sqrt(se_hco^2),
    # p_hco30.5 = pt(Z_hco30.5, 23, lower.tail = FALSE),
    # Z score and p-value for if NEP slope is = 0.13
    # Z_spc0.13 = abs((slo_nep - -0.13)) / sqrt(se_nep^2),
    # p_spc0.13 = pt(Z_spc0.13, 23, lower.tail = FALSE),
    # # Z score and p-value for if NEP slope is 0.17
    # Z_spc0.17 = abs((slo_nep - -0.17)) / sqrt(se_nep^2),
    # p_spc0.17 = pt(Z_spc0.17, 23, lower.tail = FALSE),
    # Z score and p-value for if NEP slope is < -0.04 (p < 0.01 is < -0.04)
    Z_spc0.09 = (slo_nep - -0.092) / se_nep,
    p_spc0.09 = pt(Z_spc0.09, 23, lower.tail = TRUE))

# Get all data together
df_up <- left_join(df_dep, df_test) |>
  left_join(df_tf) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      depleted > 6 & p_dic0.5 > 0.01 & r2_dic >= 0.9 ~ "CaCO3", # need to add SI here
      depleted > 6 & p_dic0.5 > 0.01 & r2_dic <= 0.9 & p_spc0.09 < 0.01 ~ "CaCO3",
      # depleted > 6 & SI > 0.4 & between(slo_nep, -0.2, -0.06) & between(slo_hco, -0.7, -0.4) ~ "CaCO3",
      # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
      # depleted > 6 & p_hco3dic > 0.01  & p_dic0.5 < 0.01 ~ "HCO3_a",
      # DIC slope is different than HCO3 slope, but NEP slope = -0.04
      depleted > 6 & p_dic0.5 < 0.01 ~ "HCO3",
      TRUE ~ "CO2"
    ))

# Look at how often slope is < -0.04
df_test %>%
  filter(year(date) > 1993, between(yday(date), 90 ,270)) |>
  summarize(n0.5 = sum(p_dic0.5 > 0.01),
            n0.09_0.5 = sum(p_spc0.09 < 0.01 & p_dic0.5 > 0.01))

df_up %>%
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  group_by(source) |>
  summarize(count = n())

# Summary of that data
df_up |>
  # data before 1994 is not reliable for SpC...no hourly data
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  group_by(source) |>
  # group_by(troph , source) |>
  # group_by(regime, source) |>
  # group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) 


# by years
df_y <- df_up |>
  filter(between(year, 1993, 2022), between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(year, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count))

# By year, then summarize
df_up |>
  filter(between(year, 1993, 2022), between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(year, troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  ungroup() %>%
  group_by(regime, troph, sourcesink, source) |> #
  summarize(mean = mean(x, na.rm = T) * 100,
           sd = sd(x, na.rm = T) * 100) |>
  mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))

 # By year, then summarize but only regime
df_y |>
  select(-count) |>
  ungroup() %>%
  group_by(regime, source) |> #
  summarize(mean = mean(x, na.rm = T) * 100,
            sd = sd(x, na.rm = T) * 100) |>
  mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime))

# Overall
df_up |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
  group_by(regime, source) |>
  # group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) #|>
  # arrange(desc(regime), troph, desc(sourcesink))

# Overall for flux states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(sourcesink, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(sourcesink))

# Overall for trophic states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(troph)

# Overall for trophlux states 
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, troph) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(troph, desc(sourcesink))

# Overall for trophlux states and regime
df_up |>
  filter(year > 1994, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # ungroup() %>%
  # group_by(regime, troph, sourcesink, source) |> #
  # summarize(mean = mean(x, na.rm = T) * 100,
  #           sd = sd(x, na.rm = T) * 100) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))


# Overall for trophlux states
df_up |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # group_by(regime, source) |>
  group_by(troph, sourcesink, regime, source) |> #
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  # ungroup() %>%
  # group_by(regime, troph, sourcesink, source) |> #
  # summarize(mean = mean(x, na.rm = T) * 100,
  #           sd = sd(x, na.rm = T) * 100) |>
  # mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
  # select(-mean, -sd) |>
  pivot_wider(names_from= source, values_from = x) |>
  arrange(desc(regime), troph, desc(sourcesink))

# Get mean pCO2 and Q for each year
df_qc <- df |>
  filter(year > 1993, between(yday(date), 90 ,270)) |>
  # mutate(wy = if_else(month > 8, year + 1, year)) |>
  group_by(year) |>
  summarize(qQ = median(Q_m3s),
            qC = median(pCO2_uatm))


mean(df$CO2_uM, na.rm = T)
mean(df$pCO2_cor_uatm, na.rm = T)

abs((-0.507 - -0.5)) / sqrt(0.023^2)

t.test(scale(1:23)*0.023 * sqrt(23) + -0.507, mu = -0.5)

t.test(scale(1:23)*0.01 + -0.05, scale(1:n2)*sd2 + mean2, ...)

left_join(df_dep, df_test) |>
  left_join(df_tf) |>
  filter(between(year, 1993, 2022), between(yday(date), 90 ,270)) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      slo_nep < -0.02 & r2_nep > 0.7 ~ "10%",
      # slo_nep < -0.042 & r2_nep > 0.7 ~ "25%", # need to add SI here
      # slo_nep < -0.097 & r2_nep > 0.7 ~ "50%",
      TRUE ~ "bad")) %>%
  group_by(source) %>%
  summarize(n = n())
