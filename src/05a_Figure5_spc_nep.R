# -------------------------------------
# Author: Jake Diamond
# Purpose: Theoretical correspondence of SpC NEP curves
# Date: 2023-04-28
# -------------------------------------
library(patchwork)
library(tidyverse)

# All hourly data, but  SpC is only houly after 1993
df <- readRDS(file.path("data", "hourly_data.RDS")) |>
  filter(year > 1992)

# Look at relationship between NEP and C25 --------------------------------
# Get the local delta for SpC and the julian day of year
df <- df |>
   mutate(spc_dif = SpC - lag(SpC),
         doy = yday(date))

# Linear regression between delta_cond and NEP
df_lm <- df |>
  group_by(year, month, date, doy) |>
  drop_na(NEP, spc_dif) |>
  # NEP is in mmol/m2/h, want in mmol/m3/h
  mutate(NEP = NEP / depth_smooth) |>
  # Want full days, at least 20/24 points
  filter(n() > 20) |>
  nest() |>
  mutate(lm = map(data, ~lm(.$spc_dif~.$NEP, data = .)),
         td = map(lm, broom::tidy),
         gl = map(lm, broom::glance))

# Summary of those daily models
df_nepspc <- df_lm |>
  hoist(gl, "r.squared") |>
  select(-gl, -lm, -data) |>
  unnest(cols = td) |>
  filter(term == ".$NEP") |>
  select(date, slo_spc = estimate, se_spc = std.error, p_spc = p.value, r2_spc =r.squared)

df_rect <- tibble(xmin = 0, xmax = 13,
                  ymin = c(-0.02, -0.042, -0.097),
                  ymax = c(-0.042, -0.097, -0.35),
                  group = factor(c("10-25 %", "25-50%", ">50%"),
                                 levels = c("10-25 %", "25-50%", ">50%")))

# Mean slopes across the year
p_slope <- #df_nepspc |>
  # filter(p_spc < 0.05, r2_spc > 0.66) |>
  ggplot() +
  geom_boxplot(data = df_nepspc, 
               aes(x = month,
                   y = slo_spc,
                   group = month),
               outlier.shape = NA) +
  stat_summary(data = df_nepspc, 
               aes(x = month,
                   y = slo_spc,
                   group = month)) +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_hline(yintercept = -0.097, color = "#117733", linewidth = 1.5) +
  geom_rect(data = df_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = group, fill = group),
    alpha = 0.4
  ) +
  scale_fill_viridis_d(name = "Calc. %") +
  # annotate(geom = "text", x = 6.5, y = -0.18, 
  #          color = "#117733", label = expression(possible~CaCO[3]~ppt.~dominance)) +
  # annotate(geom = "rect", xmin = 0, xmax = 13, ymin = -0.35, ymax = -0.097, 
  #          fill = "#117733", alpha = 0.4) +
  # geom_hline(yintercept = -0.0445*0.92, color = "#88CCEE") +
  # annotate(geom = "text", x = 1, y =- 0.032,
  #          color = "#88CCEE", label = expression(HCO[3]^{`-`}~uptake)) +
  # geom_hline(yintercept = -0.132, color = "#44AA99") +
  # annotate(geom = "text", x = 1, y = -0.12,
  #          color = "#44AA99", label = expression(CaCO[3]~ppt.)) +
  # geom_hline(yintercept = -0.173, color = "#117733") +
  # annotate(geom = "text", x = 1.2, y =- 0.19,
  #          color = "#117733", label = expression(atop(CaCO[3]~ppt.~and,
  #                                                    HCO[3]^{`-`}~uptake))) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  coord_cartesian(ylim = c(-0.3, 0.2), xlim = c(1, 12)) +
  theme(legend.position = c(0.55, 0.83)) +
  labs(x = "month",
       y = expression(Delta*C[25]-NEP~"slope, "*beta[NEP]))
p_slope

# Example plots
a <- ungroup(df_nepspc) |>
  filter(between(month,4,9),
         se_spc < 0.1*abs(slo_spc),
         r2_spc > 0.9) |>
  slice_min((abs(slo_spc - (-0.04))), n = 3) |>
  left_join(df) |>
  group_by(date) |>
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

p_25 <- mutate(df, date = as.Date(date)) |>
  semi_join(mutate(a, date = as.Date(date)), by = "date") |>
  # filter(date %in% c(ymd(20160822), ymd(19950730), ymd(20030809))) |>
  mutate(label = ifelse(hr == 3, as.character(date), NA_character_)) |>
  ggplot(aes(x = NEP / depth_smooth,
             y = spc_dif,
             color = exCO2,
             group = date)) +
  geom_path(linewidth = 2) +
  # scale_color_
  scale_color_viridis_c(name = expression(exCO[2]~"("*mu*M*")")) +
  ggpubr::stat_regline_equation(label.y.npc = 0.29) +
  ggrepel::geom_text_repel(aes(label = label), nudge_x = 1, nudge_y = 2,
                           size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.87, 0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = expression(NEP~"(mmol"~m^{-3}~h^{-1}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-1}~h^{-1}*")"))
p_25 

b <- ungroup(df_nepspc) |>
  filter(between(month,4,9),
         se_spc < 0.1*abs(slo_spc),
         r2_spc > 0.9) |>
  slice_min((abs(slo_spc - (-0.092))), n = 10) |>
  left_join(df) |>
  group_by(date) |>
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |>
  ungroup() |>
  filter(exCO2 < 30) |>
  slice_sample(n = 3)

p_50 <- mutate(df, date = as.Date(date)) |>
  dplyr::semi_join(b, by = "date") |>
  mutate(label = ifelse(hr == 0, as.character(date), NA_character_)) |>
  # filter(date %in% c(ymd(19940628), ymd(20000526), ymd(20100503))) |>
  # filter(date %in% c(ymd(20040629), ymd(20010607), ymd(20040613))) |>
  ggplot(aes(x = NEP / depth_smooth,
             y = spc_dif,
             color = exCO2,
             group = date)) +
  geom_path(linewidth = 2) +
  scale_color_viridis_c(name = expression(exCO[2]~"("*mu*M*")")) +
  ggpubr::stat_regline_equation(label.y.npc = 0.28) +
  ggrepel::geom_text_repel(aes(label = label), nudge_x = -0.5, nudge_y = 5, 
                           size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.87, 0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = expression(NEP~"(mmol"~m^{-3}~h^{-1}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-1}~h^{-1}*")"))
p_50

p <- (p_25 + p_50) / p_slope + plot_annotation(tag_levels = "a")
p
ggsave(plot = p,
       filename = file.path("results", "figure5_nepspc_v2.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 16)
