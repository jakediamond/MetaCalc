# -------------------------------------
# Author: Jake Diamond
# Purpose: Theoretical correspondence of SpC NEP curves
# Date: 2023-04-28
# -------------------------------------
library(patchwork)
library(tidyverse)

# All hourly data, but only want houlry O2 and pH (and daily discharge and depth)
df_hr <- readRDS(file.path("data", "hourly_data_final.RDS"))

# Cleaned conductivity time series
df_cond <- readRDS(file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))

# Look at relationship between NEP and C25 --------------------------------
# Get the good data for C25
df_nepspc <- df_hr %>%
  select(-SpC) %>%
  left_join(df_cond) %>%
  mutate(spc_dif = SpC - lag(SpC),
         doy = yday(date))

# Correlation between delta_cond and NEP
df_cor <- df_nepspc %>%
  group_by(year, month, date, doy) %>%
  drop_na(NEP_mmolO2m3, spc_dif) %>%
  filter(n() > 20) %>%
  nest() %>%
  mutate(#cor = map(data, ~cor(x = .$NEP_mmolO2m3, y = .$spc_dif,
          #                    method = "spearman")),
         lm = map(data, ~lm(.$spc_dif~.$NEP_mmolO2m3, data = .)),
         td = map(lm, broom::tidy))

# x <- mutate(df_cor,
#             p = map(lm, broom::glance))
df_nepspc <- df_cor %>%
  select(-cor, -lm, -data) %>%
  tidyr::unnest(cols = td) %>%
  filter(term == ".$NEP_mmolO2m3") %>%
  select(date, slo_spc = estimate, se_spc = std.error, p_spc = p.value)

# Plot monthly pattern of spearman correlation
# p_moncor <- df_cor %>%
#   unnest(cols = cor) %>%
#   ggplot(aes(x = month,
#              y = cor)) +
#   stat_summary() +
#   theme_classic(base_size = 10) +
#   scale_x_continuous(breaks = seq(1,12,1)) +
#   labs(x = "month",
#        y = expression(Delta*C[25]-NEP~"Spearman correlation, "*rho))
# p_moncor

# Mean slopes
p_slope <- df_cor %>%
  select(-data, -lm, -cor) %>%
  tidyr::unnest(cols =td) %>%
  filter(term == ".$NEP_mmolO2m3",
         p.value < 0.01) %>%
  ggplot(aes(x = month,
             y = estimate,
             group = month)) +
  geom_violin(scale = "count", draw_quantiles = 0.5, outlier.shape = NA) +
  stat_summary() +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -0.0445*0.92, color = "#88CCEE") +
  annotate(geom = "text", x = 1, y =- 0.032,
           color = "#88CCEE", label = expression(HCO[3]^{`-`}~uptake)) +
  geom_hline(yintercept = -0.132, color = "#44AA99") +
  annotate(geom = "text", x = 1, y = -0.12,
           color = "#44AA99", label = expression(CaCO[3]~ppt.)) +
  geom_hline(yintercept = -0.173, color = "#117733") +
  annotate(geom = "text", x = 1.2, y =- 0.19,
           color = "#117733", label = expression(atop(CaCO[3]~ppt.~and,
                                                     HCO[3]^{`-`}~uptake))) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  scale_y_continuous(lim = c(-0.3, 0.2)) +
  # coord_cartesian(ylim = c(-0.3, 0.2)) +
  labs(x = "month",
       y = expression(Delta*C[25]-NEP~"slope, "*beta[NEP]))
p_slope

# Example plots
a <- ungroup(df_cor) %>%
  select(-data, -lm) %>%
  tidyr::unnest(cols =td) %>%
  filter(term == ".$NEP_mmolO2m3",
         between(month,4,9),
         std.error < 0.1*abs(estimate)) %>%
  slice_min((abs(estimate - (-0.173))), n = 20) %>%
  left_join(df_hr) %>%
  group_by(date) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
str(a)
p_both <- mutate(df_nepspc, date = as.Date(date)) %>%
  semi_join(mutate(a, date = as.Date(date)), by = "date") %>%
  # filter(date %in% c(ymd(20160822), ymd(19950730), ymd(20030809))) %>%
  mutate(label = ifelse(hr == 3, as.character(date), NA_character_)) %>%
  ggplot(aes(x = NEP_mmolO2m3,
             y = spc_dif,
             color = exCO2_uM,
             group = date)) +
  geom_path(linewidth = 2) +
  scale_color_viridis_c(name = expression(exCO[2]~"("*mu*M*")")) +
  # ggpubr::stat_regline_equation(label.y.npc = 0.18) +
  # ggrepel::geom_text_repel(aes(label = label), nudge_x = 1.5, nudge_y = 1.2, 
  #                          size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.87, 0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = expression(NEP~"(mmol"~m^{-3}~h^{-1}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-1}~h^{-1}*")"))
p_both 

b <- ungroup(df_cor) %>%
  select(-data, -lm, -cor) %>%
  tidyr::unnest(cols =td) %>%
  filter(term == ".$NEP_mmolO2m3",
         between(month,4,9),
         std.error < 0.1*abs(estimate)) %>%
  slice_min((abs(estimate - (-0.0409))), n = 10) %>%
  left_join(df_nepspc) %>%
  group_by(date) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  filter(exCO2_uM < 30) %>%
  slice_head(n = 3)

p_hco3 <- mutate(df_nepspc, date = as.Date(date)) %>%
  dplyr::semi_join(b, by = "date") %>%
  mutate(label = ifelse(hr == 0, as.character(date), NA_character_)) %>%
  # filter(date %in% c(ymd(19940628), ymd(20000526), ymd(20100503))) %>%
  ggplot(aes(x = NEP_mmolO2m3,
             y = spc_dif,
             color = exCO2_uM,
             group = date)) +
  geom_path(linewidth = 2) +
  scale_color_viridis_c(name = expression(exCO[2]~"("*mu*M*")")) +
  ggpubr::stat_regline_equation(label.y.npc = 0.21) +
  ggrepel::geom_text_repel(aes(label = label), nudge_x = 1, nudge_y = 1, 
                           size = 4) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.87, 0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = expression(NEP~"(mmol"~m^{-3}~h^{-1}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-1}~h^{-1}*")"))
p_hco3

p <- (p_hco3 + p_both) / p_slope + plot_annotation(tag_levels = "a")
ggsave(plot = p,
       filename = file.path("results", "figure3.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 16)
