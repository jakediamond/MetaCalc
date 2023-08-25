# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot examples of exO2 vs exCO2 and DIC and HCO3
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
# library(plotly)
# library(htmltools)
library(patchwork)
library(ggpmisc)
library(tidyverse)
# library(tidytable)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data_final.RDS")) |>
  mutate(trophlux = str_replace(trophlux, "_", " ")) |>
  drop_na(trophlux)

# Load regression model results
df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))

# Check for uptake 
df_test <- df_mods |>
  left_join(distinct(df, date, year, trophlux)) |>
  # Just look at high quality data
  filter(r2_dic > 0.66) |>
  # Ignore winter days when GPP is too low
  filter(!(trophlux == "heterotrophic source" & slo_dic > -0.6)) |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
  mutate(Z_hco3dic = abs((slo_hco3 - slo_dic)) / sqrt(se_hco3^2 + se_dic^2),
         p_hco3dic = pt(Z_hco3dic, 46, lower.tail = FALSE),
         Z_dic0.5 = abs((slo_dic - -0.5)) / sqrt(se_dic^2),
         p_dic0.5 = pt(Z_dic0.5, 23, lower.tail = FALSE),
         Z_hco30.5 = abs((slo_hco3 - -0.5)) / sqrt(se_hco3^2),
         p_hco30.5 = pt(Z_hco30.5, 23, lower.tail = FALSE),
         Z_spc0.13 = abs((slo_spc - -0.13)) / sqrt(se_spc^2),
         p_spc0.13 = pt(Z_spc0.13, 23, lower.tail = FALSE),
         Z_spc0.17 = abs((slo_spc - -0.17)) / sqrt(se_spc^2),
         p_spc0.17 = pt(Z_spc0.17, 23, lower.tail = FALSE),
         Z_spc0.04 = abs((slo_spc - -0.04)) / sqrt(se_spc^2),
         p_spc0.04 = pt(Z_spc0.04, 23, lower.tail = FALSE)) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
      p_hco3dic > 0.05  & p_dic0.5 < 0.05 ~ "HCO3",
      p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 > 0.05 ~ "HCO3",
      p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 < 0.05 & 
        between(slo_hco3, -1.2, -0.8) & between(slo_spc, -0.06, -0.03)~ "HCO3",
      p_dic0.5 > 0.05 & r2_dic > 0.7 ~ "CaCO3",
      p_dic0.5 > 0.05 & r2_dic <= 0.7 & between(slo_spc, -0.2, -0.04) ~ "CaCO3",
      between(slo_spc, -0.2, -0.06) & between(slo_hco3, -0.7, -0.4) ~ "CaCO3",
      TRUE ~ "CO2"
    )) 

df_example <- df_test |>
  filter(year > 2015, r2_dic > 0.8, r2_hco3 > 0.8, source != "CO2") |>
  group_by(source) |>
  slice_sample(n = 20) |>
  left_join(df)


x <- filter(df_example, source == "CaCO3")

p_dic_examples <- df_example |>
  ggplot(aes(x = DIC_uM,
             y = exO2,
             color = hr,
             group = interaction(date, source))) +
  geom_path(linewidth = 1.1) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  facet_wrap(~source) +
  # stat_poly_eq(use_label(c("eq"))) +
  geom_abline(slope = -0.5, intercept = 950, linetype = "dashed", linewidth = 1.2) +
  geom_abline(slope = -1, intercept = 1900, linewidth = 1.2) +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  labs(x = expression("DIC ("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))
p_dic_examples

p_hco3_examples <- df_example |>
  ggplot(aes(x = HCO3_uM,
             y = exO2,
             color = hr,
             group = interaction(date, source))) +
  geom_path(linewidth = 1.1) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  facet_wrap(~source) +
  # ggpubr::stat_regline_equation() +
  geom_abline(slope = -0.5, intercept = 850, linetype = "dashed", linewidth = 1.2) +
  geom_abline(slope = -1, intercept = 1700, linewidth = 1.2) +
  geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  labs(x = expression(HCO[3]^{`-`} ~"("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))
p_hco3_examples


p_spc_examples <- df_example |>
  ggplot(aes(x = NEP,
             y = SpC - lag(SpC),
             color = hr,
             group = interaction(date, source))) +
  geom_path(linewidth = 1.1) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  facet_wrap(~source, scales = "free") +
  scale_y_continuous(limits = c(-10, 10)) +
  # ggpubr::stat_regline_equation() +
  geom_abline(slope = -0.04, intercept = 0, linewidth = 1.2) +
  geom_abline(slope = -0.13, intercept = 0, linewidth = 1.2, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(y = expression(Delta*C[25]~"("*mu*S*~cm^{-2}~h^{-1}*")"),
       x = expression(NEP~"("*mmol~O[2]~m^{-3}~h^{-1}*")"))
p_spc_examples


p <- p_dic_examples / p_hco3_examples / p_spc_examples + plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect")
p

ggsave(plot = p,
       filename = file.path("results", "Figure_supp_uptake_examples.png"),
       dpi = 300,
       width = 18.4,
       height = 18,
       units = "cm")

# p_co2_examples <- df_co2 |>
#   group_by(trophlux, month) |>
#   slice_sample(n = 1) |>
#   left_join(df) |>
#   ggplot(aes(x = exCO2_uM,
#              y = O2ex,
#              color = hr,
#              group = interaction(date, trophlux))) +
#   geom_path(linewidth = 1.5) +
#   scale_color_viridis_c() +
#   # scale_color_manual(name = "trophlux state",
#   #                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
#   # scale_y_continuous(limits = c(-3,3)) +
#   # scale_x_datetime(date_breaks = "5 years", date_labels = "%Y",
#   #                  name = "") +
#   theme_classic(base_size = 10) +
#   facet_wrap(~month, scales = "free") +
#   geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
#   geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   labs(x = expression(exCO[2]~"("*mu*M*")"),
#        y = expression(exO[2]~"("*mu*M*")"))
# # theme(legend.position = "bottom",
# #       legend.direction = "horizontal")
# p_co2_examples
# ggsave(filename = file.path("results", "exo2_exco2_examples.png"),
#        dpi = 300,
#        width = 18.4,
#        height = 12,
#        units = "cm")
#