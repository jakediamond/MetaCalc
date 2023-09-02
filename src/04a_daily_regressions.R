# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot exO2 vs exCO2 and exDIC
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
# library(plotly)
# library(htmltools)
library(patchwork)
library(tidyverse)
# library(tidytable)
source(file.path("src", "000_carbonate_functions.R"))

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data_final.RDS")) |>
  mutate(trophlux = str_replace(trophlux, "_", " ")) |>
  drop_na(trophlux)

ggplot(data = filter(df, date == ymd(19960616)),
       aes(x = DIC_uM,
           y = exO2)) +
  geom_point()

# hourly exO2 vs del DIC each day
df_mod_O2DIC <- df |>
  group_by(year, month, date) |>
  filter(n() > 20) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$exO2 ~ .$DIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

df_O2DIC <- ungroup(df_mod_O2DIC) |>
  hoist(gl, "r.squared") |>
  select(-gl, -mod, -data) |>
  tidyr::unnest(cols = tid) |>
  filter(term == ".$DIC_uM") |>
  select(date, slo_dic = estimate, se_dic = std.error, p_dic = p.value, r2_dic =r.squared)
df_O2DIC |>
  filter(p_dic < 0.01, r2_dic > 0.8) |>
  filter(between(slo_dic, -3, 1)) |>
# filter(df_O2DIC, p_dic < 0.01, r2_dic > 0.8) |>
  mutate(doy = yday(date),
         regime = ifelse(year(date) < 2012, "autotrophic", "benthic")) |>
  ggplot(aes(x = doy,
             y = slo_dic,
             color = regime)) +
  stat_smooth()


df_O2DIC |>
  filter(p_dic < 0.01, r2_dic > 0.8) |>
  filter(between(slo_dic, -3, 1)) |>
  mutate(doy = yday(date),
         regime = ifelse(year(date) < 2012, "autotrophic", "benthic")) |>
  group_by(regime, doy) |>
  summarize(slo = weighted.mean(slo_dic, (1/se_dic^2))) |>
  ggplot(aes(x = doy,
             y = slo,
             color = regime)) +
  stat_smooth()


# df_hco3 <- ungroup(df_mod_hco3) |>
#   hoist(gl, "r.squared") |>
#   select(-gl, -mod, -data) |>
#   tidyr::unnest(cols = tid) |>
#   filter(term == ".$HCO3_uM") |>
#   select(date, slo_hco3 = estimate, se_hco3 = std.error, p_hco3 = p.value, r2_hco3 =r.squared)
# 
# df_dic <- ungroup(df_mod) |>
#   hoist(gl, "r.squared") |>
#   select(-gl, -mod, -data) |>
#   tidyr::unnest(cols = tid) |>
#   filter(term == ".$DIC_uM") |>
#   select(date, slo_dic = estimate, se_dic = std.error, p_dic = p.value, r2_dic =r.squared)
# 
# 
# df_mods <- left_join(df_dic, df_hco3) |>
#   left_join(df_nepspc)

# saveRDS(df_mods, file.path("data", "daily_regressions_nep_dic_hco3.RDS"))



df_res_spc <- select(df_mod_spc, -data,- mod) |>
  unnest(tid)

# Plot of slope over year based on ecosystem state
p_slope_spc <- ggplot(data = df_res_spc |>
                     filter(term == ".$NEP") |>
                     mutate(doy = yday(date)) |>
                     mutate(state = if_else(year < 2012, "planktonic", "benthic")),
                   aes(x = doy,
                       y = estimate,
                       color = state)) +
  stat_smooth() + 
  theme_classic(base_size = 10) +
  ggsci::scale_color_aaas(name = "autotrophic state") +
  scale_x_continuous(breaks = seq(0, 365, 60)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  # annotate(geom = "text", x= 200, y = -2.06, label = "calcite precipitation",
  #          size = 3.5) +
  # annotate(geom = "text", x= 200, y = -0.97, label = expression(atop(CO[2]~"/"~HCO[3]^{`-`}, uptake)),  
  #          size = 3.5) +
  # geom_hline(yintercept = -1, linetype = "dashed") +
  # geom_hline(yintercept = -2) +
  labs(x = "day of year",
       y = expression(beta))
p_slope_spc

# Density plots
p_dens_hco3 <- ggplot(data = df_res_hco3 |>
                   filter(term == ".$HCO3_uM") |>
                   # filter(r.squared > 0.67,
                   #        p.value < 0.05) |>
                   mutate(doy = yday(date)) |>
                   mutate(state = if_else(year < 2012, "planktonic", "benthic")),
                 aes(x = estimate,
                     fill = state)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 10) +
  ggsci::scale_fill_aaas(name = "autotrophic regime") +
  scale_x_continuous(limits = c(-3, 2), breaks = seq(-3, 2, 0.5)) +
  theme(legend.position = c(0.12, 0.8)) +
  # facet_wrap(~year) +
  # annotate(geom = "text", x= 60, y = -1.92, label = "calcite precipitaiton") +
  # geom_vline(xintercept = -1, linetype = "dashed") +
  # geom_vline(xintercept = -2) +
  labs(x = expression(daily~exO[2]~-~DIC~slope~"("*beta*")"),
       y = "density")
p_dens_hco3
p_dens_alk
p_dens_ex
p_dens_night
p_dens_day
p_all_slope <- p_dens + annotation_custom(ggplotGrob(p_slopes),
                                       ymin = 0.35,
                                       ymax = 0.8,
                                       xmin = -0.4,
                                       xmax = 1.1)

p_all_slope

p <- p_all / p_all_slope + plot_annotation(tag_levels = "a")

ggsave(plot = p,
       filename = file.path("results", "figure3.png"),
       dpi = 300,
       units = "cm",
       height = 20,
       width = 18.4)

# See how much CO2, HCO3, and CaCO3 support -------------------------------
# First calculate slopes for exo2-exCO2
df_modco2 <- df |>
  group_by(year, month, date) |>
  drop_na(exDIC_uM, O2ex) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Summary data
df_resco2 <- select(df_modco2, -data,- mod) |>
  unnest(tid) 

# Get into tidy and long formats
df_dic <- select(df_mod, -mod, -data, -gl) |>
  unnest(tid) |>
  mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
  select(-(statistic:p.value)) |>
  pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
  left_join(select(df_mod, -mod, -data, -tid) |>
              hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
              select(-gl)) |>
  left_join(distinct(df, date, trophlux)) |>
  arrange(desc(r.squared), p.value, trophlux) |>
  select(-std.error_intercept)

df_co2 <- select(df_modco2, -mod, -data, -gl) |>
  unnest(tid) |>
  mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
  select(-(statistic:p.value)) |>
  pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
  left_join(select(df_mod, -mod, -data, -tid) |>
              hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
              select(-gl)) |>
  left_join(distinct(df, date, trophlux)) |>
  arrange(desc(r.squared), p.value, trophlux) |>
  select(-std.error_intercept)

# all model info together
df_mods <- df_dic |>
  rename(int_dic = estimate_intercept,
         slo_dic = estimate_slope,
         r2_dic = r.squared,
         p_dic = p.value,
         se_dic = std.error_slope) |>
  left_join(df_co2 |>
              rename(int_co2 = estimate_intercept,
                     slo_co2 = estimate_slope,
                     r2_co2 = r.squared,
                     p_co2 = p.value,
                     se_co2 = std.error_slope))

# Get mean daily LSI
df_lsi <- df |>
  group_by(date) |>
  summarize(lsi = mean(LSI, na.rm = T))

df_cont <- ungroup(df_mods) |>
  left_join(df_lsi) |>
  left_join(df_nepspc)
saveRDS(df_cont, file.path("data", "slopes_co2_dic_spc.R"))

x <- df_cont |>
  left_join(select(df, date, CO2ex) |>
              group_by(date) |>
              summarize(exCO2 = min(CO2ex))) |>
  filter(exCO2 > 0)
  # mutate(Z = (slo_co2 - slo_dic) / sqrt(se_co2^2 + se_dic^2),
  #        p = pnorm(Z)) |>
  # filter(p > 0.05,
  #        r2_dic > 0.7,
  #        between(slo_dic, -1.4, -0.3),
  #        between(slo_co2, -1.4, -0.3))
  # filter(exCO2 > 0,  between(slo_co2, -1.2, -0.6))

  # filter(!(trophlux == "heterotrophic source" & slo_co2 > -0.6)) |>
  # filter(slo_co2 + 1.96 * se_co2 > slo_dic - 1.96 * se_dic)
# Estimate contributions of DIC to GPP
df_cont2 <- df_cont |>
  # Just look at high quality data
  filter(r2_dic > 0.66) |>
  # Ignore winter days when GPP is too low
  filter(!(trophlux == "heterotrophic source" & slo_co2 > -0.6)) |>
  mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
  mutate(
    # determine what the source of DIC is
    source = case_when(
      # if the CO2 slope is between -1.2 and 0.6
      between(slo_co2, -1.2, -0.6) ~ "CO2",
      # If there is no distinguishable difference between DIC and CO2 slope
      slo_co2 + 1.96 * se_co2 > slo_dic - 1.96 * se_dic ~ "CO2",
      slo_co2 < -1.2 & between(slo_dic, -1.8, -0.8) ~ "HCO3",
      between(slo_dic,-5, -1.8)  & lsi > 0 ~ "CaCO3",
      between(slo_dic,-5, -1.8)  & lsi < 0 ~ "HCO3",
      TRUE ~ "CO2"
    )) 
# Summary of that data
df_cont |>
  group_by(regime, trophlux, source) |>
  summarize(count = n()) |>
  mutate(x = count  / sum(count)) |>
  select(-count) |>
  pivot_wider(names_from= source, values_from = x)


# Time series of slopes ---------------------------------------------------
# Plot of daily slopes
P_dic_time <- df_dic |>
  ggplot(aes(x = date,
             y = estimate_slope,
             # alpha = r.squared,
             color = trophlux)) + 
  geom_point() +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_y_continuous(limits = c(-3,3)) +
  scale_x_datetime(date_breaks = "5 years", date_labels = "%Y",
               name = "") +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = -1) +
  geom_hline(yintercept = -2, linetype = "dashed") +
  labs(y = expression(exO[2]*"–"*exDIC~slope)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(filename = file.path("results", "exo2_exdic_timeseries.png"),
       plot = p_dic_time,
       dpi = 300,
       width = 18.4,
       height = 12,
       units = "cm")


p_dic_examples <- df_dic |>
  group_by(trophlux, month) |>
  slice_sample(n = 1) |>
  left_join(df) |>
  ggplot(aes(x = exDIC_uM,
             y = O2ex,
             color = hr,
             group = interaction(date, trophlux))) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  facet_wrap(~month, scales = "free") +
  geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = expression("exDIC ("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))
p_dic_examples
ggsave(filename = file.path("results", "exo2_exdic_examples.png"),
       dpi = 300,
       width = 18.4,
       height = 12,
       units = "cm")

p_co2_examples <- df_co2 |>
  group_by(trophlux, month) |>
  slice_sample(n = 1) |>
  left_join(df) |>
  ggplot(aes(x = exCO2_uM,
             y = O2ex,
             color = hr,
             group = interaction(date, trophlux))) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  # scale_color_manual(name = "trophlux state",
  #                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  # scale_y_continuous(limits = c(-3,3)) +
  # scale_x_datetime(date_breaks = "5 years", date_labels = "%Y",
  #                  name = "") +
  theme_classic(base_size = 10) +
  facet_wrap(~month, scales = "free") +
  geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = expression(exCO[2]~"("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))
# theme(legend.position = "bottom",
#       legend.direction = "horizontal")
p_co2_examples
ggsave(filename = file.path("results", "exo2_exco2_examples.png"),
       dpi = 300,
       width = 18.4,
       height = 12,
       units = "cm")
# 
# 
# # Examples over time ------------------------------------------------------
# 
# 
# 
# 
# df_mods |>
#   select(trophlux, year, slo_dic, slo_co2) |>
#   pivot_longer(cols = contains("slo"), names_to = c("slope", "type"), names_sep = "_") |>
#   filter(between(value, -3, 3)) |>
#   ggplot(aes(x = value)) +
#   geom_density()+
#   facet_grid(trophlux~type)
# 
# ?wtd.var
# 
# a = dplyr::filter(df, date == ymd(20210325)) |>
#   filter(NEP_mmolO2m3 > 0 )
# summary(lm(O2ex~exCO2_uM, data = a))
# summary(MASS::rlm(O2ex~exCO2_uM, data = a))
# 
# 
# ggplot(data = dplyr::filter(df, date == ymd(20210325)),
#        aes(x = exCO2_uM,
#            y = O2ex,
#            color = NEP_mmolO2m3)) +
#   geom_path(size = 2) +
#   scale_color_viridis_c() +
#   geom_abline(slope = -1, intercept = 0) +
#   geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
#   theme_classic(base_size = 10) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0)
# 
# ggplot(data = x,
#        aes(x = trophlux,
#            y = slope,
#            fill = trophlux)) +
#   geom_violin() +
#   scale_y_continuous(limits = c(-3, 1))
# 
# ggplot(data = filter(y, r.squared > 0.9),
#        aes(x = trophlux,
#            y = slope,
#            fill = trophlux)) +
#   geom_violin() +
#   stat_summary() +
#   scale_y_continuous(limits = c(-3, 1))
# 
# 
# # AOU vs del DIC
# df_modco2_neppos <- df |>
#   group_by(year, month, date) |>
#   drop_na(exDIC_uM, O2ex) |>
#   filter(NEP_mmolO2m3 > 0) |>
#   # filter(light > 0) |>
#   nest() |>
#   mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
#          tid = map(mod, broom::tidy),
#          gl = map(mod, broom::glance))
# 
# df_res2 <- select(df_mod2, -data,- mod) |>
#   unnest(tid) 
# 
# # Plot of daily slopes
# df_res2|>
#   filter(term == ".$exCO2_uM") |>
#   ggplot(aes(x = date,
#              y = estimate)) + 
#   geom_point() +
#   # scale_y_continuous(limits = c(-3,3)) +
#   # scale_x_date(date_breaks = "1 year", date_labels = "%y",
#   #              name = "") +
#   stat_smooth() +
#   theme_bw() +
#   labs(y = "AOU vs exDIC slope")
# 
# 
# ggplot(data = df_res2_night |>
#          filter(term == ".$exCO2_uM") |>
#          mutate(doy = yday(date)) |>
#          mutate(state = if_else(year < 2012, "planktonic", "benthic")),
#        aes(x = doy,
#            y = -estimate,
#            color = state)) +
#   stat_smooth() + 
#   theme_bw()+
#   # scale_y_continuous(limits = c(0, 2.2)) +
#   geom_hline(yintercept = 1) +
#   geom_hline(yintercept = 2, linetype = "dashed")
# 
# y_night = select(df_mod2, -mod, -data, -gl) |>
#   unnest(tid) |>
#   mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
#   select(-(std.error:p.value)) |>
#   pivot_wider(names_from = term, values_from = estimate) |>
#   left_join(select(df_mod, -mod, -data, -tid) |>
#               hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
#               select(-gl)) |>
#   left_join(distinct(df, date, trophlux)) |>
#   arrange(desc(r.squared), p.value, trophlux)
# 
# seacarb::carb(flag = 21,
#               354,
#               6,
#               0,25)
# 
# ggplot(data = df_res |>
#          filter(term == ".$exDIC_uM") |>
#          mutate(doy = yday(date)),
#        aes(x = doy,
#            y = -estimate)) +
#   stat_smooth() + 
#   facet_wrap(~year) + 
#   theme_bw()+
#   scale_y_continuous(limits = c(0, 3)) +
#   geom_hline(yintercept = 1) +
#   geom_hline(yintercept = 2, linetype = "dashed")
# 
# 
# ggplot(data = df_res2 |>
#          filter(term == ".$exCO2_uM") |>
#          mutate(doy = yday(date)) |>
#          mutate(state = if_else(year < 2012, "planktonic", "benthic")),
#        aes(x = doy,
#            y = estimate,
#            linetype = state)) +
#   stat_smooth(color = "black") + 
#   theme_classic(base_size = 10) +
#   scale_linetype(name = "autotrophic state") +
#   scale_x_continuous(breaks = seq(0, 365, 30)) +
#   theme(legend.position = c(0.5, 0.8)) +
#   annotate(geom = "text", x= 60, y = -1.92, label = "calcite precipitaiton") +
#   geom_hline(yintercept = -1) +
#   geom_hline(yintercept = -2, linetype = "dashed") +
#   labs(x = "day of year",
#        y = expression(daily~exO[2]~-~exCO[2]~slope))
# 
# 
# 
# 
# z = df |>
#   filter(light > 0) |>
#   group_by(year, month, date, trophlux) |>
#   summarize(nepo = sum(NEP_mmolO2m3),
#             nepc = sum(NEP_mmolCm2hr),
#             dicrem = -sum(delDIC))
# 
# ggplot(data = z,
#        aes(y = nepc,
#            x = dicrem,
#            color = year)) +
#   geom_point() +
#   facet_wrap(~month) +
#   # geom_path(size = 2) +
#   scale_color_viridis_c() +
#   geom_abline(slope = 1, intercept = 0) +
#   # geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
#   theme_classic(base_size = 10) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   ggpubr::stat_regline_equation( aes(label =  paste(..eq.label.., 
#                                                     ..adj.rr.label.., 
#                                                     sep = "~~~~"))) +
#   labs(x = expression(DIC[removed]~"("*mmol~m^{-2}~d^{-1}*")"),
#        y = expression(NEP[C]~"("*mmol~m^{-2}~d^{-1}*")"))
# 
# 
# ggplot(data = z,
#        aes(y = nepo,
#            x = dicrem,
#            color = year)) +
#   geom_point() +
#   facet_wrap(~month) +
#   # geom_path(size = 2) +
#   scale_color_viridis_c() +
#   geom_abline(slope = 1, intercept = 0) +
#   # geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
#   theme_classic(base_size = 10) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   ggpubr::stat_regline_equation( aes(label =  paste(..eq.label.., 
#                                                     ..adj.rr.label.., 
#                                                     sep = "~~~~"))) +
#   labs(x = expression(DIC[removed]~"("*mmol~m^{-2}~d^{-1}*")"),
#        y = expression(NEP[C]~"("*mmol~m^{-2}~d^{-1}*")"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot() +
#   geom_path(data = dplyr::filter(df, date == ymd(20210325)),
#             aes(x = hr,
#                 y = HCO3,
#                 color = NEP_mmolO2m3),
#             size = 2) +
#   scale_color_viridis_c() +
#   stat_smooth(data = dplyr::filter(df, date == ymd(20210325),
#                                    between(hr, 12, 18)),
#             aes(x = hr,
#                 y = HCO3),
#             method = "lm") +
#   theme_classic(base_size = 10) +
#   ggpubr::stat_regline_equation(data = dplyr::filter(df, date == ymd(20210325),
#                                                      between(hr, 12, 18)),
#                                 aes(x = hr,
#                                     y = HCO3,
#                                     label =  paste(..eq.label.., 
#                                                    ..adj.rr.label.., 
#                                                    sep = "~~~~")))
# 
# # ALK and DIC
# t1 <- seacarb::carb(flag = 15, var1 = 0.00207, var2 = 0.00217,
#                     S = 0, T = 20, warn = "n")
# # CO2 and ALK
# t2 <- seacarb::carb(flag = 4, var1 = t1$CO2*0.001, var2 = 0.00207,
#                     S = 0, T = 20, warn = "n")
#  #HCO3 qnd ALK
# t3 <- seacarb::carb(flag = 11, var1 = t2$HCO3, var2 = 0.00207,
#                     S = 0, T = 20, warn = "n")
#   t3
#   
# 0.000104424*0.8

# Random examples ---------------------------------------------------------
df |>
  filter(date == ymd(20200703)) |>
         # between(hr, 5, 17)) |>
  ggplot(aes(x = DIC_uM,# - lag(AT)*1000,
             y = exO2,#- lag(exO2),
             color = hr)) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  ggpubr::stat_regline_equation() +
  # facet_wrap(~month, scales = "free") +
  # geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
  # geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  labs(y = expression(exO[2]~"("*mu*M*")"))

df |>
  filter(date == ymd(20200703)) |>
         # between(hr, 5, 12)) |>
  ggplot(aes(x = NEP, #- lag(HCO3),
             y = SpC - lag(SpC), #- lag(O2ex),
             color = hr)) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  ggpubr::stat_regline_equation()
  # facet_wrap(~month, scales = "free") +
  # geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
  # geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  # labs(y = expression(exO[2]~"("*mu*M*")"))


df |>
  filter(date == ymd_h(2003081000)) |>
  # between(hr, 5, 17)) |>
  ggplot(aes(x = DIC, #- lag(HCO3),
             y = O2ex, #- lag(O2ex),
             color = hr)) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  ggpubr::stat_regline_equation()
# facet_wrap(~month, scales = "free") +
# geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
# geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
# geom_hline(yintercept = 0) +
# geom_vline(xintercept = 0) +
# labs(y = expression(exO[2]~"("*mu*M*")"))

seacarb::carb(flag = 24, 
              var1 = 400,
              var2 = 1750/1E6,
              T = 25, S = 0)

seacarb::carb(flag = 24, 
              var1 = 400,
              var2 = 1850/1E6,
              T = 25, S = 0)




# old

# # Read in all the models
# df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))
# 
# ggplot(filter(df_mods, between(slo_spc, -0.2, 0)), 
#        aes(x = slo_spc,
#            y = slo_hco3)) +
#   stat_summary_bin()
# 
# # Check for uptake 
# df_test <- df_mods |>
#   left_join(distinct(df, date, year, trophlux, troph, sourcesink)) |>
#   # Just look at high quality data
#   filter(r2_dic > 0.66) |>
#   # Ignore winter days when GPP is too low
#   filter(!(trophlux == "heterotrophic source" & slo_dic > -0.6)) |>
#   mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
#   mutate(Z_hco3dic = abs((slo_hco3 - slo_dic)) / sqrt(se_hco3^2 + se_dic^2),
#          p_hco3dic = pt(Z_hco3dic, 46, lower.tail = FALSE),
#          Z_dic0.5 = abs((slo_dic - -0.5)) / sqrt(se_dic^2),
#          p_dic0.5 = pt(Z_dic0.5, 23, lower.tail = FALSE),
#          Z_hco30.5 = abs((slo_hco3 - -0.5)) / sqrt(se_hco3^2),
#          p_hco30.5 = pt(Z_hco30.5, 23, lower.tail = FALSE),
#          Z_spc0.13 = abs((slo_spc - -0.13)) / sqrt(se_spc^2),
#          p_spc0.13 = pt(Z_spc0.13, 23, lower.tail = FALSE),
#          Z_spc0.17 = abs((slo_spc - -0.17)) / sqrt(se_spc^2),
#          p_spc0.17 = pt(Z_spc0.17, 23, lower.tail = FALSE),
#          Z_spc0.04 = abs((slo_spc - -0.04)) / sqrt(se_spc^2),
#          p_spc0.04 = pt(Z_spc0.04, 23, lower.tail = FALSE)) |>
#   mutate(
#     # determine what the source of DIC is
#     source = case_when(
#       # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
#       p_hco3dic > 0.05  & p_dic0.5 < 0.05 ~ "HCO3",
#       p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 > 0.05 ~ "HCO3",
#       p_hco3dic < 0.05  & p_dic0.5 < 0.05 & p_spc0.04 < 0.05 & 
#         between(slo_hco3, -1.2, -0.8) & between(slo_spc, -0.06, -0.03)~ "HCO3",
#       p_dic0.5 > 0.05 & r2_dic > 0.7 ~ "CaCO3",
#       p_dic0.5 > 0.05 & r2_dic <= 0.7 & between(slo_spc, -0.2, -0.04) ~ "CaCO3",
#       between(slo_spc, -0.2, -0.06) & between(slo_hco3, -0.7, -0.4) ~ "CaCO3",
#       TRUE ~ "CO2"
#     )) 
# 
# x = df_test |>
#   filter(year > 1994) |>
#   filter(year > 1994, between(yday(date), 90 ,270)) |>
#   filter(source == "HCO3")
# # Summary of that data
# df_test |>
#   filter(year > 1994) |>
#   filter(year > 1994, between(yday(date), 90 ,270)) |>
#   group_by(troph , source) |>
#   # group_by(regime, source) |>
#   # group_by(troph, sourcesink, regime, source) |> #
#   summarize(count = n()) |>
#   mutate(x = count  / sum(count)) |>
#   select(-count) |>
#   pivot_wider(names_from= source, values_from = x) 
# # arrange(desc(regime), troph, desc(sourcesink))
# 
# 
# df_test |>
#   filter(year > 1994, between(yday(date), 90 ,270)) |>
#   # group_by(regime, source) |>
#   group_by(year, troph, sourcesink, regime, source) |> #
#   summarize(count = n()) |>
#   mutate(x = count  / sum(count)) |>
#   select(-count) |>
#   ungroup() %>%
#   group_by(regime, troph, sourcesink, source) |> #
#   summarize(mean = mean(x, na.rm = T) * 100,
#             sd = sd(x, na.rm = T) * 100) |>
#   mutate(x = paste(signif(mean,2), "+/-", signif(sd,2))) |>
#   select(-mean, -sd) |>
#   pivot_wider(names_from= source, values_from = x) |>
#   arrange(desc(regime), troph, desc(sourcesink))
# 
# 
# mean(df$CO2_uM, na.rm = T)
# mean(df$pCO2_cor_uatm, na.rm = T)
# 
# abs((-0.507 - -0.5)) / sqrt(0.023^2)
# 
# t.test(scale(1:23)*0.023 * sqrt(23) + -0.507, mu = -0.5)
# 
# t.test(scale(1:23)*0.01 + -0.05, scale(1:n2)*sd2 + mean2, ...)
# mutate(
#   # determine what the source of DIC is
#   source = case_when(
#     # if the DIC slope is not different than the hco3 slope, but not CaCO3 signal
#     p_hco3dic > 0.05  & p_dic0.5 < 0.05 ~ "HCO3",
#     # If there is no distinguishable difference between DIC and CO2 slope
#     between(slo_hco3, -1.2, -0.6) & between(slo_spc, -0.03, -0.05) ~ "HCO3",
#     p_spc0.04 > 0.05 & r2_spc > 0.7 ~ "HCO3",
#     p_spc0.17 > 0.05 & between(slo_hco3, -0.8, -0.1) ~ "HCO3",
#     between(slo_hco3, -0.6, 0) & between(slo_spc, -0.25, -0.05) ~ "CaCO3",
#     p_hco30.5 > 0.05 & p_spc0.13 > 0.05 ~ "CaCO3",
#     p_dic0.5 > 0.5 & r2_dic > 0.9 ~ "CaCO3",
#     p_hco30.5 > 0.05 & r2_hco3 > 0.9 ~ "CaCO3",
#     # between(slo_dic,-5, -1.8)  & lsi > 0 ~ "CaCO3",
#     # between(slo_dic,-5, -1.8)  & lsi < 0 ~ "HCO3",
#     TRUE ~ "CO2"
#   )) 