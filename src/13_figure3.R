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

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data_final.RDS")) |>
  mutate(trophlux = str_replace(trophlux, "_", " ")) |>
  drop_na(trophlux) |>
  mutate(flux_update = if_else(FCO2 > 0, "source", "sink"),
         trophic_update = if_else(NEP > 0, "autotrophic", "heterotrophic"),
         trophlux_update = paste(trophic_update, flux_update)) |>
  drop_na(trophic_update, flux_update)

# Load regression model results
df_mods <- readRDS(file.path("data", "daily_regressions_nep_dic_hco3.RDS"))

# How many samples were under vs oversaturated
sum(df$exO2 > 0, na.rm = T) / nrow(df)
sum(df$exCO2 > 0, na.rm = T) / nrow(df)

# Overall plot of exDIC vs ex o2
p_all_dic <- ggplot(data = df,
       aes(x = exDIC,
           y = exO2,
           color = trophlux_update)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.16, 0.175),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill = "transparent", color = "transparent")) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  labs(x = expression("exDIC ("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))

# p_all_dic

# Overall plot of exCO2 vs ex o2
p_all_co2 <- ggplot(data = df,
                    aes(x = exCO2,
                        y = exO2,
                        color = trophlux)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  labs(x = expression(exCO[2]~"("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))

# p_all_co2

p_all <- p_all_dic + annotation_custom(ggplotGrob(p_all_co2),
                                                  ymin = 40,
                                                  ymax = 460,
                                                  xmin = 20,
                                                  xmax = 400)


# Photosynthetic quotient -------------------------------------------------
# Get just the good data for exO2 vs DIC relationship
df_pq <- df_mods |>
  select(date, contains("dic")) |>
  filter(p_dic < 0.05, r2_dic > 0.66) |>
  filter(between(slo_dic, -3, 1)) |>
  mutate(doy = yday(date)) |>
  mutate(regime = if_else(year(date) < 2012, "planktonic", "benthic"))

# Check for bimodality
multimode::modetest(df_pq$slo_dic)

# Plot of slope (=PQ) over year based on ecosystem regime
p_pq_t <- ggplot(data = df_pq,
                   aes(x = doy,
                       y = slo_dic,
                       color = regime)) +
  stat_smooth() + 
  theme_classic(base_size = 10) +
  ggsci::scale_color_aaas(name = "autotrophic regime") +
  scale_x_continuous(breaks = seq(0, 365, 60)) +
  scale_y_continuous(breaks = seq(-1.2, 0, 0.2)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  # annotate(geom = "text", x= 200, y = -2.06, label = "calcite precipitation",
  #          size = 3.5) +
  # annotate(geom = "text", x= 200, y = -0.97, label = expression(atop(CO[2]~"/"~HCO[3]^{`-`}, uptake)),  
  #          size = 3.5) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  # geom_hline(yintercept = -2) +
  labs(x = "day of year",
       y = expression(PQ~"("*beta[DIC]*")"))
p_pq_t


# dataframe for medians
df_text <- df_pq |>
  filter(between(doy, 90, 270)) |>
  group_by(regime) |>
  summarize(pq = median(slo_dic, na.rm = T)) |>
  mutate(x = c(0, -2), y = 1.4,
         label = paste("median PQ =", signif(pq,2)))

# Density plots
p_pq_dens <- ggplot(data = filter(df_pq, between(doy, 90, 270)),
                 aes(x = slo_dic,
                     group = regime)) +
  geom_density(aes(fill = regime), alpha = 0.5) +
  stat_summary(aes(xintercept = after_stat(x), y = 0, color = regime), fun = median,
                  linetype = "dashed", geom = "vline", linewidth = 1, orientation = "y") +
  geom_text(data = df_text, aes(x = x, y = y, label = label, color = regime)) +
  theme_classic(base_size = 10) +
  guides(text = "none", color = "none") +
  ggsci::scale_fill_aaas(name = "autotrophic regime") +
  ggsci::scale_color_aaas(name = "autotrophic regime") +
  theme(legend.position = c(0.2, 0.5),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = expression(daily~exO[2]~-~DIC~slope~"("*beta[DIC]*")"),
       y = "density")
p_pq_dens
p_all_slope <- p_pq_dens + annotation_custom(ggplotGrob(p_pq_t),
                                       ymin = 0.1,
                                       ymax = 1.2,
                                       xmin = -0.7,
                                       xmax = 1)

# p_all_slope

p <- p_all / p_all_slope + plot_annotation(tag_levels = "a")
# p
ggsave(plot = p,
       filename = file.path("results", "figure3_smaller.png"),
       dpi = 300,
       units = "cm",
       height = 15,
       width = 13.5)


# 
# # See how much CO2, HCO3, and CaCO3 support -------------------------------
# # First calculate slopes for exo2-exCO2
# df_modco2 <- df |>
#   group_by(year, month, date) |>
#   drop_na(exDIC_uM, O2ex) |>
#   nest() |>
#   mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
#          tid = map(mod, broom::tidy),
#          gl = map(mod, broom::glance))
# 
# # Summary data
# df_resco2 <- select(df_modco2, -data,- mod) |>
#   unnest(tid) 
# 
# # Get into tidy and long formats
# df_dic <- select(df_mod, -mod, -data, -gl) |>
#   unnest(tid) |>
#   mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
#   select(-(statistic:p.value)) |>
#   pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
#   left_join(select(df_mod, -mod, -data, -tid) |>
#               hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
#               select(-gl)) |>
#   left_join(distinct(df, date, trophlux)) |>
#   arrange(desc(r.squared), p.value, trophlux) |>
#   select(-std.error_intercept)
# 
# df_co2 <- select(df_modco2, -mod, -data, -gl) |>
#   unnest(tid) |>
#   mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
#   select(-(statistic:p.value)) |>
#   pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
#   left_join(select(df_mod, -mod, -data, -tid) |>
#               hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
#               select(-gl)) |>
#   left_join(distinct(df, date, trophlux)) |>
#   arrange(desc(r.squared), p.value, trophlux) |>
#   select(-std.error_intercept)
# 
# # all model info together
# df_mods <- df_dic |>
#   rename(int_dic = estimate_intercept,
#          slo_dic = estimate_slope,
#          r2_dic = r.squared,
#          p_dic = p.value,
#          se_dic = std.error_slope) |>
#   left_join(df_co2 |>
#               rename(int_co2 = estimate_intercept,
#                      slo_co2 = estimate_slope,
#                      r2_co2 = r.squared,
#                      p_co2 = p.value,
#                      se_co2 = std.error_slope))
# 
# # Get mean daily LSI
# df_lsi <- df |>
#   group_by(date) |>
#   summarize(lsi = mean(LSI, na.rm = T))
# 
# df_cont <- ungroup(df_mods) |>
#   left_join(df_lsi) |>
#   left_join(df_nepspc)
# saveRDS(df_cont, file.path("data", "slopes_co2_dic_spc.R"))
# 
# x <- df_cont |>
#   left_join(select(df, date, CO2ex) |>
#               group_by(date) |>
#               summarize(exCO2 = min(CO2ex))) |>
#   filter(exCO2 > 0)
#   # mutate(Z = (slo_co2 - slo_dic) / sqrt(se_co2^2 + se_dic^2),
#   #        p = pnorm(Z)) |>
#   # filter(p > 0.05,
#   #        r2_dic > 0.7,
#   #        between(slo_dic, -1.4, -0.3),
#   #        between(slo_co2, -1.4, -0.3))
#   # filter(exCO2 > 0,  between(slo_co2, -1.2, -0.6))
# 
#   # filter(!(trophlux == "heterotrophic source" & slo_co2 > -0.6)) |>
#   # filter(slo_co2 + 1.96 * se_co2 > slo_dic - 1.96 * se_dic)
# # Estimate contributions of DIC to GPP
# df_cont2 <- df_cont |>
#   # Just look at high quality data
#   filter(r2_dic > 0.66) |>
#   # Ignore winter days when GPP is too low
#   filter(!(trophlux == "heterotrophic source" & slo_co2 > -0.6)) |>
#   mutate(regime = if_else(year < 2005, "planktonic", "benthic")) |>
#   mutate(
#     # determine what the source of DIC is
#     source = case_when(
#       # if the CO2 slope is between -1.2 and 0.6
#       between(slo_co2, -1.2, -0.6) ~ "CO2",
#       # If there is no distinguishable difference between DIC and CO2 slope
#       slo_co2 + 1.96 * se_co2 > slo_dic - 1.96 * se_dic ~ "CO2",
#       slo_co2 < -1.2 & between(slo_dic, -1.8, -0.8) ~ "HCO3",
#       between(slo_dic,-5, -1.8)  & lsi > 0 ~ "CaCO3",
#       between(slo_dic,-5, -1.8)  & lsi < 0 ~ "HCO3",
#       TRUE ~ "CO2"
#     )) 
# # Summary of that data
# df_cont |>
#   group_by(regime, trophlux, source) |>
#   summarize(count = n()) |>
#   mutate(x = count  / sum(count)) |>
#   select(-count) |>
#   pivot_wider(names_from= source, values_from = x)
# 
# 
# # Time series of slopes ---------------------------------------------------
# # Plot of daily slopes
# P_dic_time <- df_dic |>
#   ggplot(aes(x = date,
#              y = estimate_slope,
#              # alpha = r.squared,
#              color = trophlux)) + 
#   geom_point() +
#   scale_color_manual(name = "trophlux state",
#                      values = c("black", "#E69F00", "#0072B2", "#009E73")) +
#   scale_y_continuous(limits = c(-3,3)) +
#   scale_x_datetime(date_breaks = "5 years", date_labels = "%Y",
#                name = "") +
#   theme_classic(base_size = 10) +
#   geom_hline(yintercept = -1) +
#   geom_hline(yintercept = -2, linetype = "dashed") +
#   labs(y = expression(exO[2]*"–"*exDIC~slope)) +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal")
# 
# ggsave(filename = file.path("results", "exo2_exdic_timeseries.png"),
#        plot = p_dic_time,
#        dpi = 300,
#        width = 18.4,
#        height = 12,
#        units = "cm")
# 
# 
# p_dic_examples <- df_dic |>
#   group_by(trophlux, month) |>
#   slice_sample(n = 1) |>
#   left_join(df) |>
#   ggplot(aes(x = exDIC_uM,
#              y = O2ex,
#              color = hr,
#              group = interaction(date, trophlux))) + 
#   geom_path(linewidth = 1.5) +
#   scale_color_viridis_c() +
#   theme_classic(base_size = 10) +
#   facet_wrap(~month, scales = "free") +
#   geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
#   geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   labs(x = expression("exDIC ("*mu*M*")"),
#        y = expression(exO[2]~"("*mu*M*")"))
# p_dic_examples
# ggsave(filename = file.path("results", "exo2_exdic_examples.png"),
#        dpi = 300,
#        width = 18.4,
#        height = 12,
#        units = "cm")
# 
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
# # 
# # 
# # # Examples over time ------------------------------------------------------
# # 
# # 
# # 
# # 
# # df_mods |>
# #   select(trophlux, year, slo_dic, slo_co2) |>
# #   pivot_longer(cols = contains("slo"), names_to = c("slope", "type"), names_sep = "_") |>
# #   filter(between(value, -3, 3)) |>
# #   ggplot(aes(x = value)) +
# #   geom_density()+
# #   facet_grid(trophlux~type)
# # 
# # ?wtd.var
# # 
# # a = dplyr::filter(df, date == ymd(20210325)) |>
# #   filter(NEP_mmolO2m3 > 0 )
# # summary(lm(O2ex~exCO2_uM, data = a))
# # summary(MASS::rlm(O2ex~exCO2_uM, data = a))
# # 
# # 
# # ggplot(data = dplyr::filter(df, date == ymd(20210325)),
# #        aes(x = exCO2_uM,
# #            y = O2ex,
# #            color = NEP_mmolO2m3)) +
# #   geom_path(size = 2) +
# #   scale_color_viridis_c() +
# #   geom_abline(slope = -1, intercept = 0) +
# #   geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
# #   theme_classic(base_size = 10) +
# #   geom_hline(yintercept = 0) +
# #   geom_vline(xintercept = 0)
# # 
# # ggplot(data = x,
# #        aes(x = trophlux,
# #            y = slope,
# #            fill = trophlux)) +
# #   geom_violin() +
# #   scale_y_continuous(limits = c(-3, 1))
# # 
# # ggplot(data = filter(y, r.squared > 0.9),
# #        aes(x = trophlux,
# #            y = slope,
# #            fill = trophlux)) +
# #   geom_violin() +
# #   stat_summary() +
# #   scale_y_continuous(limits = c(-3, 1))
# # 
# # 
# # # AOU vs del DIC
# # df_modco2_neppos <- df |>
# #   group_by(year, month, date) |>
# #   drop_na(exDIC_uM, O2ex) |>
# #   filter(NEP_mmolO2m3 > 0) |>
# #   # filter(light > 0) |>
# #   nest() |>
# #   mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
# #          tid = map(mod, broom::tidy),
# #          gl = map(mod, broom::glance))
# # 
# # df_res2 <- select(df_mod2, -data,- mod) |>
# #   unnest(tid) 
# # 
# # # Plot of daily slopes
# # df_res2|>
# #   filter(term == ".$exCO2_uM") |>
# #   ggplot(aes(x = date,
# #              y = estimate)) + 
# #   geom_point() +
# #   # scale_y_continuous(limits = c(-3,3)) +
# #   # scale_x_date(date_breaks = "1 year", date_labels = "%y",
# #   #              name = "") +
# #   stat_smooth() +
# #   theme_bw() +
# #   labs(y = "AOU vs exDIC slope")
# # 
# # 
# # ggplot(data = df_res2_night |>
# #          filter(term == ".$exCO2_uM") |>
# #          mutate(doy = yday(date)) |>
# #          mutate(state = if_else(year < 2012, "planktonic", "benthic")),
# #        aes(x = doy,
# #            y = -estimate,
# #            color = state)) +
# #   stat_smooth() + 
# #   theme_bw()+
# #   # scale_y_continuous(limits = c(0, 2.2)) +
# #   geom_hline(yintercept = 1) +
# #   geom_hline(yintercept = 2, linetype = "dashed")
# # 
# # y_night = select(df_mod2, -mod, -data, -gl) |>
# #   unnest(tid) |>
# #   mutate(term = if_else(grepl("Int", term), "intercept", "slope")) |>
# #   select(-(std.error:p.value)) |>
# #   pivot_wider(names_from = term, values_from = estimate) |>
# #   left_join(select(df_mod, -mod, -data, -tid) |>
# #               hoist(gl, r.squared = "r.squared", p.value = "p.value") |>
# #               select(-gl)) |>
# #   left_join(distinct(df, date, trophlux)) |>
# #   arrange(desc(r.squared), p.value, trophlux)
# # 
# # seacarb::carb(flag = 21,
# #               354,
# #               6,
# #               0,25)
# # 
# # ggplot(data = df_res |>
# #          filter(term == ".$exDIC_uM") |>
# #          mutate(doy = yday(date)),
# #        aes(x = doy,
# #            y = -estimate)) +
# #   stat_smooth() + 
# #   facet_wrap(~year) + 
# #   theme_bw()+
# #   scale_y_continuous(limits = c(0, 3)) +
# #   geom_hline(yintercept = 1) +
# #   geom_hline(yintercept = 2, linetype = "dashed")
# # 
# # 
# # ggplot(data = df_res2 |>
# #          filter(term == ".$exCO2_uM") |>
# #          mutate(doy = yday(date)) |>
# #          mutate(state = if_else(year < 2012, "planktonic", "benthic")),
# #        aes(x = doy,
# #            y = estimate,
# #            linetype = state)) +
# #   stat_smooth(color = "black") + 
# #   theme_classic(base_size = 10) +
# #   scale_linetype(name = "autotrophic state") +
# #   scale_x_continuous(breaks = seq(0, 365, 30)) +
# #   theme(legend.position = c(0.5, 0.8)) +
# #   annotate(geom = "text", x= 60, y = -1.92, label = "calcite precipitaiton") +
# #   geom_hline(yintercept = -1) +
# #   geom_hline(yintercept = -2, linetype = "dashed") +
# #   labs(x = "day of year",
# #        y = expression(daily~exO[2]~-~exCO[2]~slope))
# # 
# # 
# # 
# # 
# # z = df |>
# #   filter(light > 0) |>
# #   group_by(year, month, date, trophlux) |>
# #   summarize(nepo = sum(NEP_mmolO2m3),
# #             nepc = sum(NEP_mmolCm2hr),
# #             dicrem = -sum(delDIC))
# # 
# # ggplot(data = z,
# #        aes(y = nepc,
# #            x = dicrem,
# #            color = year)) +
# #   geom_point() +
# #   facet_wrap(~month) +
# #   # geom_path(size = 2) +
# #   scale_color_viridis_c() +
# #   geom_abline(slope = 1, intercept = 0) +
# #   # geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
# #   theme_classic(base_size = 10) +
# #   geom_hline(yintercept = 0) +
# #   geom_vline(xintercept = 0) +
# #   ggpubr::stat_regline_equation( aes(label =  paste(..eq.label.., 
# #                                                     ..adj.rr.label.., 
# #                                                     sep = "~~~~"))) +
# #   labs(x = expression(DIC[removed]~"("*mmol~m^{-2}~d^{-1}*")"),
# #        y = expression(NEP[C]~"("*mmol~m^{-2}~d^{-1}*")"))
# # 
# # 
# # ggplot(data = z,
# #        aes(y = nepo,
# #            x = dicrem,
# #            color = year)) +
# #   geom_point() +
# #   facet_wrap(~month) +
# #   # geom_path(size = 2) +
# #   scale_color_viridis_c() +
# #   geom_abline(slope = 1, intercept = 0) +
# #   # geom_abline(slope = -2, intercept = 0, linetype = "dashed")+
# #   theme_classic(base_size = 10) +
# #   geom_hline(yintercept = 0) +
# #   geom_vline(xintercept = 0) +
# #   ggpubr::stat_regline_equation( aes(label =  paste(..eq.label.., 
# #                                                     ..adj.rr.label.., 
# #                                                     sep = "~~~~"))) +
# #   labs(x = expression(DIC[removed]~"("*mmol~m^{-2}~d^{-1}*")"),
# #        y = expression(NEP[C]~"("*mmol~m^{-2}~d^{-1}*")"))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ggplot() +
# #   geom_path(data = dplyr::filter(df, date == ymd(20210325)),
# #             aes(x = hr,
# #                 y = HCO3,
# #                 color = NEP_mmolO2m3),
# #             size = 2) +
# #   scale_color_viridis_c() +
# #   stat_smooth(data = dplyr::filter(df, date == ymd(20210325),
# #                                    between(hr, 12, 18)),
# #             aes(x = hr,
# #                 y = HCO3),
# #             method = "lm") +
# #   theme_classic(base_size = 10) +
# #   ggpubr::stat_regline_equation(data = dplyr::filter(df, date == ymd(20210325),
# #                                                      between(hr, 12, 18)),
# #                                 aes(x = hr,
# #                                     y = HCO3,
# #                                     label =  paste(..eq.label.., 
# #                                                    ..adj.rr.label.., 
# #                                                    sep = "~~~~")))
# # 
# # # ALK and DIC
# # t1 <- seacarb::carb(flag = 15, var1 = 0.00207, var2 = 0.00217,
# #                     S = 0, T = 20, warn = "n")
# # # CO2 and ALK
# # t2 <- seacarb::carb(flag = 4, var1 = t1$CO2*0.001, var2 = 0.00207,
# #                     S = 0, T = 20, warn = "n")
# #  #HCO3 qnd ALK
# # t3 <- seacarb::carb(flag = 11, var1 = t2$HCO3, var2 = 0.00207,
# #                     S = 0, T = 20, warn = "n")
# #   t3
# #   
# # 0.000104424*0.8
# 
# # Random examples ---------------------------------------------------------
# df |>
#   filter(date == ymd(20200703)) |>
#          # between(hr, 5, 17)) |>
#   ggplot(aes(x = DIC_uM,# - lag(AT)*1000,
#              y = exO2,#- lag(exO2),
#              color = hr)) + 
#   geom_path(linewidth = 1.5) +
#   scale_color_viridis_c() +
#   theme_classic(base_size = 10) +
#   ggpubr::stat_regline_equation() +
#   # facet_wrap(~month, scales = "free") +
#   # geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
#   # geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
#   # geom_hline(yintercept = 0) +
#   # geom_vline(xintercept = 0) +
#   labs(y = expression(exO[2]~"("*mu*M*")"))
# 
# df |>
#   filter(date == ymd(20200703)) |>
#          # between(hr, 5, 12)) |>
#   ggplot(aes(x = NEP, #- lag(HCO3),
#              y = SpC - lag(SpC), #- lag(O2ex),
#              color = hr)) + 
#   geom_path(linewidth = 1.5) +
#   scale_color_viridis_c() +
#   theme_classic(base_size = 10) +
#   ggpubr::stat_regline_equation()
#   # facet_wrap(~month, scales = "free") +
#   # geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
#   # geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
#   # geom_hline(yintercept = 0) +
#   # geom_vline(xintercept = 0) +
#   # labs(y = expression(exO[2]~"("*mu*M*")"))
# 
# 
# df |>
#   filter(date == ymd_h(2003081000)) |>
#   # between(hr, 5, 17)) |>
#   ggplot(aes(x = DIC, #- lag(HCO3),
#              y = O2ex, #- lag(O2ex),
#              color = hr)) + 
#   geom_path(linewidth = 1.5) +
#   scale_color_viridis_c() +
#   theme_classic(base_size = 10) +
#   ggpubr::stat_regline_equation()
# # facet_wrap(~month, scales = "free") +
# # geom_abline(slope = -2, intercept = 0, linetype = "dashed") +
# # geom_abline(slope = -1, intercept = 0, linetype = "dotted") +
# # geom_hline(yintercept = 0) +
# # geom_vline(xintercept = 0) +
# # labs(y = expression(exO[2]~"("*mu*M*")"))
# 
# seacarb::carb(flag = 24, 
#               var1 = 400,
#               var2 = 1750/1E6,
#               T = 25, S = 0)
# 
# seacarb::carb(flag = 24, 
#               var1 = 400,
#               var2 = 1850/1E6,
#               T = 25, S = 0)
# 
# 
# # 
# # hourly exO2 vs del DIC each day
# df_mod_spc <- df |>
#   group_by(year, month, date) |>
#   # filter(hr > 17 | hr < 8) |>
#   # drop_na(exO2, HCO3_uM) |>
#   mutate(SpC = SpC - lag(SpC)) |>
#   filter(n() > 20) %>%
#   nest() |>
#   mutate(mod = map(data, ~lm(.$SpC ~ .$NEP)),
#          tid = map(mod, broom::tidy),
#          gl = map(mod, broom::glance))
# 
# df_res_spc <- select(df_mod_spc, -data,- mod) |>
#   unnest(tid)
# 
# # p_all
# df_nepspc <- ungroup(df_mod_spc) |>
#   hoist(gl, "r.squared") |>
#   select(-gl, -mod, -data) |>
#   tidyr::unnest(cols = tid) |>
#   filter(term == ".$NEP") |>
#   select(date, slo_spc = estimate, se_spc = std.error, p_spc = p.value, r2_spc =r.squared)
# 
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
