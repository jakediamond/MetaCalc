# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot exO2 vs exCO2 and exDIC
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(patchwork)
library(tidyverse)
library(tidytable)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete.RDS")) %>%
  mutate(trophic = str_split_i(archetype, "_", 1),
         flux = str_split_i(archetype, "_", 2),
         archetype = str_replace(archetype, "_", " ")) %>%
  drop_na(archetype)

# How many samples were under vs oversaturated
sum(df$O2ex > 0, na.rm = T) / nrow(df)

# Overall plot of exDIC vs ex o2
p_all_dic <- ggplot(data = df,
       aes(x = exDIC_uM,
           y = O2ex,
           color = archetype)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.15, 0.15),
        legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  labs(x = expression("exDIC ("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))

# p_all_dic

# Overall plot of exCO2 vs ex o2
p_all_co2 <- ggplot(data = df,
                    aes(x = exCO2_uM,
                        y = O2ex,
                        color = archetype)) +
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

p_all_co2

p_all <- p_all_dic + annotation_custom(ggplotGrob(p_all_co2),
                                                  ymin = 80,
                                                  ymax = 460,
                                                  xmin = 50,
                                                  xmax = 460)

# p_all

# hourly AOU vs del DIC each day
df_mod <- df %>%
  group_by(year, month, date) %>%
  drop_na(exDIC_uM, O2ex) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$O2ex ~ .$exDIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

df_res <- select(df_mod, -data,- mod) %>%
  unnest(tid)

# Plot of slope over year based on ecosystem state
p_slopes <- ggplot(data = df_res %>%
                     filter(term == ".$exDIC_uM") %>%
                     mutate(doy = yday(date)) %>%
                     mutate(state = if_else(year < 2012, "planktonic", "benthic")),
                   aes(x = doy,
                       y = estimate,
                       color = state)) +
  stat_smooth() + 
  theme_classic(base_size = 10) +
  ggsci::scale_color_aaas(name = "autotrophic state") +
  scale_x_continuous(breaks = seq(0, 365, 60)) +
  # scale_y_continuous(breaks = seq(-2.2, 365, 60)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  annotate(geom = "text", x= 200, y = -2.06, label = "calcite precipitation",
           size = 3.5) +
  annotate(geom = "text", x= 200, y = -0.97, label = expression(atop(CO[2]~"/"~HCO[3]^{`-`}, uptake)),  
           size = 3.5) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = -2) +
  labs(x = "day of year",
       y = expression(beta))
p_slopes

# Density plots
p_dens <- ggplot(data = df_res %>%
                   filter(term == ".$exDIC_uM") %>%
                   # filter(r.squared > 0.67,
                   #        p.value < 0.05) %>%
                   mutate(doy = yday(date)) %>%
                   mutate(state = if_else(year < 2012, "planktonic", "benthic")),
                 aes(x = estimate,
                     fill = state)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 10) +
  ggsci::scale_fill_aaas(name = "autotrophic regime") +
  scale_x_continuous(limits = c(-3, 1), breaks = seq(-3, 1, 0.5)) +
  theme(legend.position = c(0.12, 0.8)) +
  # annotate(geom = "text", x= 60, y = -1.92, label = "calcite precipitaiton") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = -2) +
  labs(x = expression(daily~exO[2]~-~exDIC~slope~"("*beta*")"),
       y = "density")

p_dens

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
df_modco2 <- df %>%
  group_by(year, month, date) %>%
  drop_na(exDIC_uM, O2ex) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Summary data
df_resco2 <- select(df_modco2, -data,- mod) %>%
  unnest(tid) 

# Get into tidy and long formats
df_dic <- select(df_mod, -mod, -data, -gl) %>%
  unnest(tid) %>%
  mutate(term = if_else(grepl("Int", term), "intercept", "slope")) %>%
  select(-(statistic:p.value)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error)) %>%
  left_join(select(df_mod, -mod, -data, -tid) %>%
              hoist(gl, r.squared = "r.squared", p.value = "p.value") %>%
              select(-gl)) %>%
  left_join(distinct(df, date, archetype)) %>%
  arrange(desc(r.squared), p.value, archetype) %>%
  select(-std.error_intercept)

df_co2 <- select(df_modco2, -mod, -data, -gl) %>%
  unnest(tid) %>%
  mutate(term = if_else(grepl("Int", term), "intercept", "slope")) %>%
  select(-(statistic:p.value)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error)) %>%
  left_join(select(df_mod, -mod, -data, -tid) %>%
              hoist(gl, r.squared = "r.squared", p.value = "p.value") %>%
              select(-gl)) %>%
  left_join(distinct(df, date, archetype)) %>%
  arrange(desc(r.squared), p.value, archetype) %>%
  select(-std.error_intercept)

# all model info together
df_mods <- df_dic %>%
  rename(int_dic = estimate_intercept,
         slo_dic = estimate_slope,
         r2_dic = r.squared,
         p_dic = p.value,
         se_dic = std.error_slope) %>%
  left_join(df_co2 %>%
              rename(int_co2 = estimate_intercept,
                     slo_co2 = estimate_slope,
                     r2_co2 = r.squared,
                     p_co2 = p.value,
                     se_co2 = std.error_slope))

# Get mean daily LSI
df_lsi <- df %>%
  group_by(date) %>%
  summarize(lsi = mean(LSI, na.rm = T))

# Estimate contributions of DIC to GPP
df_cont <- ungroup(df_mods) %>%
  left_join(df_lsi) %>%
  # Just look at high quality data
  filter(r2_dic > 0.66) %>%
  # Ignore winter days when GPP is too low
  filter(!(archetype == "heterotrophic source" & slo_co2 > -0.6)) %>%
  mutate(regime = if_else(year < 2005, "planktonic", "benthic")) %>%
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
df_cont %>%
  group_by(regime, archetype, source) %>%
  summarize(count = n()) %>%
  mutate(x = count  / sum(count)) %>%
  select(-count) %>%
  pivot_wider(names_from= source, values_from = x)


# Time series of slopes ---------------------------------------------------
# Plot of daily slopes
P_dic_time <- df_dic %>%
  ggplot(aes(x = date,
             y = estimate_slope,
             # alpha = r.squared,
             color = archetype)) + 
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


p_dic_examples <- df_dic %>%
  group_by(archetype, month) %>%
  slice_sample(n = 1) %>%
  left_join(df) %>%
  ggplot(aes(x = exDIC_uM,
             y = O2ex,
             color = hr,
             group = interaction(date, archetype))) + 
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

p_co2_examples <- df_co2 %>%
  group_by(archetype, month) %>%
  slice_sample(n = 1) %>%
  left_join(df) %>%
  ggplot(aes(x = exCO2_uM,
             y = O2ex,
             color = hr,
             group = interaction(date, archetype))) + 
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
# df_mods %>%
#   select(archetype, year, slo_dic, slo_co2) %>%
#   pivot_longer(cols = contains("slo"), names_to = c("slope", "type"), names_sep = "_") %>%
#   filter(between(value, -3, 3)) %>%
#   ggplot(aes(x = value)) +
#   geom_density()+
#   facet_grid(archetype~type)
# 
# ?wtd.var
# 
# a = dplyr::filter(df, date == ymd(20210325)) %>%
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
#        aes(x = archetype,
#            y = slope,
#            fill = archetype)) +
#   geom_violin() +
#   scale_y_continuous(limits = c(-3, 1))
# 
# ggplot(data = filter(y, r.squared > 0.9),
#        aes(x = archetype,
#            y = slope,
#            fill = archetype)) +
#   geom_violin() +
#   stat_summary() +
#   scale_y_continuous(limits = c(-3, 1))
# 
# 
# # AOU vs del DIC
# df_modco2_neppos <- df %>%
#   group_by(year, month, date) %>%
#   drop_na(exDIC_uM, O2ex) %>%
#   filter(NEP_mmolO2m3 > 0) %>%
#   # filter(light > 0) %>%
#   nest() %>%
#   mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
#          tid = map(mod, broom::tidy),
#          gl = map(mod, broom::glance))
# 
# df_res2 <- select(df_mod2, -data,- mod) %>%
#   unnest(tid) 
# 
# # Plot of daily slopes
# df_res2%>%
#   filter(term == ".$exCO2_uM") %>%
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
# ggplot(data = df_res2_night %>%
#          filter(term == ".$exCO2_uM") %>%
#          mutate(doy = yday(date)) %>%
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
# y_night = select(df_mod2, -mod, -data, -gl) %>%
#   unnest(tid) %>%
#   mutate(term = if_else(grepl("Int", term), "intercept", "slope")) %>%
#   select(-(std.error:p.value)) %>%
#   pivot_wider(names_from = term, values_from = estimate) %>%
#   left_join(select(df_mod, -mod, -data, -tid) %>%
#               hoist(gl, r.squared = "r.squared", p.value = "p.value") %>%
#               select(-gl)) %>%
#   left_join(distinct(df, date, archetype)) %>%
#   arrange(desc(r.squared), p.value, archetype)
# 
# seacarb::carb(flag = 21,
#               354,
#               6,
#               0,25)
# 
# ggplot(data = df_res %>%
#          filter(term == ".$exDIC_uM") %>%
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
# ggplot(data = df_res2 %>%
#          filter(term == ".$exCO2_uM") %>%
#          mutate(doy = yday(date)) %>%
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
# z = df %>%
#   filter(light > 0) %>%
#   group_by(year, month, date, archetype) %>%
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
