# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)
library(tidytable)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete.RDS"))

# Load NEP:CO2 archetypes
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

# Fraction of year by archetype -------------------------------------------
p_doy <- df_nepco2 %>%
  mutate(wyear = if_else(month > 8, year + 1, year),
         origin = ymd(paste0(wyear-1, "1001"))) %>%
  mutate(dowy = as.numeric(interval(origin, date))/86400 + 1,
         datewy = ymd(paste(if_else(month > 8, "1970", "1971"), 
                            month, day(date), sep = "-"))) %>%
  ggplot(aes(x = datewy,
             y = (..count..),
             fill = archetype)) +
  geom_bar(position = "fill", width = 1) +
  theme_bw(base_size = 10) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_date(date_labels = "%b", date_breaks = "61 days",
               # limits 
               expand = expansion(mult = c(0, 0))) +
  # scale_x_continuous(breaks = seq(0, 365, 30), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(y = "32-yr occurence probability (-)",
       x = "") + 
  theme(legend.position = "none", axis.title.x = element_blank())
p_doy


# # summary of that data by month and year
# df_arch_sum <- df_nepco2 %>%
#   mutate(doy = yday(date)) %>%
#   group_by(doy, archetype) %>%
#   summarize(n = n(),
#             per = n()/32)
# 
# ggplot(data = df_arch_sum,
#        aes(x = doy,
#            y = per,
#            color = archetype,
#            group = archetype)) +
#   stat_smooth() +
#   # stat_summary(geom = "line", fun = "mean") +
#   # geom_line() +
#   theme_bw() +
#   # facet_wrap(~archetype) +
#   scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
#   scale_x_continuous(breaks = seq(1,12,1))+
#   labs(y = "fraction of time (-)")
# 
# ggsave(file.path("results", "archetype_month_mean_ts.png"),
#        dpi = 600,
#        units = "cm",
#        width = 18.4,
#        height = 12)


# # discharge vs ratio
# ggplot(data = df_nepco2,
#        aes(x = Q,
#            y = ..count..,
#            fill = archetype,
#            group = archetype)) +
#   geom_bar(stat = "bin") +
#   # geom_line() +
#   # geom_point() +
#   theme_bw() +
#   facet_wrap(~archetype) +
#   scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
#   scale_x_log10()+
#   labs(y = "number of days (-)",
#        x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))
# 
# ggsave(file.path("results", "archetype_discharge_bins.png"),
#        dpi = 600,
#        units = "cm",
#        width = 18.4,
#        height = 12)


# mean annual cumulative flux by discharge bin and archetype
df_cum <- df_nepco2 %>%
  mutate(wyear = ifelse(month > 9, year + 1, year),
         qbins = cut(log10(Q), 30)) %>%
  group_by(wyear) %>%
  mutate(totflux = sum(filtered_CO2_meanenh, na.rm= T)) %>%
  ungroup() %>%
  group_by(qbins, archetype, wyear) %>%
  summarize(flux = sum(filtered_CO2_meanenh, na.rm = T),
            totflux = mean(totflux))

# Proportion of annual fco2 occuring in each archetype and discharge bin
df_fco2_prop <- ungroup(df_cum) %>%
  mutate(per = flux / totflux) %>%
  group_by(qbins, archetype) %>%
  summarize(per = mean(per)) %>%
  mutate(qbinend = 10^as.numeric(str_extract(qbins, "\\b[0-9.]+")),
         archetype = str_replace(archetype, "_", " "))

p_co2_q_hist <- ggplot(data = df_fco2_prop,
       aes(x = qbinend,
           y = per * 100,
           fill = archetype,
           group = archetype)) +
  geom_col() +
  # geom_bar() +
  # geom_line() +
  # geom_point() +
  theme_classic(base_size = 10) +
  geom_vline(xintercept = 300, linetype = "dashed") +
  # facet_wrap(~archetype) +
  # scale_y_log10() +
  scale_fill_manual(name = "trophlux state",
                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  geom_hline(yintercept = 0) +
  theme(legend.position = c(0.24, 0.84),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(y = expression("contribution to annual "*CO[2]~ "flux (%)"),
       x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))
p_co2_q_hist
p_final <- p_doy / p_co2_q_hist + plot_annotation(tag_levels = "a")

ggsave(plot = p_final,
       filename = file.path("results", "figure1.png"),
       dpi = 300,
       units = "cm",
       width = 8.9,
       height = 14)

# Percent of fco2 occuring greater than mean
ungroup(df_cum) %>%
  mutate(qbinend = 10^as.numeric(str_extract(qbins, "\\b[0-9.]+"))) %>%
  mutate(more_less = if_else(qbinend > 300, "more", "less")) %>%
  group_by(wyear, more_less) %>%
  arrange(wyear, more_less) %>%
  mutate(cusum_lessmean = cumsum(flux)) %>%
  filter(cusum_lessmean == max(cusum_lessmean)) %>%
  mutate(prop = cusum_lessmean/totflux) %>%
  ungroup() %>%
  group_by(more_less) %>%
  summarize(mn = mean(prop),
            med = median(prop),
            sd = sd(prop))

mean(log(df_nepco2$Q),na.rm =T)
sd(log(df_nepco2$Q),na.rm =T)
exp(5.35 + 0.84^2/2)

ggplot(data = df_nepco2,
       aes(x = Q)) +
  geom_histogram() +
  scale_x_log10()

# Cumulative fluxes -------------------------------------------------------
# Quick look at cumulative NEP and CO2 over time
df_cumNEP <- df_nepco2 %>%
  mutate(wyear = ifelse(month > 9, year + 1, year)) %>%
  group_by(wyear) %>%
  drop_na() %>%
  mutate(dCO2 = (filtered_CO2_meanenh - filtered_CO2_2.5enh) / 1.96,
         dNEP = (filtered_NEP_mean - filtered_NEP_97.5) / 1.96) %>%
  mutate(cumNEP = cumsum(filtered_NEP_mean),
         dcumNEP = sqrt(sum(dNEP^2)),
         cumCO2 = cumsum(filtered_CO2_meanenh),
         dcumCO2 = sqrt(sum(dCO2^2)),
         jday = julian(date, origin = min(date)))

# Data ends for plot text
data_ends <- df_cumNEP %>%
  group_by(wyear) %>%
  top_n(1, jday)
data_ends

# estimates of cumnep and cumco2
df_tot <- df_cumNEP %>%
  summarize(maxNEP = max(cumNEP, na.rm = T),
            dNEP = mean(dNEP),
            maxCO2 = max(cumCO2, na.rm = T),
            dCO2 = mean(dCO2)) %>%
  mutate(NEP2.5 = maxNEP - 1.96 * dNEP,
         NEP97.5 = maxNEP + 1.96 * dNEP,
         CO22.5 = maxCO2 - 1.96 * dCO2,
         CO297.5 = maxCO2 + 1.96 *dCO2)

ungroup(df_tot) %>%
  mutate(state = if_else(wyear < 2005, "planktonic", "benthic")) %>%
  group_by(state) %>%
  summarize(nep = mean(maxNEP),
            nepsd = sd(maxNEP)  /sqrt(n),
            co2 = mean(maxCO2),
            co2sd = sd(maxCO2))
30*12
46*12/365


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


# Lorenz curves 
(ggplot(data = mutate(df_nepco2, 
                     CO2 = if_else(filtered_CO2_meanenh <0,
                                   0,
                                   filtered_CO2_meanenh),
                     wyear = ifelse(month > 9, year + 1, year)) %>%
         arrange(CO2),
       aes(x = CO2,
           color = wyear,
           group = wyear)) +
  gglorenz::stat_lorenz(desc = FALSE, size = 1.5) +
  geom_abline() +
  theme_classic() +
  scale_color_viridis_c() +
  labs(x = "cumulative % of time",
       y = "cumulative % of CO2")) %>%
  ggplotly()

df_gini_Q <- mutate(df_nepco2, 
       CO2 = if_else(filtered_CO2_meanenh <0,
                     0,
                     filtered_CO2_meanenh),
       wyear = ifelse(month > 9, year + 1, year)) %>%
  group_by(wyear) %>%
  arrange(wyear, date) %>%
  nest() %>%
  mutate(gini = map(data, ~ineq::Gini(.x$Q))) %>%
  select(-data) %>%
  unnest(gini)

ineq::Lc()
df_test <- mutate(df_nepco2, 
                  CO2 = if_else(filtered_CO2_meanenh <0,
                                0,
                                filtered_CO2_meanenh),
                  wyear = ifelse(month > 9, year + 1, year)) %>%
  group_by(wyear) %>%
  arrange(wyear, date) %>%
  drop_na(Q) %>%
  mutate(Qsum = cumsum(Q),
         CO2sum = cumsum(CO2)) %>%
  mutate(Qfrac = Qsum / max(Qsum),
         CO2frac = CO2sum / max(CO2sum))


ggplot(data = df_test,
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


x = df_test %>%
  group_by(wyear) %>%
  slice(which.min(abs(CO2frac - 0.8)))
mean(x$date)
mean(x$Qfrac)
sd(x$Qfrac)
mean(df_gini$gini)
mean(df_gini_Q$gini)
df_gini %>%
  ungroup() %>%
  mutate(state = if_else(wyear < 2006, "planktonic", "benthic")) %>%
  group_by(state) %>%
  summarize(mean = mean(gini),
            sd = sd(gini))
