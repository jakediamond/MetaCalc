
# Apply some rules
# First rule to account for when NEP > CO2, which should generally be impossible
# IF –NEP>0 AND CO2 >0 AND -NEP > CO2 flux AND their 95%CI overlap, THEN NEP == CO2 flux 
df_nepco2_clean <- df_nepco2 %>%
  mutate(NEP = if_else(archetype == "heterotrophic_source" &
                         filtered_NEP_mean > filtered_CO2_meanenh &
                         filtered_NEP_97.5 < filtered_CO2_97.5enh,
                       filtered_CO2_meanenh,
                       filtered_NEP_mean))

# Count those conditions
filter(df_nepco2, archetype == "heterotrophic_source" &
         filtered_NEP_mean > filtered_CO2_meanenh) #2723 or 23%

filter(df_nepco2, archetype == "heterotrophic_source" &
         filtered_NEP_mean > filtered_CO2_meanenh &
         filtered_NEP_97.5 < filtered_CO2_97.5enh) #2714 or 99.7% of previous condition

# Second rule(s)
# IF heterotrophic and sink THEN
# IF 95%CI for CO2 contains 0, THEN CO2 = 0; OR
# IF 95% CI for CO2 does not contain 0 AND if NEP 95% CI contains 0, THEN NEP = 0
df_nepco2_clean <- df_nepco2_clean %>%
  mutate(CO2 = if_else(archetype == "heterotrophic_sink" &
                         filtered_CO2_97.5enh > 0,
                       0,
                       filtered_CO2_meanenh)) %>%
  mutate(NEP = if_else(archetype == "heterotrophic_sink" &
                         filtered_CO2_97.5enh < 0 &
                         filtered_NEP_2.5 > 0,
                       0,
                       NEP))

(ggplot(data = df_nepco2_clean,
        aes(x = date)) +
    geom_line(aes(y = CO2), color = "orange") +
    geom_line(aes(y = NEP), color = "blue")) %>%
  ggplotly()


(ggplot(data = df_nepco2_clean,
        aes(x = date,
            y = NEP/CO2,
            color = archetype)) + 
    geom_path()) %>%
  ggplotly()

# summary of that data by month and year
df_arch_sum <- df_nepco2 %>%
  group_by(year, month, archetype) %>%
  summarize(n = n(),
            per = n()/30)

ggplot(data = df_arch_sum,
       aes(x = month,
           y = per,
           color = archetype,
           group = archetype)) +
  stat_summary() +
  stat_summary(geom = "line", fun = "mean") +
  # geom_line() +
  theme_bw() +
  # facet_wrap(~archetype) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_continuous(breaks = seq(1,12,1))+
  labs(y = "fraction of time (-)")

ggsave(file.path("results", "archetype_month_mean_ts.png"),
       dpi = 600,
       units = "cm",
       width = 18.4,
       height = 12)

ggplot(data = df_arch_sum,
       aes(x = year,
           y = n,
           color = archetype,
           group = archetype)) +
  stat_summary(fun = "sum") +
  stat_summary(geom = "line", fun = "sum") +
  # geom_line() +
  # geom_point() +
  theme_bw() +
  facet_wrap(~archetype) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_continuous(breaks = seq(1990,2022,5))+
  labs(y = "number of days (-)")

ggsave(file.path("results", "archetype_year_mean_ts.png"),
       dpi = 600,
       units = "cm",
       width = 18.4,
       height = 12)


# summary of that data by discharge bin and year
df_arch_sum <- df_nepco2 %>%
  mutate(qbin =)
group_by(year, month, archetype) %>%
  summarize(n = n(),
            per = n()/30)


# discharge vs ratio
ggplot(data = df_nepco2,
       aes(x = Q,
           y = ..count..,
           fill = archetype,
           group = archetype)) +
  geom_bar(stat = "bin") +
  # geom_line() +
  # geom_point() +
  theme_bw() +
  facet_wrap(~archetype) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_log10()+
  labs(y = "number of days (-)",
       x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))

ggsave(file.path("results", "archetype_discharge_bins.png"),
       dpi = 600,
       units = "cm",
       width = 18.4,
       height = 12)


# cumulative flux by discharge bin and archetype
df_cum <- df_nepco2 %>%
  mutate(wyear = ifelse(month > 9, year + 1, year),
         qbins = cut(log(Q, base = 10), 30)) %>%
  group_by(qbins, archetype) %>%
  summarize(totflux = sum(filtered_CO2_meanenh, na.rm= T)) %>%
  ungroup() %>%
  mutate(per = totflux / sum(totflux),
         qbinend = 10^as.numeric(str_extract(qbins, "\\b[0-9.]+")))

ggplot(data = df_cum,
       aes(x = qbinend,
           y = per * 100,
           fill = archetype,
           group = archetype)) +
  geom_col() +
  # geom_bar() +
  # geom_line() +
  # geom_point() +
  theme_bw() +
  # facet_wrap(~archetype) +
  # scale_y_log10() +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_log10()+
  theme(legend.position = c(0.2, 0.8)) +
  labs(y = "contribution to total CO2 flux (%)",
       x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))

ggsave(file.path("results", "archetype_discharge_bins_cumulative_flux.png"),
       dpi = 600,
       units = "cm",
       width = 14.3,
       height = 9.2)


ggplot(data = df_nepco2,
       aes(x = Q,
           y = ratio,
           color = archetype,
           group = archetype)) +
  stat_summary_bin(bins = 10) + 
  # geom_bar(stat = "bin") +
  # geom_line() +
  # geom_point() +
  theme_bw() +
  facet_wrap(~archetype, scales = "free") +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_log10() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # scale_y_continuous(limits = c(-10,10)) +
  labs(y = "NEP:CO2 flux (-)",
       x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))

ggsave(file.path("results", "nepco2_discharge_bins_archetype.png"),
       dpi = 600,
       units = "cm",
       width = 18.4,
       height = 12)


df_nepco2 <- df_clean %>%
  select(date, Q, filtered_enh_mean, contains(c("NEP", "CO2"))) %>%
  mutate(across(contains("NEP"), ~.*-1)) %>%
  mutate(ratio_raw = filtered_NEP_mean / filtered_CO2_mean,
         ratio_enh = filtered_NEP_mean / filtered_CO2_meanenh)

x=ggplot(data = df_nepco2, aes(x = date))+
  geom_path(aes(y = filtered_NEP_mean), color = "blue", alpha = 0.7) +
  geom_path(aes(y = filtered_CO2_meanenh), color = "orange", alpha = 0.7)
ggplotly(x)

# Account for known problems
df_nepco2 <- df_nepco2 %>%
  mutate(ratio = if_else((filtered_NEP_mean > filtered_CO2_mean) & 
                           (filtered_CO2_97.5 > filtered_NEP_97.5),
                         1,
                         ratio_raw)
  ) %>%
  mutate(ratio = if_else((filtered_NEP_mean > 0 & filtered_CO2_mean < 0 & filtered_CO2_97.5 > 0),
                         0,
                         ratio))


# Cumulative fluxes per year ----------------------------------------------
# Quick look at cumulative NEP and CO2 over time
df_cumNEP <- df_nepco2 %>%
  mutate(wyear = ifelse(month > 9, year + 1, year)) %>%
  group_by(wyear) %>%
  drop_na() %>%
  mutate(cumNEP = cumsum(filtered_NEP_mean),
         # dcumNEP = sqrt(sum(dNEP^2))*12/32*0.85,
         cumCO2 = cumsum(filtered_CO2_meanenh),
         # dcumCO2 = sqrt(sum(dCO2flux_met_mmolm2d^2)) * 12 / 1000,
         jday = julian(date, origin = min(date)))

# Data ends for plot text
data_ends <- df_cumNEP %>%
  group_by(wyear) %>%
  top_n(1, jday)
data_ends


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
ggplot(data = mutate(df_nepco2, 
                     CO2 = if_else(filtered_CO2_meanenh <0,
                                   0,
                                   filtered_CO2_meanenh),
                     wyear = ifelse(month > 9, year + 1, year)),
       aes(x = CO2,
           color = wyear,
           group = wyear)) +
  gglorenz::stat_lorenz(desc = FALSE, size = 1.5) +
  geom_abline() +
  theme_classic() +
  scale_color_viridis_c() +
  labs(x = "cumulative % of CO2",
       y = "cumulative % of time")

x = filter(df, wyear == 2020) %>%
  mutate(CO2 = if_else(CO2flux_gCm2d <= 0, NA_real_, CO2flux_gCm2d)) %>%
  imputeTS::na_kalman()


erf(sd(log(x$CO2))/2)
Xx = 0.5*(erf((log(x$CO2) - mean(log(x$CO2))) / (sd(log(x$CO2)) * sqrt(2))) + 1)
Yx = 0.5*(erf(((log(x$CO2) - mean(log(x$CO2))) / (sd(log(x$CO2)) * sqrt(2))) - sd(log(x$CO2))/sqrt(2)) + 1)
YX = cdf(cdfinv(Xx) - sd(x$CO2))
plot(Xx, YX)

z = ineq::Lc(x$CO2)
plot(z)

ggplot() +
  gglorenz::stat_lorenz(data = mutate(df, CO2 = if_else(CO2flux_gCm2d <0,0.001,CO2flux_gCm2d)) %>%
                          filter(wyear == 2020),
                        aes(x = log(CO2),
                            color = wyear,
                            group = wyear),
                        desc = FALSE, size = 1.5) +
  geom_line(data = data.frame(x = Xx, y = YX), aes(x = x, y = y)) +
  geom_abline() +
  theme_classic() +
  scale_color_viridis_c() +
  labs(x = "cumulative % of CO2",
       y = "cumulative % of time")
# Find when growing season starts  is each year-------------------------------
df_gs <- df %>%
  mutate(troph = if_else(NEP > 0, "auto", "hetero")) %>%
  group_by(autrop = {autrop = rle(troph); rep(seq_along(autrop$lengths), autrop$lengths)}) %>%
  group_by(autrop, troph, wyear) %>%
  # group_by(wyear, troph) %>%
  # mutate(ltroph = sequence(rle(troph==1)$lengths)) %>% #run length groups of trophic state
  # group_by(wyear, troph) %>%
  # trop = {trop = rle(NEP>0); rep(seq_along(autrop$lengths), autrop$lengths)}) %>%
  summarize(sumNEP = sum(NEP),
            dsumNEP = sqrt(sum(dNEP^2)),
            len = length(NEP),
            mindate = min(date),
            maxdate = max(date),
            minT = min(temp_C),
            meanQ = mean(Q_m3s)) %>%
  group_by(troph, wyear) %>%
  mutate(cumNEP = cumsum(sumNEP)) %>%
  ungroup()

ggplot(data = df_gs,
       aes(x = mindate,
           y = abs(cumNEP),
           color = troph)) +
  geom_line() +
  facet_wrap(~wyear, scales = "free")

# Calculate winter ratio of NEP to CO2 flux -------------------------------
df %>%
  filter.(NEP <= 0) %>%
  filter.(CO2flux_gCm2d > 0) %>%
  filter.(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter.(month < 5 | month > 8) %>%
  filter.(temp_C < 15, Q_m3s > 300) %>%
  filter.(wyear < 2022) %>%
  select(date, wyear, month, NEP_gCm2d, CO2flux_gCm2d) %>%
  pivot_longer(cols = c(NEP_gCm2d, CO2flux_gCm2d)) %>%
  left_join(df_q) %>%
  ggplot(aes(x = wyear_fct)) +
  stat_summary(aes(y = value, fill = name), geom = "col", fun = mean,
               position = "identity") + 
  scale_fill_manual(values = c("#1E88E5", "#FFC107"), name = "",
                    labels = c(expression(CO[2]), "NEP")) + 
  theme_classic() +
  theme(legend.position = c(0.8,0.8),
        axis.title.x = element_blank()) +
  labs(y = expression(mean~flux~to~atmosphere~"("*g~C~m^{-2}~d^{-1}*")"))

ggplot(data = df,
       aes(x = date,
           y = K600_ray)) +
  geom_point() +
  theme_classic()


# df %>%
#   filter(NEP <= 0) %>%
#   filter(CO2flux_gCm2d > 0) %>%
#   filter(NEP_gCm2d < CO2flux_gCm2d) %>%
#   filter(month < 5 | month > 8) %>%
#   filter(temp_C < 15, Q_m3s > 300) %>%
#   filter(wyear < 2022) %>%
#   select(date, wyear, month, NEP_gCm2d, CO2flux_gCm2d) %>%
#   pivot_longer(cols = c(NEP_gCm2d, CO2flux_gCm2d)) %>%
#   ggplot(aes(x = wyear)) +
#   stat_summary(aes(y = value, fill = name), geom = "col", fun = mean,
#                position = "identity") + 
#   scale_fill_manual(values = c("#1E88E5", "#FFC107"), name = "",
#                     labels = c(expression(CO[2]), "NEP")) + 
#   theme_classic() +
#   theme(legend.position = c(0.8,0.8),
#         axis.title.x = element_blank()) +
#   labs(y = expression(mean~flux~to~atmosphere~"("*g~C~m^{-2}~d^{-1}*")"))

df %>%
  filter.(NEP <= 0) %>%
  filter.(CO2flux_gCm2d > 0) %>%
  filter.(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter.(month < 5 | month > 8) %>%
  filter.(temp_C < 15, Q_m3s > 300) %>%
  filter.(wyear < 2022) %>%
  select(date, year, month, NEP_gCm2d, CO2flux_gCm2d) %>%
  mutate(wyear = if_else(month > 9, year+1, year)) %>%
  pivot_longer(cols = c(NEP_gCm2d, CO2flux_gCm2d)) %>%
  group_by(wyear, name) %>%
  summarize(sumC = sum(value, na.rm = T)) %>%
  pivot_wider(values_from = sumC) %>%
  mutate(ratio = NEP_gCm2d / CO2flux_gCm2d) %>%
  ggplot(aes(x = wyear,
             y = ratio)) +
  geom_point() +
  geom_line() +
  stat_smooth(method = MASS::rlm, formula = y ~ x) +
  ggpubr::stat_regline_equation() +
  theme_classic() +
  labs(x = "",
       y = expression(NEP*":"*CO[2]))






df_nepco2 %>%
  filter(filtered_NEP_mean <= 0) %>%
  filter(filtered_CO2_meanenh > 0) %>%
  filter(filtered_NEP_mean < filtered_CO2_meanenh) %>%
  filter(month < 5 | month > 8) %>%
  filter(Q > 100) %>% #temp_C < 15, 
  select(date, year, Q, filtered_NEP_mean, filtered_CO2_meanenh) %>%
  pivot_longer(cols = c(filtered_NEP_mean, filtered_CO2_meanenh)) %>%
  mutate(period = if_else(year < 2005, "phyto", "macro")) %>%
  # mutate(value = if_else(value ==0, 0.001, value)) %>%
  # filter(wyear < 2022) %>%
  ggplot(aes(x = Q,
             y = value,
             color = name)) +
  # linetype = period,
  # group = interaction(name, period))) +
  scale_x_log10() +
  scale_y_log10() +
  # facet_wrap(~wyear) +
  # annotation_logticks() +
  stat_summary_bin() +
  stat_smooth(method = MASS::rlm, formula = y ~ x) +
  # stat_smooth(method = "nls", formula= y~ a*x^(b),
  #             method.args = list(start= c(a = -3.5, b=1.3))) +
  scale_color_manual(values = c("#1E88E5", "#FFC107"), name = "",
                     labels = c(expression(CO[2]), "NEP")) + 
  ggpubr::stat_regline_equation() +
  theme_classic() +
  labs(x = expression(Q~"("*m^3~s^{-1}*")"),
       y = expression(flux~to~atmosphere~"("*g~C~m^{-2}~d^{-1}*")"))



df %>%
  filter.(NEP <= 0) %>%
  filter.(CO2flux_gCm2d > 0) %>%
  filter.(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter.(month < 5 | month > 8) %>%
  filter.(temp_C < 15, Q_m3s > 100) %>%
  select(date, wyear, temp_C, NEP_gCm2d, CO2flux_gCm2d) %>%
  pivot_longer(cols = c(NEP_gCm2d, CO2flux_gCm2d)) %>%
  # filter(wyear < 2022) %>%
  ggplot(aes(x = temp_C,
             y = value,
             color = name)) +
  # scale_x_log10() +
  # scale_y_log10() +
  # annotation_logticks() +
  stat_summary_bin() +
  # stat_smooth(method = MASS::rlm, formula = y ~ x) +
  scale_color_manual(values = c("#1E88E5", "#FFC107"), name = "",
                     labels = c(expression(CO[2]), "NEP")) + 
  # ggpubr::stat_regline_equation() +
  theme_classic() +
  labs(x = expression(stream~temperature~"("*degree*C*")"),
       y = expression(flux~to~atmosphere~"("*g~C~m^{-2}~d^{-1}*")"))



df_hl <- df_use5 %>%
  group_by(year) %>%
  mutate(sequence(rle(NEP < 0)$lengths))


ggplot(data = df,
       aes(x = KCO2_met_mean,
           y = CO2_flux)) +
  stat_summary_bin()
cor(df$KCO2_met_mean, y = df$CO2_flux, use = "pairwise.complete")

# mutate(temptest = if_else(temp_min < 14, 0, 1),
#        qtest = if_else(Q > 300, 0, 1),
#        temptestrun = sequence(rle(temptest==1)$lengths),
#        qtestrun = sequence(rle(qtest==1)$lengths))

select(date, NEP_mean = filtered_NEP_mean, 
       NEP_2.5 = filtered_NEP_2.5, NEP_97.5 = filtered_NEP_97.5,
       CO2_mean = filtered_CO2_meanenh,
       CO2_2.5 = filtered_CO2_2.5enh, CO2_97.5 = filtered_CO2_97.5enh) %>%
  mutate(year = year(date),
         month = month(date)) %>%
  pivot_longer(cols = -c(month, year, date), names_sep = "_", names_to = c("type", "val_type")) %>%
  group_by(year, month, type, val_type) %>%
  summarize(value = mean(value, na.rm = T))

