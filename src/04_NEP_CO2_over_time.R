# -------------------------------------
# Author: Jake Diamond
# Purpose: Compare NEP and CO2 efflux over time
# Date: 24 October 2022
# -------------------------------------

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(patchwork)
library(tidyverse)

# Source some functions ---------------------------------------------------
source(file.path("src", "000_co2_chem_enhance_functions.R"))
source(file.path("src", "000_lowpass_function.R"))
source(file.path("src", "000_error_propagation_function.R"))
# Load data---------------------------------------------------------------
# Daily C fluxes and all variables
df <- readRDS(file.path("data", "03_CO2", 
                        "CO2_with_uncertainty_dampierre_up_2023-02-21.RDS"))

# Load discharge, pH and temperature, which are needed for the chemical enhancement
df_carb <- read_csv(file.path("data", "03_CO2", 
                                    "DAM_full_correct_daily2.csv")) %>%
  select(date = datetime, Discharge, temp = `Temp (C)`, pH) %>%
  mutate(date = mdy(date))

# Calculate chemical enhancement ------------------------------------------
df_enh <- select(df, date, KCO2 = KCO2_met_mean, depth) %>%
  left_join(df_carb) %>%
  mutate(
    # the chemical enhancement factor
    enh = chem_enh_fun(temp, pH, KCO2, depth)
  )

# Quick look
ggplot(data = df_enh,
       aes(x = date,
           y = enh)) +
  geom_point() +
  scale_y_log10()

ggplot(data = df_enh,
       aes(x = enh))+
  geom_histogram() +
  scale_x_log10()

# Clean the data ----------------------------------------------------------
df_clean <- df %>%
  select(date, GPP_mean, GPP_2.5, GPP_97.5,
                   ER_mean, ER_2.5, ER_97.5, NEP_mean, NEP_2.5, NEP_97.5,
                   CO2_mean = CO2_flux, CO2_2.5 = CO2_flux_2.5, 
                   CO2_97.5 = CO2_flux_97.5) %>%
  left_join(select(df_enh, date, enh_mean = enh)) %>%
  mutate(CO2_meanenh = CO2_mean * enh_mean, 
         CO2_2.5enh = CO2_2.5 * enh_mean, 
         CO2_97.5enh = CO2_97.5 * enh_mean) %>%
  pivot_longer(cols = -date, names_sep = "_", names_to = c("type", "val_type")) %>%
  group_by(type, val_type) %>%
  arrange(type, val_type, date) %>%
  nest() %>%
  mutate(
    # Apply a lowpass function to get rid of excess noise
    filt = map(data, lowpass_fun)
    ) %>%
  unnest(cols = filt) %>%
  select(-data) %>%
  pivot_wider(names_from = c(type, val_type), values_from = c(value, filtered)) %>%
  ungroup() %>%
  drop_na() %>%
  left_join(select(df_carb, date, Q = Discharge)) %>%
  arrange(date)

# Plot an interactive time series -----------------------------------------
p_ts <- plot_ly(data = df_clean, 
        x=~date) %>%
  add_trace(y= ~ -filtered_NEP_mean, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ filtered_CO2_meanenh, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{pH}")),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{[CO_{2}]}~or~\\color{#FFC107}{[HCO_3^{-}]}~(mmol~m^{-3})"),
                      range = list(0, 460)),
         title = TeX("\\text{carbonate system}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_ts)

# try a monthly time series to show the data better
# Nice plotting format first
df_l <- df_clean %>%
  select(date, NEP_mean = filtered_NEP_mean, 
         NEP_2.5 = filtered_NEP_2.5, NEP_97.5 = filtered_NEP_97.5,
         CO2_mean = filtered_CO2_meanenh,
         CO2_2.5 = filtered_CO2_2.5enh, CO2_97.5 = filtered_CO2_97.5enh) %>%
  mutate(year = year(date),
         month = month(date),
         across(contains("NEP"), ~.*-1)) %>%
  pivot_longer(cols = -c(month, year, date), names_sep = "_", 
               names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  mutate(wts = 1/abs(`2.5` - `97.5`)^2) %>%
  group_by(type) %>%
  arrange(type, date) %>%
  mutate(rollmean = slider::slide_dbl(
    .x = cur_data(),
    .f = ~ weighted.mean(
      x = .x$mean,
      w = .x$wts
    ), .before = 90
  ))

# Weighted average time series
p_ts_mo <- ggplot(data = df_l,
                  aes(x = date,
                      color = type)) +
  geom_line(aes(y = rollmean),  linewidth = 1.2) + 
  scale_color_manual(name = "", values = c("#1E88E5", "#FFC107")) +
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  labs(x = "",
       y = expression(CO[2]~flux~"("*mmol~m^{-2}~d^{-1}*")"))
  
p_ts_mo

# Same thing for the ratio, but need to calculate with uncertainty first
df_rat <- df_clean %>%
  mutate(dNEP_mean = (filtered_NEP_mean - filtered_NEP_2.5) / 1.96,
         dCO2_mean = (filtered_CO2_meanenh - filtered_CO2_2.5enh) / 1.96) %>%
  select(date, NEP_mean = filtered_NEP_mean, dNEP_mean,
         CO2_mean = filtered_CO2_meanenh, dCO2_mean) %>%
  mutate_with_error(ratio ~ -NEP_mean/ CO2_mean) %>%
  mutate(wts = 1/dratio^2) %>%
  mutate(rollmean = slider::slide_dbl(
    .x = .,
    .f = ~ weighted.mean(
      x = .x$ratio,
      w = .x$wts
    ), .before = 365
  ))

p_ts_mo_rat <- ggplot(data = df_rat,
                   aes(x = date, y =rollmean)) +
  geom_line(aes(y = rollmean),  linewidth = 1.2) + 
  scale_color_manual(name = "", values = c("#1E88E5", "#FFC107")) +
  theme_classic() + 
  # scale_y_continuous(limits = c(0, 5)) +
  geom_hline(yintercept = 1) + 
  labs(x = "",
       y = expression(NEP~":"*CO[2]~"(-)"))

p_ts_mo_rat

p_ts <- p_ts_mo / p_ts_mo_rat + plot_annotation(tag_levels = "a")
ggsave(plot= p_ts,
       filename = file.path("results", "co2_nep_timeseries.png"),
       dpi = 600,
       units = "cm",
       height = 16,
       width = 14.8)
# Calculate NEP:CO2 -------------------------------------------------------
# Get only the necessary data
df_nepco2 <- df_clean %>%
  select(date, Q, filtered_enh_mean, contains(c("NEP", "CO2")))

# Define based on archetype
df_nepco2 <- df_nepco2 %>%
  mutate(troph = if_else(filtered_NEP_mean > 0, "autotrophic", "heterotrophic"),
         sourcesink = if_else(filtered_CO2_meanenh > 0, "source", "sink"),
         archetype = str_c(troph, sourcesink, sep = "_")) %>%
  mutate(year = year(date),
         month = month(date),
         across(contains("NEP"), ~.*-1)) # get NEP from atmosphere perspective)

# Save this data
saveRDS(df_nepco2, file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

df_nepco2 %>%
  group_by(archetype) %>%
  summarize(n = n(),
            per = n()/nrow(df_nepco2))

# What occurs just before the river becomes a heterotrophic sink?
# First get some group information about each archetype series
df_archgroup <- df_nepco2 %>%
  select(date, archetype) %>%
  group_by(archetype) %>%
  mutate(group = cumsum(c(1, diff(date) > 1))) %>% # group by archeytpe events
  mutate(archgroup = paste0(archetype, group)) %>%
  mutate(arch_l = sequence(rle(archgroup)$lengths)) # length of events

r <- rep(0,nrow(df_nepco2))
df_hetsink <- within(df_archgroup, 
                     identifier <- replace(r,
                                           sapply(which(archetype=="heterotrophic_sink" &
                                                          arch_l == 1),
                                                  `+`, -3:3),
                                           -3:3))
# 71/92 (77%) events are prefaced by autotrophic sinks
filter(ungroup(df_hetsink), identifier == -1) %>%
  group_by(archetype) %>%
  summarize(n = n())


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

ggplot(data = df_nepco2,
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
# Old ---------------------------------------------------------------------
# Estimate NEP in g C m-2 d-1 with respiratory quotient = 0.85 
# But first the problem with NEP when either ER and GPP are 0
df <- df %>%
  # mutate(GPP = if_else(GPP - dGPP < 0, 0, GPP),
  #        ER = if_else(ER + dER > 0, 0, ER),
  #        NEP = if_else(abs(NEP) - dNEP < 0, 0, NEP)) %>%
  mutate(NEP_gCm2d = -NEP * 12/ 32 * 0.85,
         NEP_mmolm2d = -NEP * 1000 / 32,
         CO2flux_gCm2d = CO2flux_met_mmolm2d * 12 / 1000)

# Rank water year by discharge
df_q <- df %>%
  group_by(wyear) %>%
  summarize(qrankmed = median(Q_m3s, na.rm = T),
            qranksum = sum(Q_m3s, na.rm = T)) %>%
  mutate(wyear_fct = fct_reorder(as.factor(wyear), qranksum))

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

# Error functions ---------------------------------------------------------
## if you want the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

# Cumulative error function
cdf <- function(x) 0.5*(erf(x/sqrt(2)) + 1)
# Inverse error function
cdfinv <- function(x) -sqrt(2) *erfcinv(2*x)
