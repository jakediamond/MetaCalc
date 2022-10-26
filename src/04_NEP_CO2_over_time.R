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

# Load data---------------------------------------------------------------
# Daily C fluxes and all variables
df <- readRDS(file.path("data", "03_CO2", "dampierre_all_daily_data.RDS"))

# Get water years
df <- df %>%
  mutate(wyear = if_else(month > 9, year + 1, year))

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
df_cumNEP <- df %>%
  group_by(wyear) %>%
  drop_na() %>%
  mutate(cumNEP = cumsum(NEP_gCm2d),
         dcumNEP = sqrt(sum(dNEP^2))*12/32*0.85,
         cumCO2 = cumsum(CO2flux_gCm2d),
         dcumCO2 = sqrt(sum(dCO2flux_met_mmolm2d^2)) * 12 / 1000,
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
  geom_ribbon(aes(x = jday, ymin = cumCO2 - 1.96 * dcumCO2, 
                  ymax = cumCO2 + 1.96 * dcumCO2,
                  fill = wyear), alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_c(name = 'water year') +
  scale_fill_viridis_c(name = 'water year') +
  theme_classic() +
  scale_color_viridis_c(name = 'water year') +
  theme(legend.position = c(0.1, 0.85)) +
  labs(x=  "day in water year",
       y = expression("cumulative"~CO[2]~"(g"~O[2]~m^{-2}~d^{-1}*")"))

p_cumco2 + p_cumnep


# Lorenz curves 
ggplot(data = mutate(df, CO2 = if_else(CO2flux_gCm2d <0,0,CO2flux_gCm2d)),
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
  filter(NEP <= 0) %>%
  filter(CO2flux_gCm2d > 0) %>%
  filter(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter(month < 5 | month > 8) %>%
  filter(temp_C < 15, Q_m3s > 300) %>%
  filter(wyear < 2022) %>%
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
  filter(NEP <= 0) %>%
  filter(CO2flux_gCm2d > 0) %>%
  filter(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter(month < 5 | month > 8) %>%
  filter(temp_C < 15, Q_m3s > 300) %>%
  filter(wyear < 2022) %>%
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






df %>%
  filter(NEP <= 0) %>%
  filter(CO2flux_gCm2d > 0) %>%
  filter(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter(month < 5 | month > 8) %>%
  filter(temp_C < 15, Q_m3s > 100) %>%
  select(date, wyear, Q_m3s, NEP_gCm2d, CO2flux_gCm2d) %>%
  pivot_longer(cols = c(NEP_gCm2d, CO2flux_gCm2d)) %>%
  # mutate(value = if_else(value ==0, 0.001, value)) %>%
  # filter(wyear < 2022) %>%
  ggplot(aes(x = Q_m3s,
             y = value,
             color = name)) +
  scale_x_log10() +
  scale_y_log10() +
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
  filter(NEP <= 0) %>%
  filter(CO2flux_gCm2d > 0) %>%
  filter(NEP_gCm2d < CO2flux_gCm2d) %>%
  filter(month < 5 | month > 8) %>%
  filter(temp_C < 15, Q_m3s > 100) %>%
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





# mutate(temptest = if_else(temp_min < 14, 0, 1),
#        qtest = if_else(Q > 300, 0, 1),
#        temptestrun = sequence(rle(temptest==1)$lengths),
#        qtestrun = sequence(rle(qtest==1)$lengths))



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
