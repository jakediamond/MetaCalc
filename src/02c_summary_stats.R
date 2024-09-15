# -------------------------------------
# Author: Jacob Diamond
# Purpose: Calculate basic summary statistics for FCO2 and NEP by trophlux
# Date: 01-03-2023
# -------------------------------------
# Load libraries
library(plotly)
library(htmltools)
library(tidyverse)
library(patchwork)

source(file.path("src", "000_carbonate_functions.R"))

# Load data ---------------------------------------------------------------
df_nepco2 <- readRDS(file.path("data", "daily_trophlux.RDS"))
df <- readRDS(file.path("data", "hourly_data.RDS")) |>
  mutate(wy = if_else(month(date) > 8, year + 1, year))


# Some summaries ----------------------------------------------------------
# Quick histogram of FCO2
x <- df_nepco2$value_FCO2_mean
y <- x[x > 0]
q <- df_nepco2$Q_m3s
hist(x)
hist(y)
hist(q)
hist(log(x))
hist(log(y))
hist(log(q))
mean(y)
exp(mean(log(y)))
shapiro.test(sample(x, size = 5000))
shapiro.test(sample(log(y), size = 5000))

# Summary of CO2
summary(df$CO2_uM)
df %>% transmute(exCO2 = CO2_uM - CO2eq_uM) %>%
  summary(exCO2)

df |>
  filter(NEP > 0, FCO2_enh < 0) |>
  summarize(p=median(pCO2_uatm, na.rm = T))

quantile(df$FCO2_enh, na.rm = T) / 1000 * 12
mean(df$FCO2_enh, na.rm = T) / 1000 * 12

quantile(df$NEP_mean, na.rm = T) / 1000 * 12
mean(df$NEP_mean, na.rm = T) / 1000 * 12

mean(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
sd(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12 / sqrt(283777)
median(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
quantile(df$FCO2_enh, na.rm = T) * 365 / 1000 * 12
mean(df$pCO2_uatm, na.rm = T)
sd(df$pCO2_uatm, na.rm = T)
median(df$pCO2_uatm, na.rm = T)



df |>
  group_by(troph) |>
  summarize(n = n())

df_y <- df_nepco2 |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  filter(n() > 360) |>
  summarize(fco2_y = sum(filtered_FCO2_mean, na.rm = T),
            nep_y = sum(filtered_NEP_mean, na.rm = T))

# Look at discharge sampling bias
df_y_est <- df_nepco2 |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  filter(n() > 360) |>
  nest() %>%
  mutate(Q_mean = map(data, ~filter(., Q_m3s <= 300)),
         Q_75   = map(data, ~filter(., Q_m3s <= 409)),
         Q_80   = map(data, ~filter(., Q_m3s <= 471)),
         Q_90   = map(data, ~filter(., Q_m3s <= 659))) %>%
  mutate(Q_50_30 = map(Q_mean, ~slice_sample(., n = 30) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_50_50 = map(Q_mean, ~slice_sample(., n = 50) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_50_100 = map(Q_mean,~slice_sample(., n = 100) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_75_30 = map(Q_75,     ~slice_sample(., n = 30) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_75_50 = map(Q_75,     ~slice_sample(., n = 50) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_75_100 = map(Q_75,    ~slice_sample(., n = 100) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_80_30 = map(Q_80,     ~slice_sample(., n = 30) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_80_50 = map(Q_80,     ~slice_sample(., n = 50) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_80_100 = map(Q_80,    ~slice_sample(., n = 100) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_90_30 = map(Q_90,     ~slice_sample(., n = 30) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_90_50 = map(Q_90,     ~slice_sample(., n = 50) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T))),
         Q_90_100 = map(Q_90,    ~slice_sample(., n = 100) %>% summarize(mean = mean(value_FCO2_mean, na.rm = T)))) %>%
  select(wy, Q_50_30:Q_90_100) %>%
  unnest(cols = -wy, names_sep = ".") %>%
  rename_with((~stringr::str_replace(.,"^Q_",""))) %>%
  rename_with((~stringr::str_replace(.,".mean",""))) %>%
  pivot_longer(cols = -wy, names_sep = "_", names_to = c("Qfrac", "samp")) %>%
  mutate(year_est = value * 365) %>%
  right_join(df_y) %>%
  ungroup() %>%
  mutate(diff_fco2 = year_est - fco2_y) %>%
  mutate(perdif_fco2 = diff_fco2 / fco2_y * 100,
         Qfrac = as.numeric(Qfrac),
         samp = as.numeric(samp))

ggplot(data = df_y_est,
       aes(x = samp,
          y = perdif_fco2,
          color = Qfrac,
          group = Qfrac)) +
  stat_summary()

mean(df_y$fco2_y) / 1000 * 12
sd(df_y$fco2_y) / 1000 * 12
quantile(df_y$fco2_y) / 1000 * 12

mean(df_nepco2$value_FCO2_mean)*365/ 1000 * 12
median(df_nepco2$value_FCO2_mean)*365/ 1000 * 12

mean(df_nepco2$filtered_FCO2_mean)*365/ 1000 * 12
median(df_nepco2$filtered_FCO2_mean)*365/ 1000 * 12

mean(df_y$nep_y) / 1000 * 12
sd(df_y$nep_y) / 1000 * 12
quantile(df_y$nep_y) / 1000 * 12



# Try with weighted mean for kco2 -----------------------------------------
# Look at discharge sampling bias
df_y_est_wt <- df_nepco2 |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  filter(n() > 360) |>
  left_join(select(df, date, KCO2 = KCO2_mean, z = depth_m) %>%
              group_by(date) %>%
              summarize(kco2 = mean(KCO2) * mean(z)) %>%
              ungroup()) %>%
  drop_na(value_FCO2_mean, kco2) %>%
  nest() %>%
  mutate(Q_mean = map(data, ~filter(., Q_m3s <= 300)),
         Q_75   = map(data, ~filter(., Q_m3s <= 409)),
         Q_80   = map(data, ~filter(., Q_m3s <= 471)),
         Q_90   = map(data, ~filter(., Q_m3s <= 659))) %>%
  mutate(Q_50_30 = map(Q_mean, ~slice_sample(., n = 30) %>% summarize(mean    = weighted.mean(value_FCO2_mean, kco2))), 
         Q_50_50 = map(Q_mean, ~slice_sample(., n = 50) %>% summarize(mean    = weighted.mean(value_FCO2_mean, kco2))), 
         Q_50_100 = map(Q_mean,~slice_sample(., n = 100) %>% summarize(mean   = weighted.mean(value_FCO2_mean, kco2))),
         Q_75_30 = map(Q_75,     ~slice_sample(., n = 30) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_75_50 = map(Q_75,     ~slice_sample(., n = 50) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_75_100 = map(Q_75,    ~slice_sample(., n = 100) %>% summarize(mean = weighted.mean(value_FCO2_mean, kco2))),
         Q_80_30 = map(Q_80,     ~slice_sample(., n = 30) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_80_50 = map(Q_80,     ~slice_sample(., n = 50) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_80_100 = map(Q_80,    ~slice_sample(., n = 100) %>% summarize(mean = weighted.mean(value_FCO2_mean, kco2))),
         Q_90_30 = map(Q_90,     ~slice_sample(., n = 30) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_90_50 = map(Q_90,     ~slice_sample(., n = 50) %>% summarize(mean  = weighted.mean(value_FCO2_mean, kco2))),
         Q_90_100 = map(Q_90,    ~slice_sample(., n = 100) %>% summarize(mean = weighted.mean(value_FCO2_mean, kco2)))) %>%
  select(wy, Q_50_30:Q_90_100) %>%
  unnest(cols = -wy, names_sep = ".") %>%
  rename_with((~stringr::str_replace(.,"^Q_",""))) %>%
  rename_with((~stringr::str_replace(.,".mean",""))) %>%
  pivot_longer(cols = -wy, names_sep = "_", names_to = c("Qfrac", "samp")) %>%
  mutate(year_est = value * 365) %>%
  right_join(df_y) %>%
  ungroup() %>%
  mutate(diff_fco2 = year_est - fco2_y) %>%
  mutate(perdif_fco2 = diff_fco2 / fco2_y * 100,
         Qfrac = as.numeric(Qfrac),
         samp = as.numeric(samp))

ggplot(data = df_y_est_wt,
       aes(x = samp,
           y = perdif_fco2,
           color = Qfrac,
           group = Qfrac)) +
  stat_summary()


# Some more summaries -----------------------------------------------------
df_nepco2 %>%
  ungroup() %>%
  filter(trophlux == "heterotrophic source") %>%
  mutate(rat = -filtered_NEP_mean / filtered_FCO2_mean) %>%
  summarize(int = quantile(rat))

df_nepco2 %>%
  filter(trophlux == "autotrophic sink") %>%
  mutate(rat = -filtered_NEP_mean / filtered_FCO2_mean) %>%
  summarize(int = quantile(rat))

df_nepco2 %>%
  filter(trophlux %in% c("autotrophic sink", "heterotrophic source")) %>%
  mutate(rat = -filtered_NEP_mean / filtered_FCO2_mean) %>%
  summarize(int = quantile(rat))

df_nepco2 %>%
  filter(trophlux %in% c("autotrophic sink", "heterotrophic sink")) %>%
  group_by(year) %>%
  summarize(int = sum(filtered_FCO2_mean)) %>%
  ungroup() %>%
  summarize(quan = quantile(int))

quantile(df_nepco2$filtered_FCO2_mean)
quantile(df_nepco2$filtered_NEP_mean)
quantile(-df_nepco2$filtered_NEP_mean / df_nepco2$filtered_FCO2_mean)

df_nepco2 |>
  rename(FCO2 = filtered_FCO2_mean, NEP = filtered_NEP_mean) |>
  mutate(regime = if_else(year(date) < 2005, "auto", "hetero"),
         troph = if_else(NEP > 0, "auto", "hetero"),
         flux = if_else(FCO2 > 0, "1source", "2sink")) |>
  group_by(regime) |>
  drop_na(FCO2, NEP) |>
  mutate(totFCO2 = sum(FCO2),
         totNEP = sum(NEP),
         len = n()) |>
  ungroup() |>
  group_by(regime, troph, flux) |>
  summarize(count = n(),
            fr = count / mean(len),
            frac = sum(FCO2) / mean(totFCO2),
            nep = median(NEP),
            fco2 = median(FCO2),
            rat = median(-NEP/(FCO2)))


# df |>
#   drop_na(FCO2, NEP) |>
#   mutate(rat = -NEP / FCO2) |>
#   summarize(rat2 = median(rat))
# 
# df |>
#   drop_na(FCO2) |>
#   summarize(meanFCO2 = mean(FCO2),
#             sdFCO2 = sd(FCO2),
#             seFCO2 = sdFCO2 / sqrt(n()),
#             medFCO2 = median(FCO2),
#             kwtmedCO2 = Hmisc::wtd.quantile(FCO2, discharge, probs = 0.5, normwt = T),
#             kwtFCO2 = weighted.mean(FCO2, discharge),
#             sekwtFCO2 = se_magwt(FCO2, discharge))

df_nepco2 |>
  rename(CO2 = filtered_FCO2_mean) |>
  mutate(#dNEP = (filtered_NEP_2.5 - filtered_NEP_mean) / 1.96,
         dCO2 = (filtered_FCO2_97.5 - CO2) / 1.96) |>
  mutate(wtCO2 = 1/dCO2^2) |>
  drop_na(CO2, wtCO2) |>
  # group_by(wy) %>%
  # filter(archetype == "autotrophic_sink") |>
  # filter(archetype == "autotrophic_sink") |>
  # filter(str_detect(archetype, "sink")) |>
  summarize(fCO2 = mean(CO2),
            sdfCO2 = sd(CO2),
            sefCO2 = sdfCO2 / sqrt(n()),
            medfCO2 = median(CO2),
            wtmedfCO2 = Hmisc::wtd.quantile(CO2, wtCO2, probs = 0.5, normwt = T),
            mwtfCO2 = weighted.mean(CO2, wtCO2),
            sewtfCO2 = sqrt(Hmisc::wtd.var(CO2, wtCO2, normwt = TRUE)))
            # sewtFCO2 = sqrt(sum(dCO2^2)))
            # totNEP = sum(filtered_NEP_mean),
            # dtotNEP = sqrt(sum(dNEP^2)))
# 
# 
#   # group_by(archetype) %>%
#   drop_na(filtered_CO2_meanenh, filtered_NEP_mean) %>%
#   summarize(mCO2 = weighted.mean(CO2, wtCO2),
#             dCO2 = sqrt(Hmisc::wtd.var(CO2, wtCO2, normwt = TRUE)),
#             mNEP = weighted.mean(NEP, wtNEP),
#             dNEP = sqrt(Hmisc::wtd.var(NEP, wtNEP, normwt = TRUE)))

# Median estimates
df_nepco2 %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic"),
         nepco2 = filtered_NEP_mean / filtered_CO2_meanenh) %>%
  # filter(archetype %in% c("heterotrophic_source", "autotrophic_sink")) %>%
  # filter(str_detect(archetype, "sink")) %>%
  # filter(archetype == "heterotrophic_source") %>%
  # group_by(archetype, state) %>%
  # group_by(year) |>
  drop_na(filtered_CO2_meanenh, filtered_NEP_mean) %>%
  summarize(
    mCO2 = weighted.mean(filtered_CO2_meanenh, abs(filtered_CO2_meanenh)),
    mNEP = weighted.mean(filtered_NEP_mean, filtered_NEP_mean),
    mNEPco2 = weighted.mean(nepco2, abs(filtered_CO2_meanenh)),
    # mCO2 = median(filtered_CO2_meanenh),
    # dCO2 = quantile(filtered_CO2_meanenh),
    # mNEP = median(filtered_NEP_mean),
    # dNEP = IQR(filtered_NEP_mean),
    medNEPco2 = median(nepco2),
    meanNEPco2 = mean(nepco2)
  )
            # dNEPco2 = IQR(nepco2))

# a <- ggplot(data = filter(df_mediandaily_wt, year != 1997),
#             aes(x = mNEPco2 * 100)) +
#   geom_histogram() + 
#   theme_classic(base_size = 16) + 
#   labs(x = expression(mean~annual~internal~CO[2]~production~"(%)"))
# a

df_y <- df |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  filter(FCO2 > 0) |>
  group_by(wy) %>%
  summarize(totCO2 = sum(FCO2/24, na.rm = T),
            totNEP = sum(NEP, na.rm = T)) |>
  mutate(rat = -totNEP / totCO2)


# Year sums of CO2 and NEP
df_y3 <- df_nepco2 |>
  mutate(dNEP = (filtered_NEP_97.5 - filtered_NEP_mean) / 1.96,
         dCO2 = (filtered_FCO2_97.5 - filtered_FCO2_mean) / 1.96) |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  # filter(value_CO2_mean > 0) |>
  group_by(wy) %>%
  drop_na(filtered_FCO2_mean, filtered_NEP_mean) |>
  # filter(archetype == "autotrophic_sink") |>
  # filter(archetype == "autotrophic_sink") |>
  # filter(str_detect(archetype, "sink")) |>
  summarize(totCO2 = sum(filtered_FCO2_mean, na.rm = T),
            dtotCO2 = sqrt(sum(dCO2^2)),
            totNEP = sum(filtered_NEP_mean),
            dtotNEP = sqrt(sum(dNEP^2)))

x <- df_y3 |>
  mutate(nepco2 = -totNEP/totCO2)


b <- ggplot(data = x,
            aes(x = totCO2 / 1000)) +
  geom_histogram(fill = "#0072B2") + 
  theme_classic(base_size = 16) + 
  labs(x = expression(annual~FCO[2]~production~"("*mol~m^{-2}~y^{-1}*")"))
b

c <- ggplot(data = x,
            aes(x = -totNEP / 1000)) +
  geom_histogram(fill = "#E69F00") + 
  theme_classic(base_size = 16) + 
  labs(x = expression(annual~NEP~"("*mol~m^{-2}~y^{-1}*")"))
c

p <- b + c + a
ggsave(plot = p,
       filename = file.path("results", "annual_budget_summary_histograms.png"),
       units = 'cm',
       dpi = 300,
       height = 12,
       width = 34)

# Yearly sums based on archeytpe and autotrophic state
df_arch <- df_nepco2 %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic")) %>%
  group_by(archetype, year, state) %>%
  summarize(archtot = sum(filtered_CO2_meanenh))


# Summary
df_arch %>%
  left_join(df_y) %>%
  ungroup() %>%
  mutate(prop = archtot / totCO2) %>%
  group_by(archetype, state) %>%
  summarize(mp = median(prop, na.rm = T),
            sp = sd(prop, na.rm = T),
            tot = median(archtot),
            tots = sd(archtot))

# Proportions of time
df_y2 <- df_nepco2 %>%
  group_by(year) %>%
  summarize(n = n())

df_arch2 <- df_nepco2 %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic")) %>%
  group_by(trophlux, state, year) %>%
  summarize(narch = n())

df_arch3 <- df_nepco2 %>%
  group_by(sourcesink, year) %>%
  summarize(narch = n())
z= filter(df_nepco2, year == 2016)
ggplot(z, aes(x = date,
             y = filtered_FCO2_mean,
             color = sourcesink,
              group = 1)) +
  geom_line()

df_arch2 %>%
  left_join(df_y2) %>%
  ungroup() %>%
  group_by(archetype, state) %>%
  mutate(prop = narch / n) %>%
  summarize(mn = median(prop, na.rm = T),
            sn = sd(prop, na.rm = T))

# Compare hourly vs mean daily -------------------------------------------------
df_daily <- df %>%
  select(date, depth_m, CO2eq_uM, KCO2_mean, enh, temp, AT_mM, pH, SpC) %>%
  group_by(date) %>%
  summarize(across(where(is.numeric), mean)) %>%
  mutate(CO2_uM = DIC(temp+273.15, AT_mM, pH, SpC)$CO2_uM) %>%
  mutate(FCO2_mean = (CO2_uM - CO2eq_uM) * KCO2_mean * depth_m *enh) %>%
  left_join(select(df_nepco2, date, troph, sourcesink, FCO2_sum = value_FCO2_mean))

ungroup(df_daily) %>%
  group_by(sourcesink) %>%
  mutate(ratio = FCO2_mean / FCO2_sum,
         bias = FCO2_mean - FCO2_sum) %>%
  summarize(mean_rat = mean(ratio, na.rm = T),
            sd_rat = sd(ratio, na.rm = T),
            mean_bi = mean(bias, na.rm = T),
            sd_bi = sd(bias, na.rm = T))


# How often is NEP 0 or CO2 0 ---------------------------------------------


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

df_meandaily3 <- df_nepco2_clean %>%
  mutate(dNEP = (filtered_NEP_2.5 - filtered_NEP_mean) / 1.96,
         dCO2 = (filtered_CO2_97.5enh - filtered_CO2_meanenh) / 1.96) %>%
  mutate(wtNEP = 1/ dNEP^2,
         wtCO2 = 1/dCO2^2) %>%
  # group_by(archetype) %>%
  drop_na(filtered_CO2_meanenh, filtered_NEP_mean) %>%
  summarize(mCO2 = weighted.mean(CO2, wtCO2),
            dCO2 = sqrt(Hmisc::wtd.var(CO2, wtCO2, normwt = TRUE)),
            mNEP = weighted.mean(NEP, wtNEP),
            dNEP = sqrt(Hmisc::wtd.var(NEP, wtNEP, normwt = TRUE)))

df_meandaily4 <- df_nepco2_clean %>%
  mutate(dNEP = (filtered_NEP_2.5 - filtered_NEP_mean) / 1.96,
         dCO2 = (filtered_CO2_97.5enh - filtered_CO2_meanenh) / 1.96) %>%
  mutate(wtNEP = 1/ dNEP^2,
         wtCO2 = 1/dCO2^2) %>%
  group_by(archetype) %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic")) %>%
  group_by(archetype, state) %>%
  drop_na(filtered_CO2_meanenh, filtered_NEP_mean) %>%
  summarize(mCO2 = Hmisc::wtd.quantile(CO2, wtCO2, probs = 0.5, normwt = T),
            dCO2 = sqrt(Hmisc::wtd.var(CO2, wtCO2, normwt = TRUE)),
            mNEP = Hmisc::wtd.quantile(NEP, NEP, probs = 0.5, normwt = T),
            dNEP = sqrt(Hmisc::wtd.var(NEP, wtNEP, normwt = TRUE)))


df_d <- df %>%
  group_by(year, archetype, date) %>%
  summarize(CO2 = sum(CO2_flux_enh_mmol_m2_hr),
             NEP = sum(NEP_mmolO2m3))

df_meandaily5 <- df_d %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic")) %>%
  group_by(archetype, state) %>%
  summarize(mCO2 = median(CO2, na.rm = T),
           mNEP = median(NEP, na.rm = T))



df_meandaily6 <- df %>%
  group_by(year, archetype, date) %>%
  summarize(CO2 = mean(CO2_daily),
            NEP = mean(NEP_daily)) %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic")) %>%
  group_by(archetype, state) %>%
  summarize(mCO2 = median(CO2, na.rm = T),
            mNEP = median(NEP, na.rm = T))


x= left_join(df_d, select(df_nepco2, date, NEP_d = value_NEP_mean))
mean(-x$NEP_d/ x$NEP, na.rm = T)
sd(x$NEP/-x$NEP_d, na.rm = T)


(ggplot() +
  geom_point(data = df_d,
             aes(x = date,
                 y = NEP / 1000 * 32)) +
  geom_line(data = df_nepco2,
            aes(x = date,
                y = -filtered_NEP_mean / 1000 * 32), color = "red")) %>%
  ggplotly()

quantile(x$NEP/-x$NEP_d, na.rm = T)


# test something different for kco2 weighting -----------------------------
# get the mean exCO2 for each year with kco2 weighting
df_k <- df %>%
  drop_na(exCO2, KCO2_mean, depth_m) %>%
  summarize(mean_exCO2 = weighted.mean(exCO2, KCO2_mean * depth_m),
            mean_kco2 = mean(KCO2_mean * depth_m)) %>%
  mutate(FCO2 = mean_exCO2 * mean_kco2)
mean(df$FCO2, na.rm = T)
