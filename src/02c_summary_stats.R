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

# Load data ---------------------------------------------------------------
df_nepco2 <- readRDS(file.path("data", "daily_trophlux.RDS"))
df <- readRDS(file.path("data", "hourly_data.RDS")) |>
  mutate(wy = if_else(month(date) > 8, year + 1, year))

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
  group_by(archetype, state, year) %>%
  summarize(narch = n())

df_arch2 %>%
  left_join(df_y2) %>%
  ungroup() %>%
  group_by(archetype, state) %>%
  mutate(prop = narch / n) %>%
  summarize(mn = median(prop, na.rm = T),
            sn = sd(prop, na.rm = T))

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
