# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(tidyverse)
library(tidytable)

# Load data ---------------------------------------------------------------
df_nepco2 <- readRDS(file.path("data", "03_CO2", "NEP_CO2_archetype.RDS"))

# Median estimates
df_mediandaily <- df_nepco2 %>%
  mutate(state = if_else(year < 2005, "planktonic", "benthic"),
         nepco2 = filtered_NEP_mean / filtered_CO2_meanenh) %>%
  group_by(archetype, state) %>%
  drop_na(filtered_CO2_meanenh, filtered_NEP_mean) %>%
  summarize(mCO2 = median(filtered_CO2_meanenh),
            # dCO2 = quantile(filtered_CO2_meanenh),
            mNEP = median(filtered_NEP_mean),
            # dNEP = IQR(filtered_NEP_mean),
            mNEPco2 = median(nepco2))
            # dNEPco2 = IQR(nepco2))

# Year sums of CO2 and NEP
df_y <- df_nepco2 %>%
  mutate(dNEP = (filtered_NEP_2.5 - filtered_NEP_mean) / 1.96,
         dCO2 = (filtered_CO2_97.5enh - filtered_CO2_meanenh) / 1.96) %>%
  group_by(year) %>%
  summarize(totCO2 = sum(filtered_CO2_meanenh, na.rm = T),
            dtotCO2 = sqrt(sum(dCO2^2)),
            totNEP = sum(filtered_NEP_mean),
            dtotNEP = sqrt(sum(dNEP^2)))

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
  group_by(archetype) %>%
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
