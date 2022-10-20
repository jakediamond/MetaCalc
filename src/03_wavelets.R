# 
# Purpose: Check out wavelets
# Author: Jake Diamond
# Date: 17 October 2022
# 

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(WaveletComp)
library(tidyverse)

# Load data---------------------------------------------------------------
# Hourly carbonate system
df_carb <- read_csv(file.path("data", "02_CO2", "DAM_carbonate_hourly.csv")) %>%
  mutate(datetime = mdy_hm(datetime),
         O2_mmolm3 = O2 * 1000 / 32,
         O2sat = LakeMetabolizer::o2.at.sat.base(Temperature, 1010),
         O2sat_mmolm3 = O2sat *1000 / 32,
         O2ex = O2_mmolm3 - O2sat_mmolm3,
         CO2ex = `CO2 (mmol/m3)` - `Atmosphere CO2 (mmol/m3)`)

my.w <- analyze.wavelet(df_carb, "O2ex",
                        loess.span = 0,
                        dt = 1/24, dj = 1/100,
                        lowerPeriod = 1/2,
                        upperPeriod = 365,
                        make.pval = FALSE, n.sim = 10)

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))


my.wc <- analyze.coherency(df_carb, my.pair = c("O2ex","CO2ex"),
                           loess.span = 0,
                           dt = 1/24, dj = 1/100,
                           lowerPeriod = 1/2,
                           upperPeriod = 365,
                           make.pval = FALSE)
# Again, we set loess.span = 0 because there is no need to detrend the series; dt = 1/24 means
# we have 24 observations per time unit (1 day, this actually defines the time unit); lowerPeriod =
#   1/2 defines the lowest period to be 12 hours. — To plot the cross-wavelet power spectrum:

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period (days)")
