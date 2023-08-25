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
library(biwavelet)
library(tidyverse)

# Load data---------------------------------------------------------------
# Hourly carbonate system
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data.RDS"))

# my.w <- analyze.wavelet(df, "O2ex",
#                         loess.span = 0,
#                         dt = 1/24, dj = 1/100,
#                         lowerPeriod = 1/2,
#                         upperPeriod = 365,
#                         make.pval = FALSE, n.sim = 10)

# wt.image(my.w, color.key = "quantile", n.levels = 250,
#          legend.params = list(lab = "wavelet power levels", mar = 4.7))

df2 <- filter(df, year(datetime) == 1994)

my.wc <- analyze.coherency(df2, my.pair = c("O2ex","CO2ex"),
                           loess.span = 0,
                           dt = 1/24, dj = 1/100,
                           lowerPeriod = 1/2,
                           upperPeriod = 365,
                           make.pval = TRUE)

bw <- wtc(d1 = as.matrix(data.frame(t = df2$datetime, co2 = df2$O2ex)), 
          d2 = as.matrix(data.frame(t = df2$datetime, co2 = df2$CO2ex)),
          nrands = 10)

# Again, we set loess.span = 0 because there is no need to detrend the series; dt = 1/24 means
# we have 24 observations per time unit (1 day, this actually defines the time unit); lowerPeriod =
#   1/2 defines the lowest period to be 12 hours. — To plot the cross-wavelet power spectrum:
at <- seq(1, nrow(df2), by = 24*31) # every active day at 00:00:00
labels <- strftime(as.POSIXct(df2$datetime[at],
                              format="%F %T", tz = "GMT"), format ="%b %d")

png(file.path("results", "coherency_exO2_exCO2_1994.png"),
    width = 18, height = 16, units = "cm", res = 1200)
wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period (days)",
         spec.time.axis = list(at = at, labels = labels),
         main = "1994–wavelet coherency exO2 and exCO2")
dev.off()
