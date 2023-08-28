# 
# Purpose: To clean DO time-series for middle Loire and Vienne, reduce outliers
# Author: Jake Diamond
# Date: 25 October 2022
# 

# Load libraries
library(zoo)
library(signal)
library(imputeTS)
library(plotly)
library(lubridate)
library(tidyverse)
library(tidytable)

# Define functions --------------------------------------------------------
# Data cleaning function, finds anomalies/jumps/replaces with NA
# the probability limit for delta_DO for detecting "bad" jumps (default 0.95),
# the multiplier for delta_DO to detect the jump (default 1.5)
# also a minimum threshold for expected DO (default 2)
clean_fun <- function(data,
                      prob = 0.95,
                      mult = 1.5,
                      minDO = 2){
  
  # First run the low pass filter on the data
  # But make sure its in order
  # Order data
  data <- data[with(data, order(datetime)), ]
  # Create hourly time series 
  hr_ts <- data.frame(datetime = seq(floor_date(min(data$datetime), 
                                                unit = "year"),
                                     ceiling_date(max(data$datetime), 
                                                  unit = "day"),
                                     by = "hour"))
  # Combine with data so that every hour has a data point
  data <- right_join(data, hr_ts)
  
  # Low pass filter to get rid of noise
  data <- lowpass_fun(data, variable = DO_mgL)
  
  # # First subset data to get rid of NAs at the front or back end
  data <- na.trim(data)
  
  # Calculate upper 95% average delta DO (and DO) by month 
  # Use this to set a bound on what to expect for big jumps
  del_do <- data %>%
    filter.(between(filtered, minDO, 22)) %>%
    mutate(month = month(datetime),
           year = year(datetime),
           ddo = filtered - lag(filtered)) %>%
    group_by(month) %>%
    summarize(ddo_lim = quantile(abs(ddo), 
                                 probs = prob, 
                                 na.rm = TRUE))
  
  minmax_do <- data %>%
    filter.(between(filtered, minDO, 22)) %>%
    mutate(month = month(datetime),
           year = year(datetime),
           date = date(datetime)) %>%
    group_by(month, date) %>%
    summarize(minDO = min(filtered, na.rm = T),
              maxDO = max(filtered, na.rm = T)) %>%
    ungroup() %>%
    group_by(month) %>%
    summarize(do_limmax = quantile(maxDO, 
                                   probs = prob, 
                                   na.rm = TRUE),
              do_limmin = quantile(minDO, 
                                   probs = 1 - prob, 
                                   na.rm = TRUE))
  
  # Make the data NA where there are big jumps, or just wrong data (anomalies)
  # If the data jump more than 1.5x the 95% value for that time period
  # If DO.obs is less than 2 mg/L (not really possible in this river)
  # Finally make DO NA if an anomaly is detected up to 2 hours ahead or behind
  # this aids in smoothing the data later
  data_natests <- data %>% 
    mutate(month = month(datetime),
           year = year(datetime)) %>%
    left_join(del_do) %>%
    left_join(minmax_do) %>%
    mutate(ddo = filtered - lag(filtered),
           ddotest = if_else(abs(ddo) > mult * ddo_lim | 
                               (ddo == 0 & DO_mgL > 19) |
                               is.na(ddo), "fail", "pass"),
           mindotest = if_else(filtered <= minDO, "fail", "pass"),
           mindomonthtest = if_else(filtered <= do_limmin, "fail", "pass"),
           maxdomonthtest = if_else(filtered >= do_limmax, "fail", "pass")) %>%
    mutate(passfail = paste(ddotest, mindotest, mindomonthtest, maxdomonthtest),
           DO = if_else(str_detect(passfail, "fail"),  NA_real_, filtered )) %>%
    mutate(DO = if_else((is.na(lead(DO)) & lead(DO_mgL) > 19.5 ) |
                          (is.na(lag(DO)) & lag(DO_mgL) > 19.5 ),
                        NA_real_,
                        DO)) %>%
    mutate(DO = if_else((is.na(lead(DO)) & lead(ddotest) == "fail") |
                          (is.na(lag(DO)) & lag(ddotest) == "fail"),
                        NA_real_,
                        DO))
  
  # Check to see where original kalman filter was too filling too much (12 h)
  data_nalengths <- data_natests %>%
    dplyr::group_by(narun = {narun = rle(is.na(DO_mgL)); rep(seq_along(narun$lengths), 
                                                             narun$lengths)}) %>%
    dplyr::mutate(namax = length(is.na(DO_mgL))) %>%
    dplyr::mutate(DO = if_else((is.na(DO_mgL) & namax > 12), NA_real_, DO))
}

# Define lowpass filter function
# Can prescribe cutoff_frequency for low pass bandwidth (default 0.15 
# for hourly data)
lowpass_fun <- function(data, variable,
                        cutoff_frequency = 0.10) {
  # Impute all NAs so that there are none with kalman filtering
  data$DO_imp <- select(data, {{ variable }}) %>%
    na_kalman()
  # Order the data, just in case
  data <- data[with(data, order(datetime)),]
  # Sampling rate [s^-1]
  sr <- 1 / (as.numeric(data$datetime[2] - data$datetime[1]) * 60)
  # Nyquist frequency = half the sampling rate
  nyq <- sr / 2
  # Cutoff frequency (hours^-1)
  cutoff <- 1 / (cutoff_frequency * 60 * 60)
  # Normalized cutoff frequency for Butterworth filter
  W <- cutoff / nyq
  # Butterworth low-pass filter, digital, 2nd order
  myfilter <- butter(2, W, type = 'low', plane = 'z')
  # Forward-reverse filter to remove phase-shift 
  # associated with Butterworth filter (must be in vector-form)
  vec <- as.vector(data$DO_imp)
  filtered <- filtfilt(myfilter, vec)
  # Filtered data
  data$filtered <- filtered
  data <- data[with(data, order(datetime)), ]
}
