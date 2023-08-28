# 
# Purpose: To smooth DO/cond./lux time-series, reduce outliers
# Author: Jake Diamond
# Date: September 5, 2019
# 

# Load libraries
library(zoo)
library(signal)
library(imputeTS)
library(lubridate)
library(tidyverse)

# Read in all DO data
df <- readRDS()

# Define functions --------------------------------------------------------
# Data cleaning function, finds anomalies/jumps/replaces with NA
# the probability limit for delta for detecting "bad" jumps (default 0.95),
# the multiplier for delta to detect the jump (default 2)
# also a minimum threshold for expected DO (default 2)
clean_fun <- function(data,
                      type,
                      prob = 0.95,
                      mult = 2){
  # Order data
  data = data[with(data, order(datetime)), ]
  
  # Get normal bounds based on type
  dat_bounds = switch(type,
                      DO = c(1, 17),
                      light = c(0, 120000),
                      conductivity = c(10, 800))
  
  # Calculate upper 95% average delta DO by month and year
  # Use this to set a bound on what to expect for big jumps
  dvalue_lim = data %>%
    mutate(dvalue = value - lag(value)) %>%
    ungroup() %>%
    summarize(dvalue_lim = quantile(abs(dvalue), 
                                    probs = prob, 
                                    na.rm = TRUE)) %>%
    as.numeric()
  
  # Make the data NA where there are big jumps, or just wrong data (anomalies)
  # If the data jump more than 2x the 95% value for that time period
  # If DO.obs is less than 1 mg/L (not really possible in this river)
  # If there was an anomaly from decomposition method
  data = data %>% 
    mutate(dvalue = value - lag(value),
           value2 = ifelse(abs(dvalue) > mult * dvalue_lim | 
                             is.na(dvalue) | 
                             !between(value, dat_bounds[1], dat_bounds[2])
                           , NA,
                           value
           )
    )
  data = data[with(data, order(datetime)), ]
}

# Define lowpass filter function
# Can prescribe cutoff_frequency for low pass bandwidth (default 0.05 15min^-1
# for 15-min data)
lowpass_fun <- function(data, 
                        cutoff_frequency = 0.05) {
  # Re-interpolate all NAs so that there are none with Stineman method
  data$value_an_int <- na_interpolation(data$value2, option = "stine")
  # Order the data, just in case
  data <- data[with(data, order(datetime)),]
  # Sampling rate [s^-1]
  sr <- 1 / (as.numeric(difftime(data$datetime[2], data$datetime[1], units = "secs")))
  # Nyquist frequency = half the sampling rate
  nyq <- sr / 2
  # Cutoff frequency (s^-1)
  cutoff <- cutoff_frequency  / (60 * 15)
  # Normalized cutoff frequency for Butterworth filter
  W <- cutoff / nyq
  # Butterworth low-pass filter, digital, 2nd order
  myfilter <- butter(2, W, type = 'low', plane = 'z')
  # Forward-reverse filter to remove phase-shift 
  # associated with Butterworth filter (must be in vector-form)
  vec <- as.vector(data$value_an_int)
  filtered <- filtfilt(myfilter, vec)
  # Filtered data
  data$filtered <- filtered
  data <- data[with(data, order(datetime)), ]
  rem <- sr / cutoff
  data <- data[-c(1:rem, (nrow(data) - rem):nrow(data)),]
}

# Apply functions ---------------------------------------------------------
# Clean the data and use the lowpass filter
df_n <- df_n %>%
  mutate(clean = map2(data, type, clean_fun),
         filt = map(clean, lowpass_fun))

# Get all in one dataframe
df2 <- df_n %>%
  ungroup() %>%
  select(-data, -clean) %>%
  unnest(cols = filt)

# Finally, add  back NAs for large chunks of missing data
# because the filter can't adequately fill these gaps
df3 <- df2 %>%
  group_by(filename, type, site) %>%
  mutate(na_check = 
           {res = rle(is.na(value2));
           rep(res$values * res$lengths, res$lengths)},
         value_use = ifelse(na_check > 11,
                            NA,
                            filtered)) %>%
  mutate(value_use = if_else(((type == "light" & (value2 == 0 & !is.na(value2))) |
                                type == "light" & value_use < 0), 
                             value2, 
                             value_use)) %>%
  dplyr::ungroup()

# Add back to data
df_final <- select(df3, site:datetime, value_clean = value_use) %>%
  left_join(df) %>%
  select(-value, value = value_clean) %>%
  mutate(temp = if_else((temp < 0 & type != "light") | (temp > 30 & type != "light"), NA_real_, temp))
