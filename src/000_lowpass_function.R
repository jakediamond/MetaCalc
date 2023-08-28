# Generic lowpass function --------------------------------------------------------
lowpass_fun <- function(data,
                        cutoff_frequency = 24 * 7) # frequency in hours
{
  # Data as dataframe
  data <- as.data.frame(data)
  
  # Name of column with time info
  time_col <- colnames(Filter(function(x) inherits(x, c("POSIXct", "Date")), data))
  
  # Name of column with numeric info
  num_col <- colnames(Filter(function(x) inherits(x, c("double", "numeric")), data))
  
  # Impute all NAs so that there are none with kalman filtering
  data <- data |>
    imputeTS::na_kalman()
  
  # Order the data, just in case
  data <- data[with(data, order(data[, time_col])),]
  
  # Sampling rate [s^-1]
  sr <- 1 / as.numeric(difftime(data[2, time_col], data[1, time_col], 
                                units = "secs"))
  
  # Nyquist frequency = half the sampling rate
  nyq <- sr / 2
  
  # Cutoff frequency (s^-1)
  cutoff <- 1 / (cutoff_frequency * 60 * 60)
  
  # Normalized cutoff frequency for Butterworth filter
  W <- cutoff / nyq
  
  # Butterworth low-pass filter, digital, 2nd order
  myfilter <- signal::butter(2, W, type = 'low', plane = 'z')
  
  # Forward-reverse filter to remove phase-shift 
  # associated with Butterworth filter (must be in vector-form)
  vec <- data[, num_col]
  filtered <- signal::filtfilt(myfilter, vec)
  
  # Filtered data
  data$filtered <- filtered
  data <- data[with(data, order(data[, time_col])),]
}

