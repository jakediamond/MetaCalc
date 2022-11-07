# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
library(tsforecast)
library(lubridate)
library(tidyverse)

df_DO <- readRDS(file.path("data", "DO_cleaned_part1.RDS"))

x = filter(df_DO, site == "belleville", year == 1993) %>%
  select(datetime, DO, pos) %>%
  pivot_wider(names_from = pos, values_from = DO) %>%
  mutate(date = seq(from = ym('199001'), length.out = 8760, by='months'))

?seq

ym(199001)
# For UNIVARIATE forecasting
ts_forecast_data_univariate <- tstools::initialize_ts_forecast_data(
  data = x, 
  # Indicate which Date column corresponds to the time component
  date_col = "datetime", 
  # Indicate which column of values should be forecasted
  col_of_interest = "up",
  # OPTIONAL: Indicate which column(s) should be used to create groups to forecast separately
  # group_cols = c("state", "oil_company")
)
head(ts_forecast_data_univariate)


# For MULTIVARIATE forecasting
ts_forecast_data_multivariate <- tstools::initialize_ts_forecast_data(
  data = x, 
  date_col = "date", 
  col_of_interest = "up", 
  # group_cols = c("state", "oil_company"),
  # For multivariate forecasting, indicate which column(s) should be used as external regressors
  xreg_cols = "down"
)
head(ts_forecast_data_multivariate)




main_forecasting_table <- update_main_forecasting_table(
  # Change the path to write to a different location (instead of the current working directory)
  file_path = file.path("data", "example_forecast_file.rds"), 
  # Using the multivariate forecast data for this example, but univariate is also possible 
  data = ts_forecast_data_multivariate, 
  # Default is quarterly (every 3 months) and yearly (every 12 months) seasonality
  seasonal_periods = c(24), 
  # Default is 24 months (2 years) of miminum training period, but 180 month is used here to limit the runtime
  min_train_periods = 24 * 4, 
  # By default the training data is not limited
  max_train_periods = Inf, 
  # 12 months ahead forecasting is the default forecast horizon
  periods_ahead = 24*60, 
  # Limit the forecast methods to only two methods to limit the runtime
  fc_methods = c("linear", "prophet"),
  # Set to TRUE to get updates on progress
  verbose = FALSE
)
head(main_forecasting_table)
