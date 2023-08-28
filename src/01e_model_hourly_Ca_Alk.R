# -------------------------------------
# Author: Jake Diamond
# Purpose: Fit regressions for Alkalinity and Ca to create timeseries
# Date: 2023-04-28
# -------------------------------------
library(randomForest)
# library(leaps)
# library(caret)
# library(patchwork)
library(plotly)
library(htmltools)
library(tidyverse)

# Clean up input data for modeling ----------------------------------------
# All daily CO2 and metabolism data
df_d_co2 <- readRDS(file.path("data", "03_CO2", 
                        "CO2_with_uncertainty_dampierre_up_2023-02-21.RDS")) |>
  mutate(date = date(date))

# Add daily mean discharge, pH, SpC, temperature, keep K600
df_d <- read_csv(file.path("data", "03_CO2", 
                              "DAM_full_correct_daily2.csv")) |>
  select(date = datetime, discharge = Discharge, temp = `Temp (C)`, pH, SpC = Conductivity) |>
  mutate(date = date(mdy(date))) |>
  left_join(select(df_d_co2, date, K600 = K600_met_mean))

# Cleaned hourly conductivity time series
df_cond <- readRDS(file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))

# Hourly O2 data
df_o2 <- readRDS(file.path("data", "02_metabolism", "hourly_inputs.RDS")) |>
  filter(site == "dampierre", pos == "up") |>
  mutate(datetime = floor_date(solar.time, "hours")) |>
  select(datetime, O2_mgL = DO.obs, O2sat_mgL = DO.sat, depth, temp = temp.water) |>
  mutate(O2_uM = O2_mgL * 1000 / 32, O2sat_uM = O2sat_mgL * 1000 / 32,
         date = date(datetime))

# Calculate hourly NEP
df_o2 <- df_o2 |>
  left_join(select(df_d, date, K600)) |>
  mutate(Sc_O2 = 1801 - 120.1 * temp + 3.782 * temp^2 - 0.0476 * temp^3,
         KO2 = K600 / (600 / Sc_O2)^-0.5) |>
  arrange(datetime) |>
  # smooth the K and depth data, 3 day moving average
  mutate(K600_smooth = slider::slide_dbl(K600, mean, .before = 3*24, na.rm = T),
         KO2_smooth = slider::slide_dbl(KO2, mean, .before = 3*24, na.rm = T),
         depth_smooth = slider::slide_dbl(depth, mean, .before = 3*24, na.rm = T),
         NEP = ((O2_uM - lag(O2_uM)) - (KO2_smooth / 24) * 
                  (O2sat_uM - O2_uM)) * depth_smooth)

# Grab sample Alk and Ca data from EDF, need to convert from French degrees
df_alkca <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                  "calcium_alkalinity_dampierre.xlsx")) |>
  mutate(Ca_molm3 = THCa_degF * 4 / 40.078,
         AT = Alk_degF / 5) |>
  rename(SpC_EDF = SpC) |>
  select(-Alk_degF, -THCa_degF)

# Remove outliers
df_alkca <- df_alkca |>
  pivot_longer(cols = c(Ca_molm3, AT, SpC_EDF)) |>
  group_by(name) |>
  mutate(
    #first scale
    scaled = scale(value)) |> 
  filter(between(scaled, -2.5, +2.5)) |> #then filter
  select(-scaled) |>
  pivot_wider(values_from = value, names_from = name)

# Get data ready for model ------------------------------------------------
# Daily mean NEP, O2, and SpC (from cleaned data set)
df_d_mean <- df_o2 |>
  select(datetime, date, O2_uM, NEP) |>
  left_join(df_cond) |>
  select(-datetime) |>
  group_by(date) |>
  summarize(across(where(is.numeric), mean))

# Add the other mean daily data, fill in for missing SpC
df_d_mean <- df_d_mean |>
  complete(date = seq.Date(min(date), max(date), by="1 day")) |>
  rows_patch(select(df_d, date, SpC), by = "date") |>
  left_join(select(df_d, -SpC))

# Get all together for model
df_mod <- left_join(df_alkca, df_d_mean) |>
  select(-SpC_EDF) |> #redundant with measured SpC
  imputeTS::na_kalman() # fill the small amount of missing values

# Modeling ----------------------------------------------------------------
# First get Calcium. Use a random forest
# Show the poor fit of calcium to SpC
summary(lm(Ca_molm3 ~ SpC, data = df_mod))

# Split the data into training and testing sets
set.seed(42) # set the seed for reproducibility
df_mod_ca <- drop_na(df_mod) |>
  select(-date, -AT)

# Training and testing data
train_indices_ca <- sample(1:nrow(df_mod_ca), size = round(0.8*nrow(df_mod_ca)), replace = FALSE)
train_data_ca <- df_mod_ca[train_indices_ca, ]
test_data_ca <- df_mod_ca[-train_indices_ca, ]
 
# Train the random forest model
rf_model_ca <- randomForest(Ca_molm3 ~ ., data = train_data_ca, ntree = 500, importance = T)
rf_model_ca
rf_model_ca$importance

# Make predictions on the testing set
pred_ca <- predict(rf_model_ca, newdata = test_data_ca)
 
# Calculate the accuracy of the model
RMSE_ca <- sqrt(mean((pred_ca - test_data_ca$Ca_molm3)^2))
cat("RMSE:", RMSE_ca)

# Get the predicted Calcium data for the model set
prediction_ca <- predict(rf_model_ca, newdata = df_mod)


# Now do modeling for AT
df_mod_AT <- df_mod |>
  mutate(Ca_rf = prediction_ca)

# function to get RMSE from lm
rmse_fun <-function(mod) {
  RSS <- c(crossprod(mod$residuals))
  # Mean squared error:
  MSE <- RSS / length(mod$residuals)
  # Root MSE:
  RMSE <- sqrt(MSE)
  return(RMSE)
}

# Quick look at simple regression
AT_mod_spc <- lm(AT~SpC, data = df_mod_AT)
summary(AT_mod_spc)
rmse_fun(AT_mod_spc)

# Fit full model
AT_mod_all <- lm(AT ~ SpC + log(discharge) + temp + Ca_rf, data = df_mod_AT)
summary(AT_mod_all)
rmse_fun(AT_mod_all)

# Get dataset for prediction of daily Ca
df_pred_d <- df_d_mean |>
  mutate(Ca_rf = predict(rf_model_ca, newdata = df_d_mean)) |>
  imputeTS::na_kalman()

# Simple model just based on SpC that varies by month
AT_mod_simp <- lm(AT ~ SpC*monthf, data = mutate(df_mod_AT, 
                                                 monthf = as.factor(month(date))))
summary(AT_mod_simp)

# Get data for hourly AT predictions
df_pred_hr <- df_o2 |>
  left_join(df_cond) |>
  select(-temp) |> #only want mean daily temp, not hourly
  complete(datetime = seq.POSIXt(min(datetime), max(datetime), by="1 hour")) |>
  mutate(date = date(datetime)) |>
  rows_patch(select(df_pred_d, date, SpC), by = "date") |>
  left_join(select(df_pred_d, date, Ca_rf, discharge, temp, pH)) |>
  arrange(date, datetime) |>
  # Smooth with a 3 day moving window for smoother transitions across days
  mutate(Ca_rf = slider::slide_dbl(Ca_rf, mean, .before = 3*24, na.rm = T),
         temp = slider::slide_dbl(temp, mean, .before = 3*24, na.rm = T),
         pH = slider::slide_dbl(pH, mean, .before = 3*24, na.rm = T),
         discharge = slider::slide_dbl(discharge, mean, .before = 3*24, na.rm = T))

# Get predictions
df_pred_hr$ATreg <- predict(AT_mod_all, newdata = df_pred_hr)
df_pred_hr$ATregsimp <- predict(AT_mod_simp, newdata =  mutate(df_pred_hr, 
                                                          monthf = as.factor(month(date))))

# Take a look at data and compare with simple SpC predictions ("red")
p <- plot_ly(data = df_pred_hr,
                  x  = ~datetime) |>
  # add_trace(y = ~ Ca_molm3, type = "scatter", mode='lines', color = I("red")) |>
  add_trace(y = ~ ATregsimp, type = "scatter", mode='lines', color = I("red")) |>
  add_trace(y = ~ ATreg, type = "scatter", mode='lines', color = I("blue")) |>
  # color = ~position, showlegend = FALSE) |>
  layout(xaxis = list(title = ""),
         title = "Total Alkalinity")

htmltools::browsable(p)

# Final dataframe to save
df_final <- df_pred_hr |>
  select(-ATregsimp, -temp, -pH) |>
  rename(Ca = Ca_rf, AT = ATreg) |>
  left_join(select(df_o2, datetime, temp))

# Save new hourly data
saveRDS(df_final, file.path("data", "03_CO2", "hourly_SpC_AT_Ca_O2system.RDS"))


