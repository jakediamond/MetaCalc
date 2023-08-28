# -------------------------------------
# Author: Jake Diamond
# Purpose: Fit regressions for Alkalinity and Ca to create timeseries
# Date: 2023-04-28
# -------------------------------------
# library(randomForest)
library(leaps)
library(caret)
library(patchwork)
library(plotly)
library(tidyverse)

# All hourly data, but only want houlry O2 and pH (and daily discharge and depth)
df_hr <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete.RDS"))

# Cleaned conductivity time series
df_cond <- readRDS(file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))

# Grab sample Alk and Ca data from EDF, need to convert from French degrees
df_alkca <- readxl::read_xlsx(file.path("data", "01_EDF", "raw",
                                  "calcium_alkalinity_dampierre.xlsx")) %>%
  mutate(Ca_molm3 = THCa_degF * 4 / 40.078,
         AT = Alk_degF / 5) %>%
  rename(SpC_EDF = SpC) %>%
  select(-Alk_degF, -THCa_degF)

# Remove outliers
df_alkca <- df_alkca %>%
  pivot_longer(cols = c(Ca_molm3, AT, SpC_EDF)) %>%
  group_by(name) %>%
  mutate(
    #first scale
    scaled = scale(value)) %>% 
  filter(between(scaled, -2.5, +2.5)) %>% #then filter
  select(-scaled) %>%
  pivot_wider(values_from = value, names_from = name)

# Look at relationship, commented out for now
# summary(MASS::rlm(AT~Ca_molm3, data = df_alkca %>%
#                     filter(between(month(date), 7, 8))))
# 
# p_AT_ca <- df_alkca %>%
#   filter(between(month(date), 7, 9)) %>%
#   ggplot(aes(x = Ca_molm3 * 1000,
#              y = AT * 1000)) +
#   geom_point() +
#   stat_smooth(method=function(formula,data,weights=weight) MASS::rlm(formula,
#                                                                 data,
#                                                                 weights=weight,
#                                                                 method="MM"),
#                fullrange=TRUE) +
#   ggpubr::stat_regline_equation(aes(label =  paste(..eq.label..,
#                                                    ..adj.rr.label..,
#                                                    sep = "~~~~"))) +
#   theme_classic(base_size = 10) +
#   labs(x = expression(Ca~"("*mu*M*")"),
#        y = expression(A[T]^{`-`}~"("*mu*M*")"))
# p_AT_ca
# Look at time series of raw data for Calcium and TA, commented out for now
# p_alk <- df_alkca %>%
#   ggplot(aes(x = date, y = AT)) +
#   theme_classic(base_size = 10) +
#   geom_point() +
#   stat_smooth(method = "loess", span = 1/15) +
#   labs(x = "", y = expression(A[T]~"("*mol~m^3*")"))
# 
# p_ca <- df_alkca %>%
#   ggplot(aes(x = date, y = Ca_molm3)) +
#   theme_classic(base_size = 10) +
#   geom_point() +
#   stat_smooth(method = "loess", span = 1/15) +
#   labs(x = "", y = expression(Ca~"("*mol~m^3*")"))
# 
# p_raw <- p_alk + p_ca + plot_annotation(tag_levels = "a")
# p_raw
# ggsave(plot = p_raw,
#        filename = file.path("results", "supp_raw_alk_ca.png"),
#        dpi = 300,
#        units = "cm",
#        height = 12,
#        width = 18.4)

# Get data ready for model ------------------------------------------------
# Daily mean and amplitude of O2, pH, temp
df_mean_amp <- df_hr %>%
  select(year, month, date, datetime, O2, pH, temp, NEP = NEP_mmolO2m3) %>%
  left_join(df_cond) %>%
  select(-datetime) %>%
  group_by(year, month, date) %>%
  summarize(across(where(is.numeric), 
                   list(mean = ~ mean(.x, na.rm = T), 
                        amp = ~ max(.x) - min(.x))))

# Get lags of amplitudes, too. Yesterday's GPP affects the measurement today
df_d <- ungroup(df_mean_amp) %>%
  select(year, month, date, O2_mean, O2_amp, pH_mean, pH_amp, temp_mean, 
         temp_amp, SpC_mean, SpC_amp, NEP_mean, NEP_amp) %>%
  mutate(O2_amp_lag = lag(O2_amp),
         pH_amp_lag = lag(pH_amp),
         temp_amp_lag = lag(temp_amp),
         SpC_amp_lag = lag(SpC_amp),
         NEP_amp_lag = lag(NEP_amp))

# Get data based on the hour most likely of sampling (8am)
df_samptime <- df_hr %>%
  select(year, month, date, datetime, discharge, depth, O2, pH, temp, NEP = NEP_mmolO2m3) %>%
  left_join(df_cond) %>%
  filter(hour(datetime) == 8) %>%
  select(-datetime) %>%
  mutate(discharge = log(discharge)) #log of discharge for modeling

# Get all together for model
df_mod <- left_join(df_alkca, df_d) %>%
  left_join(df_samptime) %>%
  mutate(monthf = as.factor(month))

# Modeling ----------------------------------------------------------------
# Now do AT
df_mod_AT <- df_mod #%>%
  # mutate(Ca_mean = prediction)

ggplot(data = filter(df_mod, abs((SpC_EDF - SpC) / SpC) < 0.02),
       aes(x = SpC_EDF,
           y = AT)) + 
  geom_point() +
  facet_wrap(~month)

ggplot(data = df_mod,
       aes(x = Ca_molm3,
           y = AT,
           color = year)) +
  geom_point() +
  facet_wrap(~month)


summary(lm(AT~SpC, 
           data = filter(df_mod_AT, month == 5)))

summary(lm(SpC~(AT + Ca_molm3)*monthf, 
           data = df_mod_AT))

summary(lm(AT~SpC*discharge*temp*monthf, 
           data = df_mod_AT))

summary(lm(AT~(SpC+discharge+temp)*monthf, 
           data = df_mod_AT))
summary(lm(AT~(SpC + Ca_molm3 + temp)*discharge, 
           data = df_mod_AT))

# Best linear subsets
models_at <- regsubsets(AT~(SpC*discharge + temp) * monthf, really.big = T, # * month
                     data = df_mod_AT, nvmax = 30)
# summary(models_at)
res.sum_at <- summary(models_at)
data.frame(
  Adj.R2 = which.max(res.sum_at$adjr2),
  CP = which.min(res.sum_at$cp),
  BIC = which.min(res.sum_at$bic)
)

get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}

get_model_formula(2, models_at, "AT")

# Get cross validation error function
# get_cv_error <- function(model.formula, data){
#   set.seed(1)
#   train.control <- trainControl(method = "cv", number = 5)
#   cv <- train(model.formula, data = data, method = "lm",
#               trControl = train.control)
#   cv$results$RMSE
# }

# # Compute cross-validation error
# model.ids <- 1:11
# cv.errors <-  map(model.ids, get_model_formula, models_at, "AT") %>%
#   map(get_cv_error, data = df_mod_AT) %>%
#   unlist()
# cv.errors
# which.min(cv.errors)
# models[4]
# lm_mod_at <- lm(AT~ SpC:Ca_molm3 + pH_mean + temp_mean + discharge + Ca_molm3 + temp_mean:monthf + 
#                   discharge:monthf + Ca_molm3:monthf,
#              data = df_mod_AT)
lm_mod_at <- lm(AT ~ SpC*discharge*monthf, data = df_mod_AT)
summary(lm_mod_at)

# 
# # Split the data into training (80%) and testing sets (20%)
# set.seed(42) # set the seed for reproducibility
# train_indices <- sample(1:nrow(df_mod_AT), size = round(0.8*nrow(df_mod_AT)), 
#                         replace = FALSE)
# train_data <- df_mod_AT[train_indices, ]
# test_data <- df_mod_AT[-train_indices, ]
# 
# # Train the random forest model
# rf_model <- randomForest(AT ~ ., data = train_data, ntree = 501, importance = T)
# rf_model
# rf_model$importance
# 
# # Make predictions on the testing set
# pred <- predict(rf_model, newdata = test_data)
# 
# # Calculate the accuracy of the model
# RMSE <- sqrt(mean((pred - test_data$AT)^2))
# cat("RMSE:", RMSE)
# 
# # Predictions from RF
# prediction <- predict(rf_model, newdata = df_mod)
# df_mod$AT_rf <- prediction
# 
# # Plot
# ggplot(data = df_mod,
#        aes(x = AT,
#            y = AT_rf)) +
#   geom_point() +
#   geom_abline()

# # Do the same thing for Calcium
# df_mod_ca <- drop_na(df_mod) %>%
#   select(-date, -year)  %>%
#   filter(between(month, 4, 9))#don't want year to affect anything
#   # mutate(HCO3 = seacarb::carb(flag = 8, pH, AT / 1000, T = temp, pHscale = "F", 
#   #                             S = 0, k1k2 = "m06")$HCO3)
# models_ca <- regsubsets(Ca_molm3~., data = df_mod_ca, nvmax = 8)
# summary(models_ca)
# res.sum_ca <- summary(models_ca)
# data.frame(
#   Adj.R2 = which.max(res.sum_ca$adjr2),
#   CP = which.min(res.sum_ca$cp),
#   BIC = which.min(res.sum_ca$bic)
# )
# 
# # Compute cross-validation error
# model.ids <- 1:8
# cv.errors_ca <-  map(model.ids, get_model_formula, models, "Ca_molm3") %>%
#   map(get_cv_error, data = df_mod_ca) %>%
#   unlist()
# cv.errors_ca
# which.min(cv.errors_ca)
# get_model_formula(4, models_ca, "Ca_molm3")
# lm_mod_ca <- lm(Ca_molm3~ AT + depth + month + SpC_mean,
#              data = df_mod_ca)
# summary(lm_mod_ca)






# # Split the data into training and testing sets
# set.seed(42) # set the seed for reproducibility
# train_indices_ca <- sample(1:nrow(df_mod_ca), size = round(0.8*nrow(df_mod_ca)), replace = FALSE)
# train_data_ca <- df_mod_ca[train_indices_ca, ]
# test_data_ca <- df_mod_ca[-train_indices_ca, ]
# 
# # Train the random forest model
# rf_model_ca <- randomForest(Ca_molm3 ~ ., data = train_data_ca, ntree = 500, importance = T)
# rf_model_ca
# rf_model_ca$importance
# # Make predictions on the testing set
# pred_ca <- predict(rf_model_ca, newdata = test_data_ca)
# 
# # Calculate the accuracy of the model
# RMSE_ca <- sqrt(mean((pred_ca - test_data_ca$Ca_molm3)^2))
# cat("RMSE:", RMSE_ca)

# prediction_ca <- predict(rf_model_ca, newdata = df_mod_ca)
# df_mod_ca$Ca_molm3_rf <- prediction_ca
# 
# ggplot(data = df_mod_ca,
#        aes(x = Ca_molm3_rf,
#            y = Ca_molm3)) +
#   geom_point() +
#   geom_abline()
# summary(lm(test_data_ca$Ca_molm3~pred_ca))

# Get dataset for prediction
df_pred <- select(df_hr, year, month, date, hr, datetime, solartime, discharge,
                  depth, O2, pH, temp, NEP = NEP_mmolO2m3) %>%
  mutate(discharge = log(discharge)) %>%
  left_join(df_cond) %>%
  left_join(df_d) %>%
  mutate(monthf = as.factor(month))

lmmod2 <- lm(AT ~ SpC*monthf, data = df_mod_AT)
summary(lmmod2)
# Do predictions
df_pred$Ca_molm3 <- predict(rf_model, newdata = df_pred)
df_pred <- df_pred %>%
  group_by(date) %>%
  mutate(Ca_mean = mean(Ca, na.rm = T)) %>%
  ungroup()
df_pred$ATreg <- predict(lm_mod_at, newdata = df_pred)
df_pred$ATregsimp <- predict(lmmod2, newdata = df_pred)

# df_pred$ATreg2 <- predict(lm_mod2, newdata = df_pred)
p <- plot_ly(data = df_pred,
                  x  = ~datetime) %>%
  # add_trace(y = ~ Ca_molm3, type = "scatter", mode='lines', color = I("red")) %>%
  add_trace(y = ~ ATregsimp, type = "scatter", mode='lines', color = I("red")) %>%
  add_trace(y = ~ ATreg, type = "scatter", mode='lines', color = I("blue")) %>%
  # color = ~position, showlegend = FALSE) %>%
  layout(xaxis = list(title = ""),
         title = "Total Alkalinity")

htmltools::browsable(p)


df_final <- df_hr %>%
  mutate(AT = df_pred$ATreg,
         SpC = df_pred$SpC)
# Save new hourly data
saveRDS(df_final, file.path("data", "03_CO2", "all_hourly_data_complete_updateSpCAT.RDS"))

ggplot(data = slice_tail(df_final, n = 72),
       aes(x = AT,
           y = Alk_molkg*1000)) +
  geom_point() +
  geom_abline()






# Use a random forest
# Get the data in good format for random forest
df_mod_ca <- df_mod %>%
  select(-date, -year, -month, -AT, -contains("lag"), -contains("O2")) %>%
  drop_na()

# Split the data into training (80%) and testing sets (20%)
set.seed(42) # set the seed for reproducibility
train_indices <- sample(1:nrow(df_mod_ca), size = round(0.8*nrow(df_mod_ca)),
                        replace = FALSE)
train_data <- df_mod_ca[train_indices, ]
test_data <- df_mod_ca[-train_indices, ]

# Train the random forest model
rf_model <- randomForest(Ca_molm3 ~ ., data = train_data, ntree = 501, importance = T)
rf_model
rf_model$importance

# Make predictions on the testing set
pred <- predict(rf_model, newdata = test_data)

# Calculate the accuracy of the model
RMSE <- sqrt(mean((pred - test_data$Ca_molm3)^2))
cat("RMSE:", RMSE)

# Predictions from RF
prediction <- predict(rf_model, newdata = df_mod)
df_ca <- df_mod %>%
  mutate(ca_rf = prediction)

# take a look
ggplot(data = df_ca,
       aes(x = ca_rf,
           y = Ca_molm3)) +
  geom_point() +
  geom_abline()


ggplot(data = df_mod,
       aes(x = SpC,
           y = AT,
           color = Ca_molm3)) +
  geom_point() +
  scale_color_viridis_c() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = '~')),
  ) 
# facet_wrap(~month)