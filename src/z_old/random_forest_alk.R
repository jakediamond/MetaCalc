# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

# Load the necessary library
library(leaps)

# Load the dataset
regsubsets.out <-
  regsubsets(Alk_molm3 ~ filtered*Q_m3s*O2+temp+pH+month,
             data = df_mod,
             nbest = 1,       # 1 best model for each number of predictors
             nvmax = NULL,    # NULL for no limit on number of variables
             force.in = NULL, force.out = NULL,
             method = "exhaustive")
regsubsets.out
summary.out <- summary(regsubsets.out)
as.data.frame(summary.out$outmat)
plot(regsubsets.out, scale = "adjr2", main = "Adjusted R^2")
car::subsets(regsubsets.out, statistic="cp", legend = FALSE, min.size = 4, main = "Mallow Cp")
abline(a = 1, b = 1, lty = 2)
which.max(summary.out$adjr2)



df_mod <- df %>%
  group_by(date) %>%
  filter(hour(datetime) == 10)

library(randomForest)
df_mod <- drop_na(df_mod) %>%
  select(-SpC, -SpC_raw, -date, -datetime, -Alk2, -Alk_molkg, -Alk_rf)

df_mod <- select(df_mod, -date)
# Split the data into training and testing sets
set.seed(42) # set the seed for reproducibility
train_indices <- sample(1:nrow(df_mod), size = round(0.7*nrow(df_mod)), replace = FALSE)
train_data <- df_mod[train_indices, ]
test_data <- df_mod[-train_indices, ]

# Train the random forest model
rf_model <- randomForest(Alk_molm3 ~ ., data = train_data, ntree = 500)

# Make predictions on the testing set
pred <- predict(rf_model, newdata = test_data)

# Calculate the accuracy of the model
accuracy <- 1 - sum((pred - test_data$Alk_molm3)^2) / 
  sum((test_data$Alk_molm3 - mean(test_data$Alk_molm3))^2)
cat("Accuracy:", accuracy)


prediction <- predict(rf_model, newdata = df)
df$Alk_rf <- prediction



p2 <- plot_ly(data = df,
             x  = ~datetime) %>%
  add_trace(y = ~ Alk_molm3, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y = ~ Alk_rf, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE) %>%
  add_trace(y = ~ Alk2*1000, type = "scatter", mode='lines', 
            color = I("red"), showlegend = FALSE) %>%
  layout(        xaxis = list(title = ""),
         yaxis = list(range = list(0,450),
                      title = TeX("\\color{#1E88E5}{raw~(mu~S~cm^{-1})}")),
         title = TeX("\\text{Specific Conductance}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p2)
