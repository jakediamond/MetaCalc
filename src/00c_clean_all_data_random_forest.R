# -------------------------------------
# Author: Jake Diamond
# Purpose: Final cleaning of data with random forest filling of NAs
# Date: 27 October 2022
# -------------------------------------

# Load libraries
library(lubdridate)
library(tidyverse)
library(missForest)
library(plotly)
library(tidytable)

# Set seed
set.seed(42)

# Load physical data
df_phys <- readRDS(file.path("data", "temp_discharge_rad_data.RDS"))

# Load pH data
df_pH <- read_csv(file.path("data", "01_EDF", "hourly EDF data", "hourly_pH.csv")) %>%
  rename(dampierre = DAM, belleville = BEL, chinon = CHB) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "pH")

# Load DO data
df_DO <- readRDS(file.path("data", "DO_cleaned_part1.RDS"))

# Combine data and clear workspace
df <- left_join(df_DO, df_phys) %>%
  ungroup() %>%
  left_join(df_pH)

rm(df_DO, df_phys, df_pH)

# Remove DO data again when NA is longer than 3 days, even na_seadesc doesn't work
df <- df %>%
  group_by(site, year, pos, narun) %>%
  mutate(DO_use = if_else((namax > 24 * 3)  & is.na(DO), NA_real_, DO_use)) %>%
  ungroup()

# Remove entire years of missing data
df <- df %>%
  group_by(site, pos, year) %>%
  filter(!(sum(is.na(DO_mgL)) > 300 * 24)) %>%
  ungroup()


# Subsample just to test
x <- filter(df, site == "chinon", between(year, 2005, 2007), pos == "down") %>%
  select(datetime, site, DO_use) %>%
  # pivot_wider(names_from = pos, values_from = DO) %>%
  left_join(filter(df, pos == "up") %>%
              select(datetime, site, temp_C, rad_Wm2, Q_m3s, pH, up = DO_use) %>%
              distinct(),
            by = c("datetime", "site")) %>%
  mutate(hour = hour(datetime),
         month = month(datetime))

mod <- lm(DO_use ~ temp_C + rad_Wm2*Q_m3s + pH, data = x)
summary(mod)
modpred <- predict(mod, newdata = newd)
# plot(mod)
# y = predict(mod, newdata = x)
# plot(modpred)
# lines(x$DO_use)
plot(xpred$DO_use)
lines(modpred)
lines(newd$up, col = "red")
mean(sqrt((y-x$down)^2))

xpred <- filter(df, site == "chinon", between(date, ymd(19990709), ymd(19990822)), 
                pos == "down")
newd <- filter(df, site == "chinon", between(date, ymd(19990729), ymd(19990822)), 
               pos == "down") %>%
  select(datetime, site, DO_use) %>%
  # pivot_wider(names_from = pos, values_from = DO) %>%
  left_join(filter(df, pos == "up") %>%
              select(datetime, site, temp_C, rad_Wm2, Q_m3s, pH, up = DO_use) %>%
              distinct(),
            by = c("datetime", "site")) %>%
  mutate(hour = hour(datetime),
         month = month(datetime))



a = filter(df, site == "dampierre", year == 1999)

xdf <- x %>%
  select(-site, -datetime) %>%
  as.data.frame()
summary(x)
x_im <- missForest(xdf, verbose = TRUE, ntree = 50)

x_im$OOBerror

y <- x_im$ximp




plot_ly(data = bind_cols(y, datetime = x$datetime) %>%
          pivot_longer(cols = c(up, down)), 
        x = ~datetime,
        y = ~value,
        color = ~name) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))



plot_ly(data = filter(df, site == "chinon", pos == "up", between(year, 2005, 2007)) %>%
          select(datetime, DO_mgL, DO_use) %>%
          left_join(select( bind_cols(y, datetime = x$datetime), 
                            datetime, DO)) %>%
          pivot_longer(cols = c(DO, DO_use, DO_mgL)), 
        x = ~datetime,
        y = ~value,
        color = ~name) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))
