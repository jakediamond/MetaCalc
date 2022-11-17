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

# Load all DO data---------------------------------------------------------------
# Data directory
data_dir <- file.path("data", "01_EDF", "raw", "dissolved_oxygen")

# files
files <- fs::dir_ls(data_dir, regexp = "\\.txt$")

# Load DO data
df <- files %>% 
  map_dfr(read.table, skip = 2, header = T,
          col.names = c("datetime", "DO_mgL", "QC", "niv_mod",
                        "no_v_data", "validity", "niv_coh"),
          .id = "filename") %>%
  mutate(datetime = dmy_hm(datetime))

# Get names of each site and position (upstream or downstream)
df <- df %>%
  mutate(pos = if_else(str_detect(filename, "3_"), "down", "up"),
         site = case_when(str_detect(filename, "511") ~ "civaux",
                          str_detect(filename, "581") ~ "belleville",
                          str_detect(filename, "582") ~ "dampierre",
                          str_detect(filename, "583") ~ "stlaurent",
                          str_detect(filename, "584") ~ "chinon")
         )

# Some simple data cleaning, DO can't be less than 0 or more than 30
df <- df %>%
  mutate(DO_mgL = if_else(DO_mgL < 0 | DO_mgL > 30, NA_real_, DO_mgL))

# Summary of the data
df %>%
  dplyr::group_by(site, pos) %>%
  dplyr::summarize(nas = sum(is.na(DO_mgL)),
            dat = sum(!is.na(DO_mgL)),
            tot = sum(dat + nas),
            per = nas / (dat + nas),
            mindat = min(datetime),
            maxdat = max(datetime))

# Plot all sites upstream and downstream
p_DO_raw <- plot_ly(data = dplyr::filter(df, site != "stlaurent"), 
                    x = ~datetime,
                    y = ~DO_mgL,
                    color = ~site,
                    linetype = ~pos,
                    colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
          #             title = TeX("\\text{pH}"),
           #            range = list(6, 11)),
         xaxis = list(title = ""),
         yaxis = list(title = "DO (mg/L)",
                      range = list(0, 30)),
         title = "raw DO data")

# browsable(p_DO_raw)

# Nest data by site and period
# Split into the three periods based on obs. 1993-2003, 2004-2010, 2011-2022
df_n <- df %>%
  select(site, pos, datetime, DO_mgL) %>%
  mutate(period = case_when(datetime < ymd_hms("2004-01-01 00:00:00") ~ 1,
                            datetime > ymd_hms("2010-12-31 00:00:00") ~ 3,
                            TRUE ~ 2)) %>%
  group_by(site, pos, period) %>%
  nest()

# Remove original data to keep workspace clean
rm(df)
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
  data <- lowpass_fun(data)
  
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
  
  # Finally fill the NA where possible, kalman less than 12 h
  data_final <- data_natests %>%
    dplyr::group_by(narun = {narun = rle(is.na(DO)); rep(seq_along(narun$lengths), narun$lengths)}) %>%
    dplyr::mutate(namax = length(is.na(DO))) %>%
    group_by(year) %>%
    filter.(!(sum(is.na(DO_mgL)) > 180 * 24)) %>%
    ungroup()# %>%
    # mutate(DO_use = na_kalman(DO, maxgap = 12))
  
  # Fill in NAs as best as possible with kalman filtering
  # do_ts <- ts(data_final$DO_use, deltat = 1/(365*24))
  # do_ts_clean <- imputeTS::na_seasplit(do_ts, algorithm = "kalman")
  # data_final$DO_use <- as.numeric(do_ts_clean)
  # data_final
}

# Define lowpass filter function
# Can prescribe cutoff_frequency for low pass bandwidth (default 0.15 
# for hourly data)
lowpass_fun <- function(data, 
                        cutoff_frequency = 0.10) {
  # Re-interpolate all NAs so that there are none with kalman filtering
  data$DO_an_int <- na_kalman(data$DO)
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
  vec <- as.vector(data$DO_an_int)
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
  filter.(site != "stlaurent") %>% #not enough data at st lauren to care
  mutate(clean = map(data, clean_fun))

# Get all in one dataframe
df_DO <- unnest(df_n, clean)
a = filter.(df_DO, site == "dampierre", year(datetime) == 1996)
# Plot all sites upstream and downstream
p_DO_clean <- plot_ly(data = dplyr::filter(df_DO2, site == "dampierre"), 
                    x = ~datetime,
                    y = ~DO,
                    color = ~pos) %>%
                    # colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)),
    title = "raw DO data")

htmltools::browsable(p_DO_clean)

# Save data
saveRDS(df_DO, file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))

# A second cleaning -------------------------------------------------------
df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS")) %>%
  select(site, pos, year, datetime, period, DO_mgL = DO_use, namax, narun, DO)

# Remove entire years of missing data (more than a half year)
# Remove DO data again when NA is longer than 3 days, initial na_kalman didn't work
df_DO <- df_DO %>%
  group_by(site, pos, year) %>%
  filter.(!(sum(is.na(DO_mgL)) > 180 * 24)) %>%
  ungroup() %>%
  group_by(site, year, pos, narun) %>%
  mutate(DO_mgL = if_else((namax > 24 * 3)  & is.na(DO), NA_real_, DO_mgL)) %>%
  ungroup() %>%
  select(-namax, -narun, -DO)

# Nest data and apply
df_n <- df_DO %>%
  group_by(site, pos, period) %>%
  nest() %>%
  mutate(clean = map(data, possibly(clean_fun, NA_real_)))

# Get all in one dataframe
df_clean <- filter.(df_n, !is.na(clean)) %>%
  unnest(clean)

x = pluck(df_n, 4,2)



# Plot all daily min and max for trends -----------------------------------
df_minmax <- read_csv(file.path("data", "01_EDF", "hourly EDF data", 
                                "hourly_DO.csv"))

df_minmax <- df_minmax %>% 
  rename(belleville = BEL, dampierre = DAM, chinon = CHB) %>%
  mutate(datetime = mdy_hm(datetime), date = date(datetime)) %>%
  pivot_longer(cols = c(belleville, dampierre, chinon), names_to = "site") %>%
  group_by(site, date) %>%
  summarize(max = max(value, na.rm = T),
            min = min(value, na.rm = T)) %>%
  pivot_longer(cols = c(min, max)) %>%
  mutate(moyr = ymd(paste(year(date), month(date), "01")),
         site_f = factor(site)) %>%
  mutate(site_f = fct_relevel(site_f, "belleville", "dampierre", "chinon"))

p_minmax <- ggplot(df_minmax,
       aes(x = moyr,
           y = value,
           color = name)) +
  stat_summary(fun = median, geom = "line", size = 1.2) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  # stat_smooth(method = "loess", span = 0.1, color = "black") +
  facet_grid(rows = vars(site_f)) +
  theme_classic(base_size = 18) +
  labs(x = "",
       y = expression("oxygène dissous (mg "*L^{-1}*")")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave(plot = p_minmax,
       filename = file.path("results", "DO_minmax.png"),
       dpi = 1200,
       units = "cm",
       height = 20,
       width = 28)

# Seasonality
df_minmax_seas <- df_minmax %>%
  mutate(month = month(date),
         season = case_when(month %in% c(3,4,5) ~ "spring",
                            month %in% c(6,7,8) ~ "summer",
                            month %in% c(9,10,11) ~ "autumn",
                            TRUE ~ "winter"),
         year = year(date)) %>%
  mutate(seaf = factor(season)) %>%
  mutate(seaf = fct_relevel(seaf, "spring", "summer", "autumn", "winter")) %>%
  group_by(year, season, name, site_f ,seaf) %>%
  summarize(val = mean(value, na.rm = T))

p_seas <- ggplot(data = df_minmax_seas,
                 aes(x = year,
                     y = val,
                     color = name)) +
  geom_point() + geom_line() +
  # stat_summary() + stat_summary(geom = "line") +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(formula=y~x) +
  ggpubr::stat_cor(aes(label = ..p.label..), label.x = 2010) +
  facet_grid(rows = vars(site_f), cols = vars(seaf)) +
  theme_classic(base_size = 18) +
  labs(x = "",
       y = expression("oxygène dissous (mg "*L^{-1}*")")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p_seas

ggsave(plot = p_seas,
       filename = file.path("results", "DO_minmax_seas_slope.png"),
       dpi = 1200,
       units = "cm",
       height = 20,
       width = 28)
# Now want to compare upstream and downstream to get relationships --------
df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))

df_reg <- ungroup(df_DO) %>%
  filter(site == "belleville") %>%
  select(datetime, pos, DO_use) %>%
  pivot_wider(names_from = pos, values_from = DO_use) %>%
  mutate(date = date(datetime),
         hour = hour(datetime)) %>%
  filter(between(date, ymd(19930101), ymd(19930707))) %>%
  group_by(date) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$down~.$up * .$hour)),
         tid = map(mod, broom::tidy)) %>%
  select(-mod, -data) %>%
  unnest(tid)

filter(df_reg, date == ymd(19930605))

ggplot(data = filter(df_reg, term == ".$up", 
                     between(date, ymd(19930601), ymd(19930707))),
       aes(x = date, y = estimate)) +
  geom_point() + geom_line()


ggplot(data = filter(d, date(datetime) == ymd(19930621)) %>%
         pivot_wider(names_from = pos, values_from = DO_use),
       aes(x = down,
           y = up)) +
  geom_point() + 
  geom_path(aes(color = hour(datetime)), size = 2) + 
  scale_color_viridis_c() +
  theme_classic()

x = filter(df_reg, date(datetime) == ymd(19930626)) 

mod <- lm(up~down*hour(datetime), data=x)
summary(mod)
# plot(mod)
y = predict(mod, newdata = x)
plot(y)
lines(x$up)
plot(x$down)
lines(y)
mean(sqrt((y-x$down)^2))

model <- hysteresis::fel(x$up,x$down,period=24,times="equal")
plot(model)
summary(model)
hist((model$residuals))
# Finally, add  back NAs for large chunks of missing data (i.e., >=1 day)
# because the filter can't adequately fill these gaps
df_DO <- df_DO %>%
  group_by(site, year) %>%
  mutate(na_check = 
           {res = rle(is.na(DO_obs));
           rep(res$values * res$lengths, res$lengths)},
         DO_use = ifelse(na_check > 24,
                         NA,
                         filtered)) %>%
  ungroup() %>%
  select(-data, -clean)

# Reduce size
head(df_DO)
df_DO <- select(df_DO, -clean, -data)
# Save data
saveRDS(df_DO, "Data/all_DO_cleaned")




# Wavelets
# x= filter.(df, site == "dampierre", year(datetime) %in% c(1995, 1996), pos == "up") %>%
#   imputeTS::na_kalman()
# my.data.a  = data.frame(date = x$datetime, DO = x$DO_mgL)
# 
# my.w.a <- WaveletComp::analyze.wavelet(my.data.a, "DO",
#                                        loess.span = 0.0, # no detrending required
#                                        dt = 1/(24), # one day has 12*24 5-minute time slots
#                                        dj = 1/24, # resolution along period axis
#                                        lowerPeriod = 1/8, # lowest period of interest: 3 hours
#                                        make.pval = TRUE, # draws white lines indicating significance
#                                        n.sim = 10) # higher number will give smoother white lines
# my.rec.a <- reconstruct(my.w.a, plot.waves = FALSE)
# transactions.rec.a <- my.rec.a$series$transactions.r
# transactions.rec.a[transactions.rec.a < 0] <- 0 # some values are negative
# 
# my.data <- ungroup(df_DO) %>%
#   filter(site == "belleville", year <= 1994) %>%
#   select(datetime, pos, DO_use) %>%
#   pivot_wider(names_from = pos, values_from = DO_use) %>%
#   as.data.frame()
# 
# my.wc <- WaveletComp::analyze.coherency(my.data, my.pair = c("up","down"),
#                            loess.span = 0,
#                            dt = 1/24, dj = 1/100,
#                            lowerPeriod = 1/4,
#                            make.pval = TRUE, n.sim = 10)
# 
# 
# my.rec.b <- WaveletComp::reconstruct(my.wc, my.series = "down",
#                                      plot.waves = FALSE)
# 
# rec_test <- my.rec.b$series



# Plot all sites upstream and downstream
p <- plot_ly(data = pivot_longer(data, cols = -datetime), 
                    x = ~datetime,
                    y = ~value,
                    color = ~name,
                    # linetype = ~pos,
                    colors = c("#1E88E5", "#FFC107", "black")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))

browsable(p)

# select(data4, datetime, DO_mgL, filtered, DO_use) %>% pivot_longer(cols = -datetime), 
# Plot all sites upstream and downstream
p <- plot_ly(data = select(data5, datetime, DO, filtered) %>% pivot_longer(cols = -datetime), 
             x = ~datetime,
             y = ~value,
             color = ~name,
             # linetype = ~pos,
             colors = c("#1E88E5", "#FFC107", "black")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))

htmltools::browsable(p)

x = filter.(data3, between(date(datetime), ymd(19950729), ymd(19970805)))
x$DOstin = na_interpolation(x$DO, option = "stine")
x$DOsea = na_seadec(ts(x$DO, frequency = 24))
x$DOkal = na_kalman(x$DO)
plot_ly(data = select(x, datetime, DO_mgL, filtered, DOsea, DOstin, DOkal) %>% pivot_longer(cols = -datetime), 
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



rm(mod,model)
z = filter.(df, between(date(datetime), ymd(19930101), ymd(20220805)),
            site == "chinon", pos == "down")
# x$DOstin = na_interpolation(x$DO, option = "stine")
# x$DOsea = na_seadec(ts(x$DO, frequency = 24))
# x$DOkal = na_kalman(x$DO)
plot_ly(data = select(ungroup(z), datetime, DO_mgL, DO_use, -pos, -period, -site) %>% 
          pivot_longer(cols = -datetime), 
        x = ~datetime,
        y = ~value,
        color = ~name,
        colors = c("blue", "orange")) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))

d <- ungroup(df_DO) %>%
  filter(site == "chinon") %>%
  select(datetime, pos, DO_use, DO_mgL)

plot_ly(data = filter(d, between(date(datetime), ymd(20000110), ymd(20220707))), 
        x = ~datetime,
        y = ~DO_use,
        color = ~pos) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))



e <- ungroup(df_DO) %>%
  filter(site %in% c("belleville", "dampierre"), pos == "up") %>%
  select(datetime, site, DO_use)

plot_ly(data = filter(e, between(date(datetime), ymd(20100110), ymd(20100707))), 
        x = ~datetime,
        y = ~DO_use,
        color = ~site) %>%
  add_trace(type = "scatter", mode='lines') %>%
  layout(#yaxis2 = list(overlaying = "y", side = "right",
    #             title = TeX("\\text{pH}"),
    #            range = list(6, 11)),
    xaxis = list(title = ""),
    yaxis = list(title = "DO (mg/L)",
                 range = list(0, 30)))

l = filter(df_DO,
           site == "belleville", pos == "up", year == 2010)

plot_ly(data = left_join(my.data,rec_test) %>%
          pivot_longer(cols = -datetime), 
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
LakeMetabolizer::o2.at.sat()