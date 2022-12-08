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

# Load functions
source(file.path("src", "000_DO_cleaning_functions.R"))

# Load all DO data---------------------------------------------------------------
# Data directory
df <- readRDS(file.path("data", "01_EDF", "do_timeseries.RDS"))
  
# Plot all sites upstream and downstream
# p_DO_raw <- plot_ly(data = dplyr::filter(df, site != "stlaurent"), 
#                     x = ~datetime,
#                     y = ~DO_mgL,
#                     color = ~site,
#                     linetype = ~pos,
#                     colors = c("#1E88E5", "#FFC107", "black", "#D55E00")) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#           #             title = TeX("\\text{pH}"),
#            #            range = list(6, 11)),
#          xaxis = list(title = ""),
#          yaxis = list(title = "DO (mg/L)",
#                       range = list(0, 30)),
#          title = "raw DO data")

# browsable(p_DO_raw)

# Nest data by site and period
# Split into 5 year periods
df_n <- df %>%
  distinct() %>%
  select(-QC) %>%
  mutate(yr = year(datetime) - min(year(datetime))) %>%
  group_by(site, pos) %>%
  nest() %>%
  mutate(period = map(data, ~as.numeric(cut_interval(.$yr, length = 5)))) %>%
  unnest(c(data, period)) %>%
  group_by(site, pos, period) %>%
  nest()

# Remove original data to keep workspace clean
rm(df)

# Apply functions ---------------------------------------------------------
# Clean the data and use the lowpass filter
df_n <- df_n %>%
  filter.(site != "stlaurent") %>% #not enough data at st laurent to care
  mutate(clean = map(data, clean_fun))

# Get all in one dataframe
df_DO <- unnest(df_n, clean)

# do a kalman imputation, max 12 hours
df_DO <- df_DO %>%
  group_by(site, pos) %>%
  arrange(datetime) %>%
  mutate(DO_use = na_kalman(DO, maxgap = 12))

# Plot all sites upstream and downstream
p_DO_clean <- plot_ly(data = dplyr::filter(df_DO, site == "belleville"), 
                    x = ~datetime,
                    y = ~DO_use,
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
# 
# 
# x = filter(df_DO, site == "belleville", pos == "up", between(datetime,
#                                                              ymd(20050701),
#                                                              ymd(20110814)))
# 
# xts = ts(x$DO_use, frequency = 24)
# xcl = na_seasplit(xts, "kalman")
# plot(xcl)
# 
# # A second cleaning -------------------------------------------------------
# df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS")) %>%
#   select(site, pos, year, datetime, period, DO_mgL = DO_use, namax, narun, DO)
# 
# # Remove entire years of missing data (more than a half year)
# # Remove DO data again when NA is longer than 3 days, initial na_kalman didn't work
# df_DO <- df_DO %>%
#   group_by(site, pos, year) %>%
#   filter.(!(sum(is.na(DO_mgL)) > 180 * 24)) %>%
#   ungroup() %>%
#   group_by(site, year, pos, narun) %>%
#   mutate(DO_mgL = if_else((namax > 24 * 3)  & is.na(DO), NA_real_, DO_mgL)) %>%
#   ungroup() %>%
#   select(-namax, -narun, -DO)
# 
# # Nest data and apply
# df_n <- df_DO %>%
#   group_by(site, pos, period) %>%
#   nest() %>%
#   mutate(clean = map(data, possibly(clean_fun, NA_real_)))
# 
# # Get all in one dataframe
# df_clean <- filter.(df_n, !is.na(clean)) %>%
#   unnest(clean)
# 
# x = pluck(df_n, 4,2)
# 
# 
# 
# # Plot all daily min and max for trends -----------------------------------
# df_hourly_an <- read_csv(file.path("data", "01_EDF", "hourly EDF data", 
#                                 "hourly_DO.csv")) %>% 
#   rename(belleville = BEL, dampierre = DAM, chinon = CHB) %>%
#   mutate(datetime = mdy_hm(datetime))
# 
# # Load my data
# df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))
# 
# # Get both together
# df_DO_all <- df_hourly_an %>%
#   pivot_longer(cols = c(belleville, dampierre, chinon), names_to = "site") %>%
#   mutate(pos = "up") %>%
#   rename(DO_use = value) %>%
#   anti_join(df_DO, by = c("datetime", "site", "pos")) %>%
#   bind_rows(df_DO) %>%
#   mutate(date = date(datetime))
# 
# df_minmax <- df_DO_all %>%
#   group_by(site, date, pos) %>%
#   summarize(max = max(DO_use, na.rm = T),
#             min = min(DO_use, na.rm = T)) %>%
#   pivot_longer(cols = c(min, max))
# 
# # Get Chinon min max
# df_chi <- readxl::read_xlsx(file.path("data", "01_EDF", "Chinon_1993_2005_FM.xlsx"),
#                             sheet = 4) %>%
#   rename(date = DATE, min = `Mini Amont`, max = `Max amont`) %>%
#   pivot_longer(-date) %>%
#   mutate(site = "chinon", pos = "up", date = as.Date(date))
# 
# # Replace this in the data
# df_minmax <- filter(df_minmax, !(site == "chinon" & pos == "up" & 
#                                  between(year(date), 1993, 2005))) %>%
#   bind_rows(df_chi) %>%
#   mutate(moyr = ymd(paste(year(date), month(date), "01")),
#          site_f = factor(site)) %>%
#   mutate(site_f = fct_relevel(site_f, "belleville", "dampierre", "civaux", "chinon"))
# 
# # Get temperature
# df_t <- readRDS(file.path("data", "05_hourly_data_clean", 
#                           "temp_discharge_rad_data.RDS")) %>%
#   group_by(site, date) %>%
#   summarize(max = max(temp_C, na.rm = T),
#             min = min(temp_C, na.rm = T)) %>%
#   pivot_longer(cols = c(min, max), values_to = "temp")
# 
# # Add temp and calculate %sat
# df_minmax <- left_join(df_minmax, df_t) %>%
#   mutate(DO_sat = ifelse(temp <= 0,
#                          0,
#                          14.652 - 0.41022 * temp + 0.007991 * 
#                            temp^2 - 0.000077774 * temp^3),
#          DO_per = value / DO_sat)
#   
# 
# p_minmax <- ggplot(filter(df_minmax, pos == "up"),
#        aes(x = moyr,
#            y = DO_per * 100,
#            color = name)) +
#   stat_summary(fun = median, geom = "line", size = 1.2) +
#   scale_color_manual(values = c("#0072B2", "#D55E00")) +
#   # stat_smooth(method = "loess", span = 0.1, color = "black") +
#   facet_grid(rows = vars(site_f)) +
#   theme_classic(base_size = 18) +
#   labs(x = "",
#        y = expression("oxygène dissous (% sat.)")) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "none")
# p_minmax
# 
# ggsave(plot = p_minmax,
#        filename = file.path("results", "DO_minmax_amont_sat.png"),
#        dpi = 1200,
#        units = "cm",
#        height = 20,
#        width = 28)
# 
# # Seasonality
# df_minmax_seas <- df_minmax %>%
#   mutate(month = month(date),
#          season = case_when(month %in% c(3,4,5) ~ "spring",
#                             month %in% c(6,7,8) ~ "summer",
#                             month %in% c(9,10,11) ~ "autumn",
#                             TRUE ~ "winter"),
#          year = year(date)) %>%
#   mutate(seaf = factor(season)) %>%
#   mutate(seaf = fct_relevel(seaf, "spring", "summer", "autumn", "winter")) %>%
#   group_by(year, season, name, site_f, seaf, pos) %>%
#   summarize(val = median(value, na.rm = T),
#             sat = median(DO_per * 100, na.rm = T))
# 
# p_seas <- ggplot(data = filter(df_minmax_seas, pos == "up",
#                                !(site_f == "belleville" & year == 1995)),
#                  aes(x = year,
#                      y = sat,
#                      color = name)) +
#   geom_point() + geom_line() +
#   # stat_summary() + stat_summary(geom = "line") +
#   scale_color_manual(values = c("#0072B2", "#D55E00")) +
#   stat_smooth(method = "lm") +
#   ggpubr::stat_regline_equation(formula=y~x) +
#   ggpubr::stat_cor(aes(label = ..p.label..), label.x = 2010) +
#   facet_grid(rows = vars(site_f), cols = vars(seaf)) +
#   theme_classic(base_size = 18) +
#   labs(x = "",
#        y = expression("oxygène dissous (% sat.)")) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "none")
# 
# p_seas
# 
# ggsave(plot = p_seas,
#        filename = file.path("results", "DO_minmax_seas_slope_amont_sat.png"),
#        dpi = 1200,
#        units = "cm",
#        height = 20,
#        width = 28)
# # Now want to compare upstream and downstream to get relationships --------
# df_DO <- readRDS(file.path("data", "05_hourly_data_clean", "DO_cleaned_part1.RDS"))
# 
# df_reg <- ungroup(df_DO) %>%
#   filter(site == "belleville") %>%
#   select(datetime, pos, DO_use) %>%
#   pivot_wider(names_from = pos, values_from = DO_use) %>%
#   mutate(date = date(datetime),
#          hour = hour(datetime)) %>%
#   filter(between(date, ymd(19930101), ymd(19930707))) %>%
#   group_by(date) %>%
#   nest() %>%
#   mutate(mod = map(data, ~lm(.$down~.$up * .$hour)),
#          tid = map(mod, broom::tidy)) %>%
#   select(-mod, -data) %>%
#   unnest(tid)
# 
# filter(df_reg, date == ymd(19930605))
# 
# ggplot(data = filter(df_reg, term == ".$up", 
#                      between(date, ymd(19930601), ymd(19930707))),
#        aes(x = date, y = estimate)) +
#   geom_point() + geom_line()
# 
# 
# ggplot(data = filter(d, date(datetime) == ymd(19930621)) %>%
#          pivot_wider(names_from = pos, values_from = DO_use),
#        aes(x = down,
#            y = up)) +
#   geom_point() + 
#   geom_path(aes(color = hour(datetime)), size = 2) + 
#   scale_color_viridis_c() +
#   theme_classic()
# 
# x = filter(df_reg, date(datetime) == ymd(19930626)) 
# 
# mod <- lm(up~down*hour(datetime), data=x)
# summary(mod)
# # plot(mod)
# y = predict(mod, newdata = x)
# plot(y)
# lines(x$up)
# plot(x$down)
# lines(y)
# mean(sqrt((y-x$down)^2))
# 
# model <- hysteresis::fel(x$up,x$down,period=24,times="equal")
# plot(model)
# summary(model)
# hist((model$residuals))
# # Finally, add  back NAs for large chunks of missing data (i.e., >=1 day)
# # because the filter can't adequately fill these gaps
# df_DO <- df_DO %>%
#   group_by(site, year) %>%
#   mutate(na_check = 
#            {res = rle(is.na(DO_obs));
#            rep(res$values * res$lengths, res$lengths)},
#          DO_use = ifelse(na_check > 24,
#                          NA,
#                          filtered)) %>%
#   ungroup() %>%
#   select(-data, -clean)
# 
# # Reduce size
# head(df_DO)
# df_DO <- select(df_DO, -clean, -data)
# # Save data
# saveRDS(df_DO, "Data/all_DO_cleaned")
# 
# 
# 
# 
# # Wavelets
# # x= filter.(df, site == "dampierre", year(datetime) %in% c(1995, 1996), pos == "up") %>%
# #   imputeTS::na_kalman()
# # my.data.a  = data.frame(date = x$datetime, DO = x$DO_mgL)
# # 
# # my.w.a <- WaveletComp::analyze.wavelet(my.data.a, "DO",
# #                                        loess.span = 0.0, # no detrending required
# #                                        dt = 1/(24), # one day has 12*24 5-minute time slots
# #                                        dj = 1/24, # resolution along period axis
# #                                        lowerPeriod = 1/8, # lowest period of interest: 3 hours
# #                                        make.pval = TRUE, # draws white lines indicating significance
# #                                        n.sim = 10) # higher number will give smoother white lines
# # my.rec.a <- reconstruct(my.w.a, plot.waves = FALSE)
# # transactions.rec.a <- my.rec.a$series$transactions.r
# # transactions.rec.a[transactions.rec.a < 0] <- 0 # some values are negative
# # 
# # my.data <- ungroup(df_DO) %>%
# #   filter(site == "belleville", year <= 1994) %>%
# #   select(datetime, pos, DO_use) %>%
# #   pivot_wider(names_from = pos, values_from = DO_use) %>%
# #   as.data.frame()
# # 
# # my.wc <- WaveletComp::analyze.coherency(my.data, my.pair = c("up","down"),
# #                            loess.span = 0,
# #                            dt = 1/24, dj = 1/100,
# #                            lowerPeriod = 1/4,
# #                            make.pval = TRUE, n.sim = 10)
# # 
# # 
# # my.rec.b <- WaveletComp::reconstruct(my.wc, my.series = "down",
# #                                      plot.waves = FALSE)
# # 
# # rec_test <- my.rec.b$series
# 
# 
# 
# # Plot all sites upstream and downstream
# p <- plot_ly(data = pivot_longer(data, cols = -datetime), 
#                     x = ~datetime,
#                     y = ~value,
#                     color = ~name,
#                     # linetype = ~pos,
#                     colors = c("#1E88E5", "#FFC107", "black")) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# browsable(p)
# 
# # select(data4, datetime, DO_mgL, filtered, DO_use) %>% pivot_longer(cols = -datetime), 
# # Plot all sites upstream and downstream
# p <- plot_ly(data = select(data5, datetime, DO, filtered) %>% pivot_longer(cols = -datetime), 
#              x = ~datetime,
#              y = ~value,
#              color = ~name,
#              # linetype = ~pos,
#              colors = c("#1E88E5", "#FFC107", "black")) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# htmltools::browsable(p)
# 
# x = filter.(data3, between(date(datetime), ymd(19950729), ymd(19970805)))
# x$DOstin = na_interpolation(x$DO, option = "stine")
# x$DOsea = na_seadec(ts(x$DO, frequency = 24))
# x$DOkal = na_kalman(x$DO)
# plot_ly(data = select(x, datetime, DO_mgL, filtered, DOsea, DOstin, DOkal) %>% pivot_longer(cols = -datetime), 
#         x = ~datetime,
#         y = ~value,
#         color = ~name) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# 
# 
# rm(mod,model)
# z = filter.(df, between(date(datetime), ymd(19930101), ymd(20220805)),
#             site == "chinon", pos == "down")
# # x$DOstin = na_interpolation(x$DO, option = "stine")
# # x$DOsea = na_seadec(ts(x$DO, frequency = 24))
# # x$DOkal = na_kalman(x$DO)
# plot_ly(data = select(ungroup(z), datetime, DO_mgL, DO_use, -pos, -period, -site) %>% 
#           pivot_longer(cols = -datetime), 
#         x = ~datetime,
#         y = ~value,
#         color = ~name,
#         colors = c("blue", "orange")) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# d <- ungroup(df_DO) %>%
#   filter(site == "chinon") %>%
#   select(datetime, pos, DO_use, DO_mgL)
# 
# plot_ly(data = filter(d, between(date(datetime), ymd(20000110), ymd(20220707))), 
#         x = ~datetime,
#         y = ~DO_use,
#         color = ~pos) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# 
# 
# e <- ungroup(df_DO) %>%
#   filter(site %in% c("belleville", "dampierre"), pos == "up") %>%
#   select(datetime, site, DO_use)
# 
# plot_ly(data = filter(e, between(date(datetime), ymd(20100110), ymd(20100707))), 
#         x = ~datetime,
#         y = ~DO_use,
#         color = ~site) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# 
# l = filter(df_DO,
#            site == "belleville", pos == "up", year == 2010)
# 
# plot_ly(data = left_join(my.data,rec_test) %>%
#           pivot_longer(cols = -datetime), 
#         x = ~datetime,
#         y = ~value,
#         color = ~name) %>%
#   add_trace(type = "scatter", mode='lines') %>%
#   layout(#yaxis2 = list(overlaying = "y", side = "right",
#     #             title = TeX("\\text{pH}"),
#     #            range = list(6, 11)),
#     xaxis = list(title = ""),
#     yaxis = list(title = "DO (mg/L)",
#                  range = list(0, 30)))
# LakeMetabolizer::o2.at.sat()