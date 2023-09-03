df <- readRDS(file.path("data", "hourly_data.RDS"))

df |>
  summarize(sum((CO2_uM - CO2eq_uM) < 0, na.rm = T) / n(),
            sum(exO2 > 0, na.rm = T) / n()) 

# Median CO2
df |>
  summarize(co2 = median(CO2_uM, na.rm = T),
            pco2 = median(pCO2_uatm, na.rm = T))

df <- readRDS(file.path("data", "daily_trophlux.RDS"))

df |>
  rename(FCO2 = filtered_FCO2_mean, NEP = filtered_NEP_mean) |>
  # drop_na(FCO2, NEP) |>
  filter(FCO2 > 0, NEP < 0) |>
  mutate(rat = -NEP / FCO2 * 100) |>
  summarize(mean_nep = mean(NEP, na.rm = T),
            med_nep = median(NEP, na.rm = T),
            mean_rat = mean(rat, na.rm = T),
            # wt_rat = weighted.mean(rat, )
            q_rat = quantile(rat, na.rm = T))




df |>
  drop_na(FCO2) |>
  summarize(
    K600_m = mean(FCO2),
    K600_med = median(FCO2),
    K600_Wt = weighted.mean(FCO2, (1/(FCO2*0.3)^2)),
    K600_kwt = weighted.mean(FCO2, (KCO2*depth))
  )


# Daily means
df_meand <- df |>
  group_by(date) |>
  summarize(across(where(is.numeric), mean)) |>
  left_join(select(df_nepco2, date, co2_old = filtered_CO2_meanenh,
                   nep_old = filtered_NEP_mean))

# daily sums
df_sumd <- df |>
  group_by(date) |>
  select(date, NEP, FCO2) |>
  summarize(across(where(is.numeric), sum)) |>
  left_join(select(df_nepco2, date, year, month, co2_old = filtered_CO2_meanenh,
                   nep_old = filtered_NEP_mean))

# estimate CO2 based on daily average AT, pH, SpC
df_d <- df_meand |>
  mutate(CO2_avg = DIC(temp+273.15, AT, pH, SpC)$CO2,
         FCO2_avg = enh * KCO2 * depth * (CO2_avg-CO2eq_uM)) |>
  select(date, FCO2_avg) |>
  right_join(df_sumd)

# Comparison data frame
df_comp <- select(df_d, date, FCO2_avg, FCO2_sum = FCO2 , FCO2_old = co2_old) |>
  # right_join(df) |>
  mutate(FCO2_sum = FCO2_sum / 24) |>
  drop_na(FCO2_avg, FCO2_sum)

# RMSE and % diff from sum of hourly and mean daily
sqrt(mean((df_comp$FCO2_sum - df_comp$FCO2_avg)^2))
filter(df_comp, FCO2_sum < 0, FCO2_avg < 0) |>
  mutate(perdif = (FCO2_avg - FCO2_sum) / FCO2_sum * 100) |>
  summarize(pdf = mean(perdif))
mean((df_comp$FCO2_avg - df_comp$FCO2_sum) / df_comp$FCO2_sum)*100



# Test if CO2 is the same with seacarb/CO2sys
seacarb::carb(flag = 8,
              var1 = 8.26,
              var2 = 1.91 / rhow(22.6),
              T = 22.6,
              S = 1.6*10^-5 * 267.5 * 59,
              pH = "F",
              k1k2 = "m06"
              )




# Overall mean FCO2 and NEP
df_sumd |>
  drop_na(FCO2, NEP) |>
  summzarize

# Overall internal CO2 production
df_sumd |>
  drop_na(FCO2, NEP) |>
  # filter(FCO2 > 0, NEP < 0) |>
  mutate(rat = nep_old / (FCO2 / 24) * 100) |>
  summarize(mean_nep = mean(-nep_old, na.rm = T),
            med_nep = median(-nep_old, na.rm = T),
            mean_rat = mean(rat, na.rm = T),
            wt_rat = weighted.mean(rat, )
            q_rat = quantile(rat, na.rm = T))

# Yearly FCO2
df_y <- ungroup(df_sumd) |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  drop_na(FCO2, NEP) |>
  filter(FCO2 > 0, NEP < 0) |>
  summarize(FCO2_y = sum(FCO2/24, na.rm = T),
            FCO2_old = sum(co2_old, na.rm = T),
            NEP_y = sum(NEP, na.rm = T),
            NEP_old = sum(-nep_old, na.rm = T)) |>
  mutate(rat = round(-NEP_old / FCO2_y * 100, 2))

# Look at distribution of trophlux by year
ungroup(df_sumd) |>
  mutate(wy = if_else(month(date) > 8, year + 1, year)) |>
  group_by(wy) |>
  mutate(totFCO2 = sum(FCO2/24, na.rm = T),
         totNEP = sum(-nep_old, na.rm = T),
         len = n()) |>
  ungroup() |>
  mutate(regime = if_else(year(date) <2005, "auto", "hetero"),
         troph = if_else(-nep_old > 0, "auto", "hetero"),
         flux = if_else(FCO2 > 0, "1source", "2sink")) |>
  group_by(wy, regime, troph, flux) |>
  drop_na(FCO2, nep_old) |>
  # filter(FCO2 > 0, NEP < 0) |>
  summarize(freq = n() / len) |>
  distinct() |>
  ungroup() |>
  group_by(regime, troph, flux) |>
  summarize(f = median(freq))



ungroup(df_sumd) |>
  mutate(regime = if_else(year(date) <2005, "auto", "hetero"),
         troph = if_else(-nep_old > 0, "auto", "hetero"),
         flux = if_else(FCO2 > 0, "1source", "2sink")) |>
  group_by(regime) |>
  drop_na(FCO2, nep_old) |>
  mutate(totFCO2 = sum(FCO2/24),
         totNEP = sum(-nep_old),
         len = n()) |>
  ungroup() |>
  group_by(regime, troph, flux) |>
  summarize(count = n(),
            fr = count / mean(len),
            frac = sum(FCO2/24) / mean(totFCO2),
            nep = median(-nep_old),
            fco2 = median(FCO2/24),
            rat = median(nep_old/(FCO2/24)))
 
# median(df_y$rat)
# # yearly           
# x <- df |>
#   drop_na(pH, AT, temp, SpC) |>
#   mutate(pCO2carb = seacarb::carb(flag= 8,
#                                  var1 = pH,
#                                  var2 = AT / rhow(temp),
#                                  T = temp,
#                                  S = 1.6E-5* 53.974* SpC)$pCO2)#CO2*rhow(temp)*1000)
# y <- median(x$CO2carb/x$CO2_uM)
# z <- filter(x, between(datetime, ymd_h(2004092013), ymd_h(2004092113)))
# ee <- dd |>
#   select(date, datetime, FCO2_inst = FCO2, FCO2_avg, FCO2_sum, FCO2_d = co2_d) |>
#   pivot_longer(cols = contains("FCO2"))
# 
# plot_ly(data = ee,
#         x=~datetime,
#         y = ~value,
#         color = ~name) |>
#   add_trace(type = "scatter", mode='lines') |>
#   layout(xaxis = list(title = ""),
#          yaxis = list(title = TeX("\\color{#1E88E5}FCO2~(mM~m^{-2})"))
#          # title = TeX("\\text{Conductivity upstream downstream}"))
#   ) |>
#   config(mathjax = "cdn")


# df_met_err |>
#   drop_na(NEP_mean) |>
#   summarize(
#     GPP_m = mean(NEP_mean) * 1000/32,
#     GPP_med = median(NEP_mean) * 1000/32,
#     GPP_Wt = weighted.mean(NEP_mean, (1/dNEP_mean^2)) * 1000/32,
#     GPP_kwt = weighted.mean(NEP_mean, K600_mean) *1000/32
#   )