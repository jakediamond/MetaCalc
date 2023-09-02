x <- filter(df, between(month, 6,9))
quantile(x$enh)
# Look at chemical enhancement by trophlux state
ggplot(data = df,
       aes(x = enh)) +
  stat_density(aes(y = after_stat(scaled))) +
  theme_classic() +
  facet_grid(trophic~flux) + 
  scale_x_continuous(limits = c(0, 10)) +
  labs(x = expression(alpha[enh]~"("*`-`*")"),
       y = "scaled density")

df_final <- readRDS(file.path("data", "hourly_data_final.RDS"))
df_com <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete_updateSpCAT.RDS"))
x <- select(df_final, datetime, final_enh = enh,
            final_temp = temp, final_pH = pH, final_KCO2 = KCO2, final_depth = depth) |> 
  left_join(select(df_com, datetime, temp, pH, KCO2, depth, enh))

