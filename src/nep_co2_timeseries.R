
# Plot an interactive time series -----------------------------------------
p_ts <- plot_ly(data = df_clean, 
                x=~date) %>%
  add_trace(y= ~ -filtered_NEP_mean, type = "scatter", mode='lines', 
            color = I("#1E88E5"), showlegend = FALSE) %>%
  add_trace(y= ~ filtered_FCO2_mean, type = "scatter", mode='lines', 
            color = I("#FFC107"), showlegend = FALSE) %>%
  layout(yaxis2 = list(overlaying = "y", side = "right",
                       title = TeX("\\text{pH}")),
         xaxis = list(title = ""),
         yaxis = list(title = TeX("\\color{#1E88E5}{NEP}~or~\\color{#FFC107}{FCO_{2}}~(mmol~m^{-3})")),
         title = TeX("\\text{carbonate system}"),
         margin = list(r = 50)) %>%
  config(mathjax = "cdn")

htmltools::browsable(p_ts)

# try a monthly time series to show the data better
# Nice plotting format first
df_l <- df_clean %>%
  select(date, NEP_mean = filtered_NEP_mean, 
         NEP_2.5 = filtered_NEP_2.5, NEP_97.5 = filtered_NEP_97.5,
         CO2_mean = filtered_CO2_meanenh,
         CO2_2.5 = filtered_CO2_2.5enh, CO2_97.5 = filtered_CO2_97.5enh) %>%
  mutate(year = year(date),
         month = month(date),
         across(contains("NEP"), ~.*1)) %>%
  pivot_longer(cols = -c(month, year, date), names_sep = "_", 
               names_to = c("type", "val_type")) %>%
  pivot_wider(names_from = val_type, values_from = value) %>%
  mutate(wts = 1/abs(`2.5` - `97.5`)^2) %>%
  group_by(type) %>%
  arrange(type, date) %>%
  mutate(rollmean = slider::slide_dbl(
    .x = cur_data(),
    .f = ~ weighted.mean(
      x = .x$mean,
      w = .x$wts
    ), .before = 90
  ))

# Weighted average time series
p_ts_mo <- ggplot(data = df_l,
                  aes(x = date,
                      color = type)) +
  geom_line(aes(y = rollmean),  linewidth = 1.2) + 
  scale_color_manual(name = "", values = c("#1E88E5", "#FFC107")) +
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  labs(x = "",
       y = expression(CO[2]~flux~"("*mmol~m^{-2}~d^{-1}*")"))

p_ts_mo

# Same thing for the ratio, but need to calculate with uncertainty first
df_rat <- df_clean %>%
  mutate(dNEP_mean = (filtered_NEP_mean - filtered_NEP_2.5) / 1.96,
         dCO2_mean = (filtered_CO2_meanenh - filtered_CO2_2.5enh) / 1.96) %>%
  select(date, NEP_mean = filtered_NEP_mean, dNEP_mean,
         CO2_mean = filtered_CO2_meanenh, dCO2_mean) %>%
  mutate_with_error(ratio ~ -NEP_mean/ CO2_mean) %>%
  mutate(wts = 1/dratio^2) %>%
  mutate(rollmean = slider::slide_dbl(
    .x = .,
    .f = ~ weighted.mean(
      x = .x$ratio,
      w = .x$wts
    ), .before = 365
  ))

p_ts_mo_rat <- ggplot(data = df_rat,
                      aes(x = date, y =rollmean)) +
  geom_line(aes(y = rollmean),  linewidth = 1.2) + 
  scale_color_manual(name = "", values = c("#1E88E5", "#FFC107")) +
  theme_classic() + 
  # scale_y_continuous(limits = c(0, 5)) +
  geom_hline(yintercept = 1) + 
  labs(x = "",
       y = expression(NEP~":"*CO[2]~"(-)"))

p_ts_mo_rat

p_ts <- p_ts_mo / p_ts_mo_rat + plot_annotation(tag_levels = "a")
ggsave(plot= p_ts,
       filename = file.path("results", "co2_nep_timeseries.png"),
       dpi = 600,
       units = "cm",
       height = 16,
       width = 14.8)