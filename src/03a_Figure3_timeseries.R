# -------------------------------------
# Author: Jake Diamond
# Purpose: plot Figure 3 (or 2, depending) of archetypal timeseries for trophlux states
# Date: 2023-05-07
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS"))
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS"))

# Get hourly averages of everything 
df_avg <- df |>
  mutate(exCO2 = CO2_uM - CO2eq_uM,
         calc_mMm2hr = -0.5 * (AT_mM - lag(AT_mM)) * depth_m * 1000) |>
  mutate(
    calc_mMm2hr = if_else(SI < 0.4 & calc_mMm2hr > 0, 0, calc_mMm2hr),
    # calcification shouldn't be positive if SI is < 0
    calc_mMm2hr = if_else(SI < 0 & calc_mMm2hr > 0, 0, calc_mMm2hr),
    # And the inverse
    calc_mMm2hr = if_else(SI > 0 & calc_mMm2hr < 0, 0, calc_mMm2hr)
  ) |>
  left_join(select(df_tf, date, troph, sourcesink, trophlux)) |>
  pivot_longer(cols = -c(year, month, troph, sourcesink, trophlux, date, 
                         datetime, hr, depth_m, Q_m3s)) |>
  group_by(trophlux, name, hr) |>
  summarize(mean_val = mean(value, na.rm = T),
            se_val = sd(value, na.rm = T) / sqrt(n()) * 1.96) |>
  ungroup() |>
  mutate(trophlux = str_replace(trophlux, "_", " "))

# x <- filter(df, trophlux == "autotrophic_source", exCO2 < 0)
# Overall plot theme ------------------------------------------------------
th <- list(geom_line(aes(y = mean_val,
                         color = trophlux)),
           geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                           fill = trophlux),
                       alpha = 0.4),
           theme_classic(base_size = 10),
           theme(legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major.x = element_line(),
                 axis.ticks.y = element_line(color = "black"),
                 axis.line.x = element_blank(),
                 panel.border = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()),
           scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")),
           scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")),
           guides(fill = "none"))


# Layout of three panels vertically
layout <- c(area(t = 1, l = 1, b = 3, r = 5),
            area(t = 4, l = 1, b = 6, r = 5),
            area(t = 7, l = 1, b = 9, r = 5))
# Figures -----------------------------------------------------------------
p_o2ex <- ggplot(data = filter(df_avg, name == "exO2") |>
                   drop_na(trophlux),
                 aes(x = hr)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  scale_y_continuous(limits = c(-110, 200)) +
  th +
  labs(y = expression(exO[2]~"("*mu*M*")"))
p_o2ex

# Specific conductivity
p_C <- ggplot(data = filter(df_avg, name == "SpC") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  th +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(250, 290)) + 
  labs(y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_C

# pH
p_pH <- ggplot(data = filter(df_avg, name == "pH") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  th +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(7.5, 9.5)) +
  theme(axis.line.x = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) +
  labs(y = "pH",
       x = "solar hour")
p_pH

# Whole plot of measurements
p_meas <- p_o2ex + p_C + p_pH + plot_layout(design = layout)

# Plot carbonate system ---------------------------------------------------
p_co2ex <- ggplot(data = filter(df_avg, name == "exCO2") |>
                   drop_na(trophlux),
                 aes(x = hr)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  th +
  labs(y = expression(exCO[2]~"("*mu*M*")"))
p_co2ex

# bicarbonate
p_hco3 <- ggplot(data = filter(df_avg, name == "HCO3_uM") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  th + 
  labs(y = expression(HCO[3]^{`-`}~"("*mu*M*")"))
p_hco3

# calcite saturation index
p_SI <- ggplot(data = filter(df_avg, name == "SI") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  th + 
  theme(axis.line.x = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) +
  labs(y = "SI (-)",
       x = "solar hour")
p_SI
p_carb <- p_co2ex + p_hco3 + p_SI + plot_layout(design = layout)

# Plot Fluxes ------------------------------------------------------------------
p_nep <- ggplot(data = filter(df_avg, name == "NEP") |>
                    drop_na(trophlux),
                  aes(x = hr)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  th +
  labs(y = expression(NEP~"("*mmol~m^{-2}~h^{-1}*")"))
p_nep

# Co2 flux
p_fco2 <- ggplot(data = filter(df_avg, name == "FCO2") |>
                   drop_na(trophlux),
                 aes(x = hr)) +
  geom_hline(yintercept = 0) + 
  th +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  labs(y = expression(FCO[2]~"("*mmol~m^{-2}~h^{-1}*")"))
p_fco2

# Calcite precipitation
p_calc <- ggplot(data = filter(df_avg, name == "calc_mMm2hr") |>
                  drop_na(trophlux),
                aes(x = hr)) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  th +
  geom_hline(yintercept = 0) +
  theme(axis.line.x = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) +
  labs(y = expression(CaCO[3]~ppt.~"("*mmol~m^{-2}~h^{-1}*")"),
       x = "solar hour")
p_calc
p_flux <- p_nep + p_fco2 + p_calc + plot_layout(design = layout)

# Overall plot
ptot <- (p_meas | p_carb | p_flux) + plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") & 
  guides(color = "none", fill = guide_legend("trophlux state",
                                             override.aes = list(alpha=1))) &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.key.size = unit(0.25, "cm"),
        plot.tag.position = c(0.27, 0.95))
# ptot
ggsave(plot = ptot,
       filename = file.path("results", "figure3_timeseries_v2.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 15.2)
