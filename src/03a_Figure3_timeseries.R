# -------------------------------------
# Author: Jake Diamond
# Purpose: plot Figure 3 (or 2, depending) of archetypal timeseries for trophlux states
# Date: 2023-05-07
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)
source(file.path("src", "000_error_propagation_function.R"))

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS"))
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS"))
y <- filter(df, date == ymd(20100719)) # yes
z <- filter(df, date == ymd(20150803)) #maybe
a <- filter(df, date == ymd(20180917))
# Take a look at the hourly data to get an idea of how to plot it
df %>%
  left_join(select(df_tf, date, troph, sourcesink, trophlux)) |>
  filter(trophlux == "heterotrophic source") %>%
  ggplot(aes(x = hr,
             y = exO2)) +
  geom_line(aes(group = date), alpha = 0.01) +
  stat_summary(aes(group = hr)) +
  theme_classic()

# Get hourly averages of everything 
df_avg <- df |>
  left_join(select(df_tf, date, troph, sourcesink, trophlux)) |>
  pivot_longer(cols = -c(year, month, troph, sourcesink, trophlux, date, 
                         datetime, hr, depth_m, Q_m3s)) |>
  group_by(trophlux, name, hr) %>%
  summarize(mean_val = mean(value, na.rm = T),
            # sd_val = if_else(grepl("d", name), sqrt(sum(value^2)), sd(value))) %>%
            se_val = sd(value, na.rm = T) / sqrt(n()),
            q25 = quantile(value, 0.25, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T)) |>
  ungroup() |>
  mutate(trophlux = str_replace(trophlux, "_", " "))

# Do error propagation for variables
df_errprop <- df |>
  # Forgot to initially propagate error to exCO2 and hourly NEP, SI
  mutate(dCO2eq_uM = 0,
         dO2_uM = if_else(year < 2008, 0.04*O2_uM, 0.01*O2_uM), # = 0.4% and 1% uncertainty
         dpH = if_else(year < 2008, 0.1, 0.05),
         dKO2_smooth = dK600_mean,
         ddepth_smooth = 0,
         dO2sat_uM = 0,
         lagO2 = lag(O2_uM),
         dlagO2 = lag(dO2_uM),
         dFCO2 = dFCO2 * sqrt(10000), # these were saved as standard errors of 10000 run monte carlo
         dFCO2_enh = dFCO2 * sqrt(10000),
         gamma = act(temp+273, SpC),
         act = 10^(4*-gamma),
         Ca = Ca_mM * act / 1000,
         CO3 = CO3_uM * act/ 1E6,
         dCO3 = 0.08 * CO3_uM / 1E6, #error is 8% for the carb system, based on CO2 and DIC
         dCa = 0.07 / 1E3, # from RMSE of random forest model
         ksp = Ksp(temp+273, SpC * 1.6E-5*53.974),
         dksp = 0,
         dcalc_mMm2hr = 133E-3 #from the Alkalinity model
         ) %>%
  # Propoagate error to exCO2, exO2, hourly NEP, and SI
  mutate_with_error(exCO2 ~ CO2_uM - CO2eq_uM) %>%
  mutate_with_error(exO2 ~ O2sat_uM - O2_uM) %>%
  mutate_with_error(NEP ~ ((O2_uM - lagO2) - (KO2_smooth / 24) * 
           (O2sat_uM - O2_uM)) * depth_smooth) %>%
  mutate_with_error(SI ~ log10(Ca * CO3 / ksp)) %>%
  select(date, datetime, hr, dCO2_uM, dDIC_uM, dexCO2, dFCO2_enh, 
         dNEP, dFCO2, dexO2, dSI, dcalc_mMm2hr) %>%
  left_join(select(df_tf, date, troph, sourcesink, trophlux)) %>%
  pivot_longer(cols = -c(troph, sourcesink, trophlux, date, 
                         datetime, hr)) |>
  group_by(trophlux, name, hr) %>%
  # Get the standard error of the errors across all measurements, gaussian
  summarize(se_val = sqrt(sum(value^2, na.rm = T) / n())) %>%
  ungroup() %>%
  mutate(name = str_remove(name, "d"))

df_avg <- df_avg %>%
  rows_update(df_errprop, by = c("trophlux", "name", "hr"))
# Overall plot theme ------------------------------------------------------
th <- list(geom_line(aes(y = mean_val,
                         color = trophlux),
                     linewidth = 1.5),
           geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                           fill = trophlux),
                       alpha = 0.1),
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
# For observations
th_ob <- list(geom_line(aes(y = mean_val,
                         color = trophlux),
                     linewidth = 1.5),
           geom_ribbon(aes(ymin = q25, ymax = q75,
                           fill = trophlux),
                       alpha = 0.1),
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
# Observed data panels, show with 25-75% bands
p_o2ex <- ggplot(data = filter(df_avg, name == "exO2") |>
                   drop_na(trophlux),
                 aes(x = hr)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  scale_y_continuous(limits = c(-200, 350)) +
  th_ob +
  labs(y = expression(exO[2]~"("*mu*M*")"))
p_o2ex

# Specific conductivity
p_C <- ggplot(data = filter(df_avg, name == "SpC") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  th_ob +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  # scale_y_continuous(limits = c(250, 290)) + 
  labs(y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_C

# pH
p_pH <- ggplot(data = filter(df_avg, name == "pH") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  th_ob +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(7., 9.5)) +
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
  # scale_y_continuous(limits = c(-25, 100), breaks = seq(-25, 100, 25)) +
  th +
  labs(y = expression(exCO[2]~"("*mu*M*")"))
p_co2ex

# bicarbonate
p_DIC <- ggplot(data = filter(df_avg, name == "DIC_uM") |>
                 drop_na(trophlux),
               aes(x = hr)) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  th + 
  labs(y = expression(DIC~"("*mu*M*")")) #HCO[3]^{`-`}
p_DIC

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

p_carb <- p_co2ex + p_DIC + p_SI + plot_layout(design = layout)

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
ptot
ggsave(plot = ptot,
       filename = file.path("results", "figure3_timeseries_v4.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 15.2)
