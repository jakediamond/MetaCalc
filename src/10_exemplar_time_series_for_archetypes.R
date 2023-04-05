# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(tidyverse)
library(tidytable)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete.RDS"))
df <- mutate(df,
             delDIC = DIC - lag(DIC),
             NEP_mmolCm2hr = -(delDIC + calc_mmol_m2_hr - CO2_flux_enh_mmol_m2_hr),
             totC = calc_mmol_m2_hr + NEP_mmolCm2hr) 

# Plotting ----------------------------------------------------------------
# Get into plot types
fluxes <- c("CO2_flux_enh_mmol_m2_hr", "NEP_mmolO2m3", "calc_mmol_m2_hr", "delDIC")
meas <- c("O2ex", "pH", "SpC")
carbs <- c("CO2", "CO3", "HCO3", "DIC")
calc <- c("LSI", "CCPP", "calc")
oxyco <- c("O2ex", "exCO2_uM", "exDIC_uM")

# Get hourly averages of everything 
df_avg <- df %>%
  pivot_longer(cols = -c(year, month, archetype, date, datetime, solartime,
                         hr, depth,discharge)) %>%
  group_by(archetype, name, hr) %>%
  summarize(mean_val = mean(value, na.rm = T),
            se_val = sd(value, na.rm = T) / sqrt(n()) * 1.96) %>%
  ungroup() %>%
  mutate(trophic = str_split_i(archetype, "_", 1),
         flux = str_split_i(archetype, "_", 2),
         archetype = str_replace(archetype, "_", " "))

# Get a plot of hourly mean fluxes for different archetypes
df_p_fluxes <- filter(df_avg, name %in% fluxes) %>%
  drop_na(archetype)

p_fluxes <- ggplot(data = df_p_fluxes,
       aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = name)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = name), alpha = 0.4) +
  facet_wrap(~archetype, ncol = 4) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  # scale_color_manual(values = c("black", "#E69F00", "#0072B2")) +#, "#009E73")) +
  # scale_fill_manual(values = c("black", "#E69F00", "#0072B2")) +#, "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.85, 0.15), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  labs(x = "solar hour",
       y = expression(flux~"("*mmol~C~m^{-2}~hr^{-1}*")"))
p_fluxes

# Get a plot of hourly mean CO2 for different archetypes
df_p_meas <- filter(df_avg, name %in% meas) %>%
  drop_na(archetype)

ggplot(data = df_p_meas,
       aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
                # color = name)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
                  # fill = name), 
              alpha = 0.4) +
  facet_wrap(~name, scales = "free") +
  # facet_grid(flux~trophic) +
  # geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  # scale_color_manual(values = c("#648FFF", "#DC267F", "#FFB000")) +
  # scale_fill_manual(values = c("#648FFF", "#DC267F", "#FFB000")) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.85), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  labs(x = "solar hour",
       y = expression(excess~O[2]*","~CO[2]*","~or~DIC~"("*mu*M*")"))
       # y = expression(flux~"("*mmol~C~m^{-2}~hr^{-1}*")"))

quantile(df$pH, na.rm = T)

ca_rat <- filter(df,
                 archetype == "autotrophic_sink") %>%
  mutate(rat = calc_mmol_m2_hr / NEP_mmolO2m3,
         calc2 =if_else(LSI > 0.4,
                       -0.5 * alkdif_uM * depth,
                       0)) %>%
  ungroup() %>%
  group_by(date) %>%
  filter(light >0) %>%
  summarize(calc = sum(calc2),
            nep = sum(NEP_mmolO2m3),
            nepC = sum(NEP_mmolCm2hr),
            nep_d = mean(NEP_daily),
            rat = mean(rat)) %>%
  mutate(rat2 = calc / nep_d,
         ratC = calc / nepC)

filter(ca_rat, calc > 0) %>%
  ungroup() %>%
  mutate(rat3 = calc / nep) %>%
  summarize(r = median(rat),
            r2 = -median(rat2),
            r3 = median(rat3),
            rC = median(ratC))

(ca_rat %>%
  ungroup() %>%
  # mutate(doy = )
  ggplot(aes(x = date))+
  geom_line(aes(y = nepC), color = "blue") +
    geom_line(aes(y = calc), color = "red") #+ 
  # scale_x_log10() +
  # scale_y_log10()
) %>%
  ggplotly()

mean(df$discharge, na.rm = T)
# Get a plot of hourly mean CO2 for different archetypes
df_p_fluxes <- select(df, hr, archetype, any_of(fluxes)) %>%
  pivot_longer(cols = contains(fluxes)) %>%
  drop_na(archetype)

ggplot(data = df_p_fluxes,
       aes(x = hr,
           y = value,
           color = name)) +
  stat_smooth() +
  facet_wrap(~archetype) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2"))+#, "#009E73")) +
  theme_bw()+
  labs(x = "solar hour",
       y = expression(flux~"("*mmol~C~m^{-2}~hr^{-1}*")"))

# Get a plot of hourly mean measurements for different archetypes
df_p_meas <- select(df, hr, archetype, any_of(calc)) %>%
  pivot_longer(cols = contains(calc)) %>%
  # mutate(value = if_else(name == "CO2_flux_enh", value / 24, value),
  #        value = if_else(name == "NEP", -1 * value, value)) %>%
  drop_na(archetype)

ggplot(data = df_p_meas,
       aes(x = hr,
           y = value,
           color = archetype)) +
  stat_smooth() +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "solar hour",
       y = expression(measurement~"("*mmol~m^{-3}*")"))


df %>%
  group_by(date, archetype) %>%
  summarize(NEP_daily = mean(NEP_daily),
            CO2_daily = mean(CO2_daily),
            CO2_hrs = sum(CO2_flux_enh),
            NEP_hrs = sum(-NEP),
            calc = sum(calc)) %>%
  ungroup() %>%
  group_by(archetype) %>%
  summarize(across(everything(), mean, na.rm = T))



p <- plot_ly(data  = df) %>%
  add_trace(
    x = ~datetime, y = ~ alkalinity, type = "scatter", mode='lines', 
    color = I("#1E88E5")) %>%
  add_trace(data = df_calc,
            x = ~datetime, y = ~ Alk_molkg*1000, type = "scatter", mode='lines', 
            color = I("#FFC107")) %>%
  layout(
    xaxis = list(title = ""),
    yaxis = list(title = TeX("\\color{#1E88E5}{LSI~(-)}")),
    title = TeX("\\text{CO_{2}}")) %>%
  config(mathjax = "cdn")
htmltools::browsable(p)



df_NEP <- df %>%
  mutate(delDIC = DIC - lag(DIC)) %>%
  mutate(NEP_C = delDIC + calc + CO2_flux_enh,
         NEP_O = NEP) %>%
  filter(light > 0) %>%
  group_by(date, archetype, year, month) %>%
  summarize(NEP_c = sum(NEP_C, na.rm = T),
            NEP_o = sum(NEP_O, na.rm = T),
            DIC_rem = sum(delDIC, na.rm = T))

df_NEP %>%
  filter(year > 1992, month == 8) %>%
  mutate(ratio = NEP_c/ NEP_o,
         ratio2 = DIC_rem / NEP_c) %>%
  ungroup() %>%
  summarize(m = mean(ratio, na.rm = T),
            m2 = mean(ratio2, na.rm = T))



x = df_co2 %>%
  mutate(date = date(datetime)) %>%
  group_by(date) %>%
  summarize(co2_flux = sum(co2_flux_mmol_m2_d/24)) %>%
  left_join(select(df_nepco2, date, value_CO2_mean))

y = df_o2co2 %>%
  mutate(date = date(datetime)) %>%
  group_by(date) %>%
  summarize(nep = mean(NEP*24)) %>%
  left_join(select(df_nepco2, date, filtered_NEP_mean))

l = left_join(df_NEP, df_nepco2)
plot(l$NEP_o, -l$value_NEP_mean)

df$exDIC = ifelse(is.na(df$DIC), NA_real_, exDIC*1000)
p_co2o2 <- plot_ly(data = df,
                   x = ~datetime) %>%
  add_trace(
    y = ~ NEP_mmolCm2hr, type = "scatter", mode='lines', 
    color = I("#1E88E5")) %>%
  add_trace(y = ~calc_mmol_m2_hr, type = "scatter", mode='lines', 
            color = I("#FFC107")) %>%
  # add_trace(y = ~exDIC_uM, type = "scatter", mode='lines', 
  #           color = I("red")) %>%
  layout(
    xaxis = list(title = ""),
    yaxis = list(title = TeX("\\color{#1E88E5}{LSI~(-)}")),
    title = TeX("\\text{CO_{2}}")) %>%
  config(mathjax = "cdn")
htmltools::browsable(p_co2o2)






# Daily comparisons -------------------------------------------------------
df_d <- df %>%
  select(-datetime, -solartime, -light) %>%
  group_by(year, month, date, archetype) %>%
  summarize(NEP_O2 = mean(NEP_mmolO2m3) *24,
            NEP_met = mean(NEP_daily),
            NEP_C = mean(NEP_mmolCm2hr) * 24,
            exDIC = mean(exDIC_uM),
            exDO = mean(O2ex),
            exCO2 = mean(exCO2_uM))


# positive relationship HCO3 -> CO2
seacarb::carb(flag = 1,
              8,
              6E-5,
              S = 0)



# Figures -----------------------------------------------------------------
p_o2ex <- ggplot(data = filter(df_avg, name == "O2ex") %>%
                   drop_na(archetype),
                 aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0, 23, 4)) +
  scale_y_continuous(limits = c(-110, 200)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  guides(fill = "none") + 
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        # legend.direction = "horizontal",
        # legend.background = element_rect(fill = "transparent"),
        # legend.key.size = unit(0.25, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black"),
        # plot.tag = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "solar hour",
       y = expression(exO[2]~"("*mu*M*")"))
p_o2ex

# Specific conductivity
p_C <- ggplot(data = filter(df_avg, name == "SpC") %>%
                 drop_na(archetype),
               aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(250, 290), position = "right") +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        # plot.tag = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_C

# pH
p_pH <- ggplot(data = filter(df_avg, name == "pH") %>%
                 drop_na(archetype),
               aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(7.5, 9.5)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  labs(y = "pH",
       x = "solar hour")
p_pH

# Layout
layout <- c(area(t = 1, l = 1, b = 5, r = 5),
            area(t = 5.5, l = 1, b = 9.5, r = 5),
            area(t = 9.9, l = 1, b = 14, r = 5))

p <- p_o2ex + p_C + p_pH + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "measurement_archetypes.png"),
       dpi = 300,
       units = "cm",
       width = 9.2,
       height = 15)

# Plot carbonate system ---------------------------------------------------
p_co2ex <- ggplot(data = filter(df_avg, name == "exCO2_uM") %>%
                   drop_na(archetype),
                 aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(-25, 80)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black"),
        # plot.tag = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "solar hour",
       y = expression(exCO[2]~"("*mu*M*")"))
p_co2ex


p_HCO3 <- ggplot(data = filter(df_avg, name == "HCO3") %>%
                 drop_na(archetype),
               aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous(limits = c(1500, 1900), position = "right") +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        # plot.tag = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = expression(HCO[3]^{`-`}~"("*mu*M*")"))
p_HCO3


p_LSI <- ggplot(data = filter(df_avg, name == "LSI") %>%
                 drop_na(archetype),
               aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  labs(y = "LSI (-)",
       x = "solar hour")
p_LSI

# Plot Layout
layout <- c(area(t = 1, l = 1, b = 3, r = 5),
            area(t = 4, l = 1, b = 6, r = 5),
            area(t = 7, l = 1, b = 9, r = 5))

p2 <- p_co2ex + p_HCO3 + p_LSI + plot_layout(design = layout)
p2


p_all <- p | p2
p_all
ggsave(plot = p_all,
       filename = file.path("results", "timeseries_archetypes.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 15)



# Fluxes ------------------------------------------------------------------
p_nep <- ggplot(data = filter(df_avg, name == "NEP_mmolO2m3") %>%
                    drop_na(archetype),
                  aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  # scale_y_continuous(limits = c(-25, 80)) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black"),
        axis.line.x = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "solar hour",
       y = expression(NEP~"("*mmol~m^{-2}~h^{-1}*")"))
p_nep


p_fco2 <- ggplot(data = filter(df_avg, name == "CO2_flux_enh_mmol_m2_hr") %>%
                   drop_na(archetype),
                 aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  scale_y_continuous( position = "right") +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = expression(fCO[2]~"("*mmol~m^{-2}~h^{-1}*")"))
p_fco2

p_calc <- ggplot(data = filter(df_avg, name == "calc_mmol_m2_hr") %>%
                  drop_na(archetype),
                aes(x = hr)) +
  geom_line(aes(y = mean_val,
                color = archetype)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                  fill = archetype),
              alpha = 0.4) +
  scale_x_continuous(breaks = seq(0,23, 4)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major.x = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  labs(y = expression(CaCO[3]~ppt.~"("*mmol~m^{-2}~h^{-1}*")"),
       x = "solar hour")
p_calc

# Plot with kuramoto
layout <- c(area(t = 1, l = 1, b = 3, r = 5),
            area(t = 4, l = 1, b = 6, r = 5),
            area(t = 7, l = 1, b = 9, r = 5))

p3 <- p_nep + p_fco2 + p_calc + plot_layout(design = layout)
p3

ptot <- (p | p2 | p3) + plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") & 
  guides(color = "none", fill = guide_legend("trophlux state",
                                             override.aes = list(alpha=1))) &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.key.size = unit(0.25, "cm"),
        plot.tag.position = c(0.25, 0.95))
ptot
ggsave(plot = ptot,
       filename = file.path("results", "timeseries_archetypes_overall.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 15.2)


df %>%
  group_by(date) %>%
  summarize(CP = mean(CCPP),
            SI = mean(LSI)) %>%
ggplot(aes(x = SI,
           y = CP)) +
  geom_point() +
  # scale_x_log10() +
  scale_y_log10()

summary(lm(log(CP)~SI, data =ll))
exp(-3.18)*exp(2.53*-1)

(ggplot(data = df,
       aes(x = datetime,
           y = Alk_molkg)) +
  geom_line()) %>%
  ggplotly()
