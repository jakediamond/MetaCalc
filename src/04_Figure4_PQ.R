# -------------------------------------
# Author: Jake Diamond
# Purpose: to plot exO2 vs exCO2 and exDIC
# Date: 10 mars 2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)
# library(tidytable)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS"))
df_tf <- readRDS(file.path("data", "daily_trophlux.RDS"))

df <- df |>
  left_join(select(df_tf, date, trophlux)) |>
  mutate(exCO2 = CO2_uM - CO2eq_uM,
         exDIC = DIC_uM - DICeq_uM,
         flux = if_else(FCO2 > 0, "source", "sink"),
         troph = if_else(NEP > 0, "autotrophic", "heterotrophic"),
         trophlux_hr = paste(troph, flux)) |>
  drop_na(troph)
            
# How many samples were under vs oversaturated
sum(df$exO2 > 0, na.rm = T) / nrow(df)
sum(df$exCO2 > 0, na.rm = T) / nrow(df)

# Overall plot of exDIC vs ex o2
p_all_dic <- ggplot(data = df,
       aes(x = exDIC,
           y = exO2,
           color = trophlux_hr)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  scale_x_continuous(limits = c(-550, 760), breaks = seq(-500, 750, 250)) +
  theme(legend.position = c(0.155, 0.17),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  labs(x = expression("exDIC ("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))

p_all_dic

# Overall plot of exCO2 vs ex o2
p_all_co2 <- ggplot(data = df,
                    aes(x = exCO2,
                        y = exO2,
                        color = trophlux_hr)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(name = "trophlux state",
                     values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  labs(x = expression(exCO[2]~"("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")"))

p_all_co2

p_all <- p_all_dic + annotation_custom(ggplotGrob(p_all_co2),
                                                  ymin = 10,
                                                  ymax = 460,
                                                  xmin = 50,
                                                  xmax = 760)

# p_all

# Photosynthetic quotients ------------------------------------------------
# hourly exO2 vs del DIC each day
df_mod <- df |>
  group_by(year, month, date) |>
  drop_na(exDIC, exO2) |>
  nest() |>
  mutate(mod = map(data, ~lm(.$exO2 ~ .$DIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

df_res <- df_mod |>
  ungroup() |>
  hoist(gl, "r.squared") |>
  select(-data, -mod, -gl) |>
  unnest(tid) |>
  filter(term == ".$DIC_uM")

# Get just the good data for exO2 vs DIC relationship
df_pq <- df_res |>
  filter(p.value < 0.05, r.squared > 0.66) |>
  filter(between(estimate, -3, 1)) |>
  mutate(doy = yday(date)) |>
  mutate(regime = if_else(year(date) < 2012, "planktonic", "benthic"))

# Check for bimodality
multimode::modetest(df_pq$estimate)

# Plot of slope (=PQ) over year based on ecosystem regime
p_pq_t <- ggplot(data = df_pq,
                 aes(x = doy,
                     y = estimate,
                     color = regime)) +
  stat_smooth() + 
  theme_classic(base_size = 10) +
  ggsci::scale_color_aaas(name = "autotrophic regime") +
  scale_x_continuous(breaks = seq(0, 365, 90)) +
  scale_y_continuous(breaks = seq(-1.4, 0, 0.2)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  # annotate(geom = "text", x= 200, y = -2.06, label = "calcite precipitation",
  #          size = 3.5) +
  # annotate(geom = "text", x= 200, y = -0.97, label = expression(atop(CO[2]~"/"~HCO[3]^{`-`}, uptake)),  
  #          size = 3.5) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  # geom_hline(yintercept = -2) +
  labs(x = "day of year",
       y = expression(PQ~"("*beta[DIC]*")"))
p_pq_t

# dataframe for medians
df_text <- df_pq |>
  filter(between(doy, 90, 270)) |>
  group_by(regime) |>
  summarize(pq = median(estimate, na.rm = T)) |>
  mutate(x = c(0, -2), y = 1.4,
         label = paste("median PQ =", signif(pq,2)))

# Density plots
p_pq_dens <- ggplot(data = filter(df_pq, between(doy, 90, 270)),
                    aes(x = estimate,
                        group = regime)) +
  geom_density(aes(fill = regime), alpha = 0.5) +
  stat_summary(aes(xintercept = after_stat(x), y = 0, color = regime), fun = median,
               linetype = "dashed", geom = "vline", linewidth = 1, orientation = "y") +
  geom_text(data = df_text, aes(x = x, y = y, label = label, color = regime)) +
  theme_classic(base_size = 10) +
  guides(text = "none", color = "none") +
  ggsci::scale_fill_aaas(name = "autotrophic regime") +
  ggsci::scale_color_aaas(name = "autotrophic regime") +
  theme(legend.position = c(0.17, 0.5),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = expression(daily~exO[2]~-~DIC~slope~"("*beta[DIC]*")"),
       y = "density")
p_pq_dens
p_all_slope <- p_pq_dens + annotation_custom(ggplotGrob(p_pq_t),
                                             ymin = 0.1,
                                             ymax = 1.2,
                                             xmin = -0.7,
                                             xmax = 1)

# p_all_slope

p <- p_all / p_all_slope + plot_annotation(tag_levels = "a")
# p
ggsave(plot = p,
       filename = file.path("results", "figure4_PQ.png"),
       dpi = 300,
       units = "cm",
       height = 15,
       width = 13.5)



