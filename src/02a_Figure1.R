# -------------------------------------
# Author: Jacob Diamond
# Purpose: to plot Figure 1 for the manuscript
# Date: 04-04-2023
# -------------------------------------
# Load libraries
library(patchwork)
library(tidyverse)

# Load data ---------------------------------------------------------------
df <- readRDS(file.path("data", "hourly_data.RDS"))

# Get into daily and get trophlux states
df_d <- df |>
  select(date, datetime, FCO2_enh) |>
  mutate(FCO2_hr = FCO2_enh /24) |>
  group_by(date) |>
  summarize(FCO2 = sum(FCO2_hr, na.rm = T)) |>
  left_join(distinct(df, date, Q_m3s,
                     NEP = NEP_mean)) |>
  mutate(troph = if_else(NEP > 0, "autotrophic", "heterotrophic"),
         flux = if_else(FCO2 > 0, "source", "sink"),
         trophlux = paste(troph, flux))

# Fraction of year by trophlux state -------------------------------------------
p_doy <- df_nepco2 %>%
  mutate(wyear = if_else(month > 8, year + 1, year),
         origin = ymd(paste0(wyear-1, "1001"))) %>%
  mutate(dowy = as.numeric(interval(origin, date))/86400 + 1,
         datewy = ymd(paste(if_else(month > 8, "1970", "1971"), 
                            month, day(date), sep = "-"))) %>%
  ggplot(aes(x = datewy,
             y = (..count..),
             fill = archetype)) +
  geom_bar(position = "fill", width = 1) +
  theme_bw(base_size = 10) +
  scale_fill_manual(values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_date(date_labels = "%b", date_breaks = "61 days",
               # limits 
               expand = expansion(mult = c(0, 0))) +
  # scale_x_continuous(breaks = seq(0, 365, 30), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(y = "32-yr occurence probability (-)",
       x = "") + 
  theme(legend.position = "none", axis.title.x = element_blank())
p_doy

# flux by discharge bin ---------------------------------------------------
# mean annual cumulative flux by discharge bin and archetype
df_cum <- df_nepco2 %>%
  mutate(wyear = ifelse(month > 9, year + 1, year),
         # Create 30 bins in log10 space
         qbins = cut(log10(Q), 30)) %>%
  group_by(wyear) %>%
  # Add a column for total annual FCO2 = the daily sum
  mutate(totflux = sum(filtered_CO2_meanenh, na.rm= T)) %>%
  ungroup() %>%
  group_by(qbins, archetype, wyear) %>%
  # Now get total annual flux by discharge bin and trophlux
  summarize(flux = sum(filtered_CO2_meanenh, na.rm = T),
            # Total annual flux stays the same
            totflux = mean(totflux))

# Proportion of annual fco2 occuring in each trophlux and discharge bin
df_fco2_prop <- ungroup(df_cum) %>%
  # Proportion by water year, discharge bin, and trophlux
  mutate(per = flux / totflux) %>%
  group_by(qbins, archetype) %>%
  # now mean across all water years
  summarize(per = mean(per)) %>%
  mutate(qbinend = 10^as.numeric(str_extract(qbins, "\\b[0-9.]+")),
         archetype = str_replace(archetype, "_", " "))

p_co2_q_hist <- ggplot(data = df_fco2_prop,
       aes(x = qbinend,
           y = per * 100,
           fill = archetype,
           group = archetype)) +
  geom_col() +
  # geom_bar() +
  # geom_line() +
  # geom_point() +
  theme_classic(base_size = 10) +
  geom_vline(xintercept = 300, linetype = "dashed") +
  # facet_wrap(~archetype) +
  # scale_y_log10() +
  scale_fill_manual(name = "trophlux state",
                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  geom_hline(yintercept = 0) +
  theme(legend.position = c(0.24, 0.84),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_rect(fill = "transparent")) +
  labs(y = expression("contribution to annual "*CO[2]~ "flux (%)"),
       x = expression("mean daily discharge"~"("*m^3~s^{-1}*')'))
p_co2_q_hist
p_final <- p_doy / p_co2_q_hist + plot_annotation(tag_levels = "a")

ggsave(plot = p_final,
       filename = file.path("results", "figure1.png"),
       dpi = 300,
       units = "cm",
       width = 8.9,
       height = 14)
