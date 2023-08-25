# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Gwen's measure data
df <- readxl::read_xlsx(file.path("data", "observed_pCO2_Gwen_compare.xlsx"),
                        sheet = "rawdata")

# Conductivity from EDF
df_cond <- readRDS(file.path("data", "05_hourly_data_clean", "cond_damup.RDS"))

# Get into plotting format
df <- df %>%
  mutate(hr = hour(datetime),
         month = month(datetime, label = TRUE, abbr = TRUE),
         group = paste(site,month, sep = "-")) %>%
  group_by(site, month) %>%
  mutate(dt = datetime,
         datetime = if_else(group == "Bonny-May", datetime - hours(6), datetime)) %>%
  left_join(df_cond) %>%
  mutate(day = ifelse(date(dt) == min(date(dt)), 1, 2),
         time = ymd_h(paste0("2000050", day, hr)))

# data for light
df_light <- df %>%
  mutate(light = streamMetabolizer::calc_light(dt, latitude = 47.7, longitude = 2.57)) %>%
  filter(light == 0) %>% 
  group_by(site, month, group) %>%
  dplyr::summarize(xmin = min(time),
                   xmax = max(time)) %>%
  ungroup()

# Measured data
p_meas <- ggplot(data = df,
       aes(x = time)) +
  geom_point(aes(y = Alka)) +
  geom_rect(data = df_light,
            aes(NULL, NULL, xmin= xmin, xmax = xmax, group = group),
            fill = "black",
            ymin = 0, ymax = 2.5, alpha=0.2) +
  facet_wrap(~group, scales = "free") +
  scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
  theme_classic(base_size = 10) +
  theme(panel.spacing.x = unit(2, "lines")) +
  labs(x = "hour of day",
       y = expression(A[T]~"("*mol~m^{-3}*")"))

# Conductivity data
p_cond <- ggplot(data = df,
                 aes(x = time,
                     y = SpC)) +
  geom_point(color = "red") +
  facet_wrap(~group, scales = "free") +
  scale_y_continuous(position = "right") +
  scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
  theme_classic(base_size = 10) +
  labs(x = "",
       y = expression(C[25]~"("*mu*S~cm^{-1}*")")) +
  theme(panel.spacing.x = unit(2, "lines"),
        rect = element_rect(fill = "transparent"),
        axis.text.y = element_text(color = "red"),
        axis.ticks.y = element_line(color = "red"),
        axis.title.y = element_text(color = "red"),
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        legend.position = "none",
        plot.tag = element_blank(),
        panel.border = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"))
# p_cond


layout <- c(area(t = 1, l = 1, b = 5, r = 5),
            area(t = 1, l = 1, b = 5, r = 5))
p <- p_meas + p_cond + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "supp_gwen_measurements_vs_cond.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12)

str(df)
mod <- lm(Alka ~ SpC+month, data = df)
summary(mod)
plot(mod)
ggplot(data = df,
       aes(x = SpC,
           y = Alka,
           # color = group,
           group = month)) +
  geom_point() +
  stat_summary(method = "lm") +
  facet_wrap(~month) +
  ggpubr::stat_regline_equation()
