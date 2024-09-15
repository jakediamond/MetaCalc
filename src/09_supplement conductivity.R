df <- readRDS(file.path("data", "hourly_data.RDS"))


df_p <- df %>%
  group_by(date) %>%
  filter((max(SpC) - min(SpC)) > 10) %>%
  filter(month(date) == 7,
         Q_m3s < 150,
         year == 2020)
  
ggplot(df_p,
       aes(x = datetime,
           y = SpC,
           color = NEP)) +
  geom_path() +
  scale_color_viridis_c() +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%d") +
  theme_classic()+
  theme(panel.grid.major.x = element_line()) +
  coord_cartesian(expand = FALSE) +
  labs(x = "day in July 2020",
       y = expression(C[25]~"("*mu*S~cm^{-1}*")"))

ggsave(filename = file.path("results", "supplement_cond_july2020.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12.2)
