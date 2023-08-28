x <- df %>%
  select(date, datetime, hr, temp, Ca, HCO3, CO3, CO2_uM, Alk_molkg, LSI, pH, 
         SpC, exDIC_uM, exCO2_uM, O2ex, NEP_mmolO2m3)

x <- x %>%
  mutate(act = 10^-pfm_fun(temp + 273.15, SpC),
         Ksp = Ksp_fun(temp + 273.15),
         omega = ((Ca * act^4/ 40.078/1E3) * CO3 * act^4/1E6) / Ksp,
         SI = log10(omega))

ggplot(data= x,
       aes(x = LSI,
           y = omega)) +
  # geom_point()
  stat_summary_bin() +
  geom_abline() +
  scale_y_log10()


dplyr::filter(x, date == ymd(20210301)) %>%
  pivot_longer(cols = c(exDIC_uM, exCO2_uM)) %>%
  ggplot(aes(x = value,
             y = O2ex,
             shape = name,
             color = NEP_mmolO2m3)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  # scale_color_manual(name = "trophlux state",
  #                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  # geom_abline(intercept = 0, slope = -2, linetype = "dotted") +
  # geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.15, 0.15),
        legend.key.size = unit(0.3, "cm")) +
  # guides(colour = guide_legend(override.aes = list(alpha=1))) +
  labs(x = expression("exDIC or "*CO[2]*"("*mu*M*")"),
       y = expression(exO[2]~"("*mu*M*")")) +
  ggpubr::stat_regline_equation(label.x.npc = "center") +
  stat_smooth(method = 'lm')

dplyr::filter(x, date == ymd(20050201)) %>%
  mutate(Alk_uM = Alk_molkg * 1E6) %>%
  # pivot_longer(cols = c(SpC, Alk_uM)) %>%
  ggplot(aes(x = datetime,
             y = SpC,
             color = pH)) +
  geom_path(size = 2) +
  scale_color_viridis_c() +
  # scale_color_manual(name = "trophlux state",
  #                    values = c("black", "#E69F00", "#0072B2", "#009E73")) +
  # geom_abline(intercept = 0, slope = -2, linetype = "dotted") +
  # geom_abline(intercept = 0, slope = -1, linetype = "dashed") +
  # geom_hline(yintercept = 0) +
  # geom_vline(xintercept = 0) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.15, 0.15),
        legend.key.size = unit(0.3, "cm")) +
  # guides(colour = guide_legend(override.aes = list(alpha=1))) +
  labs(y = expression("Specific conductance ("*mu*S~cm^{-1}*")"),
       x = "")


b = dplyr::filter(x, date == ymd(20050701))
mean(b$pH)


colnames(df)
SPC_fun <- function(temp, cond) {
  cond / (1 - ((20 - temp) * 0.021))
}

b$sc <- SPC_fun(b$temp, 150)
b$fm <- 10^-pfm_fun(b$temp, 1.6*10^-5 * 300)
plot(b$hr, b$temp)
plot(b$temp, b$fm)


ggplot(data = df,
       aes(x = enh)) +
  stat_density(aes(y = after_stat(scaled))) +
  theme_classic() +
  facet_grid(trophic~flux) + 
  scale_x_continuous(limits = c(0, 10)) +
  labs(x = expression(alpha[enh]~"("*`-`*")"),
       y = "scaled density")


c <- dplyr::filter(df, date %in% c(ymd(20040920), ymd(20040921)))

df_d <- readRDS(file.path("data", "03_CO2", "dampierre_all_daily_data.RDS"))

df_p <- select(df_d, date, alk_molm3) %>%
  mutate(hr = 9,
         datetime = ymd_h(paste(date,hr))) %>%
  right_join(df)
(ggplot(data= df_p,
        aes(x = datetime)) +
    geom_line(aes(y = Alk_molkg * 1000)) +
    geom_point(aes(y = alk_molm3), color = "red", size = 2)) %>%
  ggplotly()



(ggplot(data= dplyr::filter(df, year == 2005),
        aes(x = datetime)) +
    geom_line(aes(y = Alk_rf)) +
    geom_line(aes(y = Alk_molkg * 1000), color = "blue") +
    geom_point(aes(y = Alk_molm3), color = "red", size = 2)) %>%
  ggplotly()

df2 <-drop_na(df)

1 - sum(((df2$Alk_molkg*1000) - df2$Alk_molm3)^2) / 
  sum((df2$Alk_molm3 - mean(df2$Alk_molm3))^2)


x = dplyr::filter(df_test, between(date, ymd(20070928), ymd(20071030))) %>%
  select(date, SpC, Ca, Alk_molkg) %>%
  group_by(date) %>%
  summarize(across(where(is.numeric), mean))
