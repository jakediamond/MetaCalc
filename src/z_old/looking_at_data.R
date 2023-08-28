# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

x  = filter(df, date == ymd_h(2020072600))
yday(x$date)
colnames(df_met)
tail(y)
df_rat <- df_met %>%
  filter(grepl("dampierre_up", source)) %>%
  group_by(date) %>%
  summarize(pr = GPP_mean / abs(ER_mean))

# hourly AOU vs del DIC each day
df_mod <- df %>%
  group_by(year, month, date) %>%
  drop_na(exDIC_uM, O2ex) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$O2ex ~ .$exDIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

df_res <- select(df_mod, -data,- mod) %>%
  unnest(tid)

df_res %>%
  filter(term == ".$exDIC_uM") %>%
  mutate(doy = yday(date)) %>%
  ggplot(aes(x = doy, y = estimate, group = doy)) +
  stat_summary() +
  geom_hline(yintercept = 0.5)


z <- df_res %>%
  mutate(date = as.Date(date)) %>%
  left_join(df_rat)

z %>%
  filter(term == ".$exDIC_uM",
         pr > 0) %>%
  ggplot(aes(x = pr, y = estimate)) +
  stat_summary_bin() +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(-3, 0)) +
  geom_abline(slope = -1)



df_mod2 <- df %>%
  group_by(year, month, date) %>%
  mutate(ALK_uM = 1e6*Alk_molkg) %>%
  drop_na(exDIC_uM, Alk_molkg) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$ALK_uM ~ .$exDIC_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))


df_res2 <- select(df_mod2, -data,- mod) %>%
  unnest(tid)


z2 <- df_res2 %>%
  mutate(date = as.Date(date)) %>%
  left_join(df_rat)

df_res2 %>%
  filter(term == ".$exDIC_uM") %>%
  mutate(doy = yday(date)) %>%
  ggplot(aes(x = doy, y = estimate, group = doy)) +
  stat_summary() +
  geom_hline(yintercept = 0.5)
  # scale_x_continuous(limits = c(0, 2)) +
  # scale_y_continuous(limits = c(-3, 0)) +
  # geom_abline(slope = -1)




df_mod3 <- df %>%
  group_by(year, month, date) %>%
  mutate(ALK_uM = 1e6*Alk_molkg) %>%
  drop_na(O2ex, Alk_molkg) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$O2ex ~ .$ALK_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))


df_res3 <- select(df_mod3, -data,- mod) %>%
  unnest(tid)


z3 <- df_res3 %>%
  mutate(date = as.Date(date)) %>%
  left_join(df_rat)

df_res3 %>%
  filter(term == ".$ALK_uM") %>%
  mutate(doy = yday(date)) %>%
  ggplot(aes(x = doy, y = estimate, group = doy)) +
  stat_summary() +
  geom_hline(yintercept = -1)


z3 %>%
  filter(term == ".$ALK_uM",
         pr > 0) %>%
  ggplot(aes(x = pr, y = estimate)) +
  stat_summary_bin() +
  scale_x_continuous(limits = c(0, 2)) +
  # scale_y_continuous(limits = c(-3, 0)) +
  geom_abline(slope = -1)



df_met %>%
  filter(date == ymd(19990824))


df %>%
  filter(date == ymd_h(1999082400)) %>%
  ggplot(aes(y = SpC - lag(SpC),
             x = NEP_mmolO2m3,
             color = hr)) +
  scale_color_viridis_c() +
  geom_path(size = 2) +
  ggpubr::stat_regline_equation()

df %>%
  filter(date == ymd_h(1999082400)) %>%
  ggplot(aes(x = exDIC_uM - lag(exDIC_uM),
             y = Alk_molkg*1E6 - lag(Alk_molkg*1E6))) +
  geom_point() +
  ggpubr::stat_regline_equation()


df %>%
  filter(date == ymd_h(1999082400)) %>%
  ggplot(aes(x = Alk_molkg*1E6 - lag(Alk_molkg*1E6),
             y = NEP_mmolO2m3,
             color = hr)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()


df_met %>%
  filter(date == ymd(20030825))

df %>%
  filter(date == ymd_h(2003082500)) %>%
  ggplot(aes(x = exDIC_uM - lag(exDIC_uM),
             y = NEP_mmolO2m3)) +
  geom_point() +
  ggpubr::stat_regline_equation()


df %>%
  filter(date == ymd_h(2003082500)) %>%
  ggplot(aes(x = Alk_molkg*1E6 - lag(Alk_molkg*1E6),
             y = NEP_mmolO2m3,
             color = hr)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()


df %>%
  filter(date == ymd_h(2003082500)) %>%
  ggplot(aes(x = Alk_molkg*1E6,
             y = O2ex,
             color = NEP_mmolO2m3)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()




a = df %>%
  filter(date == ymd_h(2005070900))

df_met %>%
  filter(date == ymd(20050709))

df %>%
  filter(date == ymd_h(2005070900)) %>%
  ggplot(aes(y = alkdif_uM,
             x = NEP_mmolO2m3,
             color = hr)) +
  scale_color_viridis_c() +
  geom_path(size = 2) +
  ggpubr::stat_regline_equation()




df %>%
  filter(date == ymd_h(2005070900)) %>%
  ggplot(aes(x = exDIC_uM,
             y = O2ex,
             color = hr)) +
  scale_color_viridis_c() +
  geom_path(size = 2) +
  ggpubr::stat_regline_equation()

df %>%
  filter(date == ymd_h(2005070900)) %>%
  ggplot(aes(x = exDIC_uM,
             y = O2ex,
             color = NEP_mmolO2m3)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()


df %>%
  filter(date == ymd_h(2005070900)) %>%
  ggplot(aes(x = Alk_molkg*1E6 - lag(Alk_molkg*1E6),
             y = NEP_mmolO2m3,
             color = hr)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()


df %>%
  filter(date == ymd_h(2005070900)) %>%
  ggplot(aes(x = Alk_molkg*1E6,
             y = O2ex,
             color = NEP_mmolO2m3)) +
  geom_path(size = 1.5) +
  ggpubr::stat_regline_equation() +
  scale_color_viridis_c()
