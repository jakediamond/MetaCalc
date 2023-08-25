# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
df %>%
  filter(date == ymd(20200531)) %>%
  # between(hr, 5, 17)) %>%
  ggplot(aes(x = NEP, #- lag(HCO3),
             y = SpC - lag(SpC), #- lag(O2ex),
             color = hr)) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  ggpubr::stat_regline_equation()

df %>%
  filter(date == ymd(20200531)) %>%
  # between(hr, 5, 17)) %>%
  ggplot(aes(x = DIC_uM, #- lag(HCO3),
             y = exO2, #- lag(O2ex),
             color = hr)) + 
  geom_path(linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  ggpubr::stat_regline_equation()

# plot(x$CO2_uM, df$CO2_uM)


# First calculate slopes for exo2-exCO2
df_modco2 <- df %>%
  group_by(year, month, date) %>%
  drop_na(exDIC_uM, O2ex) %>%
  nest() %>%
  mutate(mod = map(data, ~lm(.$O2ex ~ .$exCO2_uM)),
         tid = map(mod, broom::tidy),
         gl = map(mod, broom::glance))

# Summary data
df_resco2 <- select(df_modco2, -data,- mod) %>%
  unnest(tid) 

# Get into tidy and long formats
df_dic <- select(df_mod, -mod, -data, -gl) %>%
  unnest(tid) %>%
  mutate(term = if_else(grepl("Int", term), "intercept", "slope")) %>%
  select(-(statistic:p.value)) %>%
  pivot_wider(names_from = term, values_from = c(estimate, std.error)) %>%
  left_join(select(df_mod, -mod, -data, -tid) %>%
              hoist(gl, r.squared = "r.squared", p.value = "p.value") %>%
              select(-gl)) %>%
  left_join(distinct(df, date, archetype)) %>%
  arrange(desc(r.squared), p.value, archetype) %>%
  select(-std.error_intercept)


seacarb::carb(flag = 11,
              var1 = 0.0015,
              var2 = 0.0017,
              T = 15,
              S = 0)


seacarb::carb(flag = 11,
              var1 = 0.0013,
              var2 = 0.0017,
              T = 15,
              S = 0)

seacarb::carb(flag = 11,
              var1 = 0.0013,
              var2 = 0.0015,
              T = 15,
              S = 0)


seacarb::carb(flag = 11,
              var1 = 0.0014,
              var2 = 0.0017,
              T = 15,
              S = 0)


x = seq(0.0015, 0.0012, -0.00001)
z = seacarb::carb(flag = 11, var1 = x, var2 = 0.0017, T = 15, S = 0)$DIC
plot(x, z)
lm(x~z)


o2 = seq(0.0012, 0.00135, 0.000005)
x2 = seq(0.0015, 0.0012, -0.00001)
z2 = seacarb::carb(flag = 11, var1 = x2, var2 = x2+0.0002, T = 15, S = 0)$DIC
plot(x2, z2)
lm(o2~z2)




x3 = seq(0.0015, 0.0012, -0.00001)
y3 = seq(0.0017, 0.0017 - 0.00001*2/3 * (length(x) - 1), -0.00001 * 2/3)
z3 = seacarb::carb(flag = 11, var1 = x3, var2 = y3, T = 15, S = 0)$DIC
plot(x3, z3)
lm(x3~z3)

h3 = seacarb::carb(flag = 11, var1 = x3, var2 = y3, T = 15, S = 0)$HCO3
lm(x3~h3)
lm(x3~y3)



x2 = seq(0.0015, 0.0012, -0.00001)
z2 = seacarb::carb(flag = 11, var1 = x2, var2 = x2+0.0002, T = 15, S = 0)$DIC
plot(x2, z2)
lm(x2~z2)


# Cco2 uptake
x4 = seq(0.000015, 0.000012, -0.0000001)
z4 = seacarb::carb(flag = 4, var1 = x4, var2 = 0.0017, T = 15, S = 0)$DIC
plot(x4, z4)
lm(-x4~z4)


h4 = seacarb::carb(flag = 4, var1 = x4, var2 = 0.0017, T = 15, S = 0)$HCO3
lm(-x4~h4) 
