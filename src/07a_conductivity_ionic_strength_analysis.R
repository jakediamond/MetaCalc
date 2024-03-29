# 
# Purpose: To estimate mean molar concentrations of important WQ constituents
# Author: Jake Diamond
# Date: 28 February 2023
# 

# Load libraries
library(plotly)
library(broom)
library(patchwork)
library(tidyverse)

# Load functions ----------------------------------------------------------
# pfm function, fm is activity coefficient for monovalent species
# I is ionic strength 
# temp is temperature in Celsius
fm_fun <- function(temp, I){
  TK = temp + 273.15 #temp in Kelvin
  # I = 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  E = 60954/(TK+116) - 68.937 #dielectric constant
  # Davies estimate
  pfm = 1.82 * 10^6 * (E*TK)^-1.5 * ((sqrt(I)/(1+sqrt(I))) - 0.3 * I)
  fm = 10^-pfm
}

# Load data---------------------------------------------------------------
# Load water chemistry data, most in mg/L
df_wq <- readxl::read_xlsx(file.path("data", "00_water chemistry", 
                                     "WQ_Loire_Moyenne.xlsx"))

# Quick clean
df_wq <- dplyr::mutate(df_wq,
                dplyr::across(everything(), ~ifelse(. == -1, NA_real_, .)))

# molecular weights (g/mol) and chem info, use just Si for SiO2
df_chem <- tibble(
  solute = c("NH4", "NO2", "NO3", "SiO2", "PO4", "SO4", "HCO3", "Na", "K", "Cl", "Ca", "Mg"),
  MW = c(18.04, 46.005, 62.0049, 28.0855, 94.9714, 96.06, 61.0168, 22.9898, 39.0983, 35.453, 40.078, 24.305),
  charge = c(1, 1, 1, 0, 3, 2, 1, 1, 1, 1, 2, 2),
  # Equivalent charge at infinite dilution at 25C (S cm^2 / mol) from Miller et al 
  # (USGS) 1988 via  Harned and Owen (1964, p. 231), also in CRC Handbook 2000
  eq_cond = c(73.5, 71.8, 71.42, NA_real_, 92.8, 80.0, 44.5, 50.08, 73.48, 76.31, 59.47, 53.0)
)

# Seasonal summary of concentrations --------------------------------------
# Get seasonal summaries in uM
df_wq_sum <- df_wq %>%
  select(-ChOD) %>% #no data here
  mutate(season = case_when(month %in% c(12, 1, 2) ~ "winter",
                            month %in% c(3, 4, 5) ~ "spring",
                            month %in% c(6, 7, 8) ~ "summer",
                            month %in% c(9, 10, 11) ~ "autumn")) %>%
  pivot_longer(cols = pH:Mg, names_to = "solute", values_to = "conc") %>%
  group_by(season, solute) %>%
  summarize(mean_conc = mean(conc, na.rm = T),
            sd_conc = sd(conc, na.rm = T)) %>%
  left_join(df_chem) %>%
  mutate(mean_uM = mean_conc * 1000 / MW,
         sd_uM = sd_conc * 1000 / MW)

# Table format
tab <- df_wq_sum %>%
  drop_na() %>%
  select(season, solute, mean_uM, sd_uM) %>%
  dplyr::mutate(season = as.factor(season),
         solute = as.factor(solute)) %>%
  mutate(season = fct_relevel(season, "winter", "spring", "summer", "autumn"),
         solute = fct_relevel(solute, "Ca", "Mg", "NH4", "K", "Na", 
                             "Cl", "NO2", "NO3", "HCO3", "SO4", "PO4"),
         text = paste0(round(mean_uM), "±", round(sd_uM))) %>%
  select(-mean_uM, -sd_uM) %>%
  arrange(solute) %>%
  pivot_wider(names_from = solute, values_from = text) %>%
  arrange(season)
# write_csv(tab, file.path("results", "wq_summary.csv"))

# Estimate specific conductivity and look at calcium contribution ---------
# Estimate conductivity with ions
df_cond <- df_wq %>%
  select(date, site_no, year, month, day, temp, pH, SpC_m = SC, NH4,
         NO2, NO3, PO4, SO4, HCO3, Na, K, Cl, Ca, Mg) %>%
  # Only want days when all ions are measured
  filter(complete.cases(select(., -c(NH4)))) %>%
  pivot_longer(NH4:Mg, names_to = "solute", values_to = "conc") %>%
  left_join(df_chem) %>%
  mutate(
    # Concentration in mol/L
    conc_M = conc / MW / 1000) %>%
  group_by(site_no, year, month, day, date, temp) %>%
  mutate(
    # Ionic strength
    I = 0.5 * sum(conc_M * charge^2),
    # monovalent activity
    fm = fm_fun(temp, I),
    # divalent activity
    fd = fm^4,
    # trivalent activity,
    ft = fm^9,
    # activity, this case_when is super slow
    act = case_when(
      charge == 1 ~ fm * conc_M,
      charge == 2 ~ fd * conc_M,
      charge == 3 ~ ft * conc_M,
      TRUE ~ NA_real_),
    # Theoretical conductivity based on concentration (uS/cm)
    SpC_t = sum(conc_M * charge * eq_cond) * 1000,
    # Theoretical conductivity based on activity (uS/cm)
    SpC_t_act = sum(act * charge * eq_cond) * 1000,
    # Each ion contribution to conductivity
    SpC_cont = act * charge * eq_cond * 1000) %>%
  mutate(SpC_per = SpC_cont / SpC_t_act * 100) # ion cont. as a percent

# Plot of these contributions over a year
# For the labels
df_lab <- ungroup(df_cond) %>%
  filter(month == 12) %>%
  group_by(solute) %>% 
  dplyr::summarize(x = 12,
            y = mean(SpC_per, na.rm = T))

p_cont <- ggplot(data = df_cond,
       aes(x = month,
           y = SpC_per,
           color = solute)) +
  stat_summary() +
  stat_summary(geom = "line") +
  ggrepel::geom_text_repel(data = df_lab, 
            aes(label = solute, x = x, y = y, color = solute), size = 2.5,
            nudge_x = 1,)+ 
  guides(color = "none") +
  scale_color_viridis_d() +
  theme_classic(base_size = 10) +
  theme(plot.tag.position = c(0.15, 0.97)) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  labs(x = "month", 
       y = expression("contribution to "*C[25]~"(%)"))
p_cont

ggsave(plot = p_cont,
       filename = file.path("results", "contribution_to_cond.png"),
       width = 13.5,
       height = 13.5,
       units = "cm",
       dpi = 300)

# Plot ionic strength vs conductivity to get relationship
p_i <- ggplot(data = df_cond,
              aes(x = SpC_t_act,
                  y = I)) +
  geom_point() +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = '~')),
  ) +
  theme_classic(base_size = 10) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  labs(y = "Ionic strength", 
       x = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_i

# This does a very good job at estimating measured SpC
p_cond <- distinct(df_cond, SpC_m, SpC_t, SpC_t_act, date, site_no, temp) %>%
  filter(site_no %in% c("4045900",
                        "4046800",
                        "4048000")) %>%
  ggplot(aes(x = SpC_t_act,
             y = SpC_m)) +
  geom_point() +
  # facet_wrap(~site_no) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(method = "lm") +
  scale_x_continuous(limits = c(100, 350)) +
  scale_y_continuous(limits = c(100, 350)) +
  ggpubr::stat_regline_equation( aes(label =  paste(..eq.label.., 
                                                    ..adj.rr.label.., 
                                                    sep = "~~~~"))) +
  theme_classic(base_size = 10) +
  theme(plot.tag.position = c(0.18, 0.97)) +
  labs(x = expression(C[25]~from~charge~balance~"("*mu*S~cm^{-1}*")"),
       y = expression(measured~C[25]~"("*mu*S~cm^{-1}*")"))
p_cond

ggsave(plot = p_cond, filename = file.path("results", "measured_vs_estimated_SpC.png"),
       dpi = 300,
       units = "cm",
       width = 9.2,
       height = 9.2)

p_cond_all <- p_cond + p_cont + plot_annotation(tag_levels = "a")
p_cond_all
ggsave(plot = p_cond_all, 
       filename = file.path("results", "supplement_conductivity.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 9.2)

# Look at how much of SpC comes from Residual alkalinity and NO3
df_res <- df_cond %>%
  filter(year > 1990) %>%
  mutate(
    # Amount of conductivity from ions in Residual Alkalinity, Alk',
    # Groleau et al. 2000
    # Alk' = [Na+] + [K+] + 2[Mg+] - [Cl-] - 2[SO4--]
    # First get 1 or 0 for those ions
    Alk_res_ion = if_else(solute %in% c("Na", "K", "Mg", "Cl", "SO4"),
                          1,
                          0),
    # Then recalculate with just those ions
    SpC_Alk_res = sum(act * charge * eq_cond * Alk_res_ion) * 1000,
    # Then see what percentage of that contributes to SpC
    Alk_res_cont = SpC_Alk_res / SpC_t_act,
    # Total nitrate contribution
    NO3_ion = if_else(solute == "NO3",
                          1,
                          0),
    SpC_NO3 = sum(act * charge * eq_cond * NO3_ion) * 1000,
    # Total SpC coming not from HCO3 Ca
    epsilon = SpC_Alk_res + SpC_NO3
    )

# Look at mean and sd of residual alkalinity contribution to SpC
a=distinct(df_res, Alk_res_cont, SpC_Alk_res, date, site_no) %>%
  ungroup() %>%
  drop_na() %>%
  filter(year > 1990, between(month, 4,9)) %>%
  summarize(mean_res = mean(SpC_Alk_res),
            sd_res = sd(SpC_Alk_res),
            mean_cont = mean(Alk_res_cont),
            sd_cont = sd(Alk_res_cont))
# In general residual alkalinity contributes 42.3+/-0.1% to SpC
# So, most of the alkalinity is due to Ca and HCO3

# Look at how much ions in Alk' vary over time period (CV)
df_wq %>%
  filter(site_no %in% c("4045900",
                        "4046800",
                        "4048000")) %>%
  filter(year > 1990, between(month, 4,9)) %>%
  # select(Ca, HCO3) %>%
  select(Na, K, Mg, Cl, SO4) %>%
  # filter(solute %in% c("Na", "K", "Mg", "Cl", "SO4")) %>%
  summarize(across(where(is.numeric), ~sd(., na.rm = T)/mean(., na.rm = T)))
  # mutate(CV = sd_uM / mean_uM)
# In general, they vary 20-30%, especially in growing season

# What is the remainder of non HCO3- and Ca- derived SpC on average 
distinct(df_res, date, site_no, epsilon) %>%
  ungroup() %>%
  drop_na(epsilon) %>%
  summarize(mean = mean(epsilon),
            sd = sd(epsilon))
# Epsilon is 110 +/- 22.7 uM
