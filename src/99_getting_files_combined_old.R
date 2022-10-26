df_x <- read_csv(file.path("data", "03_CO2", "DAM_full_correct_daily2.csv")) %>%
  mutate(date = mdy(datetime),
         year = year(date),
         month = month(date)) %>%
  select(-datetime) %>%
  filter(year > 1992)


x = left_join(df_d, df_x)


df_use <- select(df_x,
                 year,
                 month,
                 date,
                 Q_m3s = Discharge,
                 depth_m = Depth,
                 temp_C = `Temp (C)`,
                 alk_molm3 = Alkalinity,
                 cond_uscm = Conductivity,
                 pH,
                 DO_gm3 = Oxy,
                 DO_mmolm3 = `O2 (mmol/m3)`,
                 DOsat = `O2 sat (mmol/m3)`,
                 pCO2_ppmv = `pCO2 (uatm)`,
                 CO2_mmolm3 = `CO2 (mmol/m3)`,
                 dCO2_mmolm3 = `CO2 (mmol/m3) std`,
                 CO2atm_ppmv = CO2_atm,
                 CO2atm_mmolm3 = `CO2_atm (mmol/m3)`,
                 GPP_mmolm2d = )

df_use2 <- df_use %>%
  left_join(select(df_met_err,
                   date,
                   GPP = GPP_mean,
                   dGPP = dGPP_mean,
                   ER = ER_mean,
                   dER = dER_mean,
                   NEP = NEP_mean,
                   dNEP = dNEP_mean))

df_use3 <- df_use2 %>%
  left_join(select(df_CO2,
                   date,
                   K600_met = K600_met_mean,
                   dK600_met = dK600_met_mean,
                   K600_ray = K600_ray_mean,
                   dK600_ray = dK600_ray_mean,
                   Sc_CO2,
                   dSc_CO2,
                   KCO2_met = KCO2_met_mean,
                   dKCO2_met = dKCO2_met_mean,
                   KCO2_ray = KCO2_ray_mean,
                   dKCO2_ray = dKCO2_ray_mean))


df_use4 <- df_use3 %>%
  mutate(ddepth_m = 0, dCO2atm_mmolm3 = 0) %>% #no uncertainty in CO2_atm and depth
  mutate_with_error(CO2flux_ray_mmolm2d ~ depth_m * 
                      (CO2_mmolm3 - CO2atm_mmolm3) * 
                      KCO2_ray) %>%
  mutate_with_error(CO2flux_met_mmolm2d ~ depth_m * 
                      (CO2_mmolm3 - CO2atm_mmolm3) * 
                      KCO2_met)

saveRDS(df_use4,
        file.path("data", "03_CO2", "dampierre_all_daily_data.RDS"))

write_csv(df_use4,
          file.path("data", "03_CO2", "dampierre_all_daily_data.csv"))