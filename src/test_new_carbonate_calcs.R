# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

df <- readRDS(file.path("data", "03_CO2", "all_hourly_data_complete_updateSpCAT.RDS"))


df2 <- select(df, year, month, date, datetime, solartime, discharge, depth,
              light, temp, pH, AT, SpC, Ca,
              CO2_atm_uatm, CO2_atm_mmolm3, KCO2, KO2,
              NEP_mmolO2m3, enh, DIC_eq_molkg, CO2_eq_molkg,
              O2_mmolm3 = O2_mmol_m3, O2_mgL = O2, exO2_uM = O2ex, O2per = O2per,
)

Ksp_fun <- function (TK = 298.15, sal = 0) {
  # this is the -log(thermodynamic solubility constant) as a function of temp. 
  # in distilled water according to Plummer and Busenberg 1982
  pKsp_0 = -171.9065 - 0.077993 * TK + 2839.319/TK + 71.595 * log10(TK)
  # These are how salinity affects the constant according to 
  # Mucci 1983, Table 7
  B = +(-0.77712 + 0.0028426 * TK + 178.34/TK) * sqrt(sal)
  C = -0.07711 * sal + 0.0041249 * sal^1.5
  log10Kspc = pKsp_0 + B + C
  Kspc = 10^(log10Kspc)
}

df2 <- df2 %>%
  mutate(carb = DIC(temp+273.15, AT, pH, SpC)) %>%
  unnest(carb) %>%
  mutate(SI = log10((10 ^ (4 * -act(temp+273.15, SpC)) * CO3_uM / 1E6) * 
                      (10 ^ (4 * -act(temp+273.15, SpC)) * Ca / 40.078 / 1000) /
                      Ksp_fun(temp+273.15)))

tail(df2$SI, n = 10) - tail(df$LSI, n = 10)  
