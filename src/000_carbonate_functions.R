# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate carbonate system
# Date: 2023-02-01
# -------------------------------------

# Density of water function kg/m3
rhow <- function (temp = 25) {
  # Martin and McCutcheon (1999)
  999.842594 + 6.793952*10^(-2)*temp - 9.095290*10^(-3)*temp^2 + 
    1.001685*10^(-4)*temp^3 - 1.120083*10^(-6)*temp^4 + 6.536335e-9*temp^5
}

# Schmidt number function
Sc <- function (gas = "CO2", temp = 25) {
  # From Table 1 in Wanninkhof L&O:Methods 2014
  co <- switch(gas,
               "CO2" =  c(1923.6, -125.06, 4.3773, -0.085681, 0.00070284),
               "O2" = c(1745.1, -124.34, 4.8055, -0.10115, 0.00086842)
  )
  co[1] + co[2] * temp + co[3] * temp^2 + co[4] * temp^3 + co[5] * temp^4
}

# Function for Henry's constant at given temperature and pressure (mol/kg/atm) Weiss 1974
henry <- function (TK = 298, S = 0) {
  KH <- exp(-60.2409 + 93.4517 * (100 / TK) + 23.3585 * log(TK / 100) +
              S * (0.023517 - 0.023656 * (TK/100) + 0.0047036 * (TK/100)^2))
}
# Davies equation, pfm is -log(activity coefficient) for monovalent species
act <- function(TK, cond){
  I <- 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  E <- 60954/(TK+116) - 68.937 #dielectric constant
  # Davies approximation for activity coefficient
  gamma <- 1.82 * 10^6 * (E*TK)^-1.5 * ((sqrt(I)/(1+sqrt(I))) - 0.3 * I)
}

# Kw function, dissociation constant for water; The Physical Chemistry of Electrolytic Solutions
Kw <- function(TK) {
  pKw <- 4471 / TK + 0.01706 * TK -6.0875
  Kw <- 10^-pKw
}

# K1 function, first dissociation constant
K1 <- function(TK){
  # [H+][HCO3-]/[CO2]
  # Plummer and Busenberg 1982
  pK1 <- 356.3094 + 0.06091964*TK - 21834.37 / TK - 126.8339 * 
    log10(TK) + 1684915/TK^2
  K1 <- 10^-pK1
}

# K2 function, second dissociation constant
K2 <- function(TK){
  # [H+][CO3--]/[HCO3-]
  # via Plummer and Busenberg 1982
  pK2 <- 107.8871 + 0.03252849*TK - 5151.79 / TK - 38.92561 * 
    log10(TK) + 563713.9/TK^2
  K2 <- 10^-pK2
}

# Calculate DIC
DIC <- function(TK, AT, pH, cond) {
  AT <- AT / rhow(TK-273.15) # mol/kg
  gamma <- act(TK, cond) # -log(monovalent act. coef)
  I <- 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  S <- 53.974 * I #salinity from ionic strength, estimated
  ACH <- 10 ^ -gamma # activity coefficient for H+
  ACOH <- 10 ^ -gamma # activity coefficient for OH-
  ACHCO3 <- 10 ^ -gamma # activity coefficient for HCO3-
  ACCO3 <- 10 ^ (4 * -gamma) # activity coefficient for CO32-
  KH <- henry(TK, S) #henry's constant corrected for salinity (mol/kg-atm)
  K1 <- K1(TK)
  K2 <- K2(TK)
  KWa <- Kw(TK) / (ACH * ACOH) # apparent dissociation coefficient 
  K1a <- K1(TK) / (ACH * ACHCO3) # apparent dissociation coefficient 
  K2a <- K2(TK) / (ACH * ACCO3 / ACHCO3) # apparent dissociation coefficient 
  KHa <- KH * exp(S * (0.023517 - 0.023656 * (TK/100) + 0.0047036 * (TK/100)^2)) #apparent Henry constant, Weiss 1974
  H <- 10^(-pH) # hydrogen ion concentration
  OH <- KWa / H # hydroxide ion concentration
  alpha1 <- (H * K1a)/(H^2 + K1a * H + K1a * K2a)  # HCO3 ionization fraction 
  alpha2 <- (K1a * K2a) / (H^2 + K1a * H + K1a * K2a) # CO3 ionization fraction
  CT  <- (AT- OH + H) / (alpha1 + 2*alpha2) #total carbon, DIC (mol/kg)
  CO2 <- CT * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a) #CO2 concentration (mol/kg)
  HCO3 <- CT * alpha1 #HCO3 concentration (mol/kg)
  CO3 <- CT * alpha2 #CO3 concentration (mol/kg)
  # Get concentrations in uM
  CO2_uM <- CO2 * rhow(TK - 273.15) * 1000
  HCO3_uM <- HCO3 * rhow(TK - 273.15) * 1000
  CO3_uM <- CO3 * rhow(TK - 273.15) * 1000
  DIC_uM <- CT * rhow(TK - 273.15) * 1000

  #Uses Henry's Law constant and converts from atm to uatm 
  pCO2 <- (CO2 / KHa) * 1000000
  
  # Result as dataframe of carbonate species
  data.frame(CO2_uM = CO2_uM, pCO2_uatm = pCO2,
             DIC_uM = DIC_uM, HCO3_uM = HCO3_uM, CO3_uM = CO3_uM)
}

# Apparent, or stoichiometric, solubility constant for calcite
# Mucci 1983
# Ca2++ + CO3-- = CaCO3 (s)
Ksp <- function (TK = 298.15, sal = 0) {
  # this is the -log(thermodynamic solubility constant) as a function of temp. 
  # in distilled water according to Plummer and Busenberg 1982
  pKsp_0 <- -171.9065 - 0.077993 * TK + 2839.319/TK + 71.595 * log10(TK)
  # These are how salinity affects the constant according to 
  # Mucci 1983, Table 7
  B <- +(-0.77712 + 0.0028426 * TK + 178.34/TK) * sqrt(sal)
  C <- -0.07711 * sal + 0.0041249 * sal^1.5
  log10Kspc <- pKsp_0 + B + C
  Kspc <- 10^(log10Kspc)
}

# Calcite saturation index function
SI <- function(TK, cond, Ca, CO3){
  rho_w <- rhow(TK - 273.15) #density of water
  gamma <- act(TK, cond) # -log(monovalent act. coef)
  SI <- log10(((10^(4 * -gamma) * Ca / rho_w)) * (10^(4 * -gamma) * CO3 / rho_w ) / 
          Ksp(TK))
}
