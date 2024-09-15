# -------------------------------------
# Author: Jake Diamond
# Purpose: Functions to solve the carbonate system for use in model
# Date: 2023-04-10
# -------------------------------------

# Helper functions --------------------------------------------------------
# Photosynthetically active radiation function (umol/m2/s)
par_fun <- function(hour, day = 180, latitude = 0, max.insolation = 2326) {
  declin <- (23.45 * sin((2 * pi / 365) * (284 + day))) * pi / 180
  hour_angle <- (360/24) * (hour - 12) * pi / 180
  lat <- latitude * pi / 180
  zenith <- acos(sin(lat) * sin(declin) +
                   cos(lat) * cos(declin) * cos(hour_angle))
  insolation <- max.insolation * cos(zenith)
  insolation <- pmax(insolation, 0)
}

# Schmidt number function
Sc <- function(gas = "CO2", temp = 25) {
  # From Table 1 in Wanninkhof L&O:Methods 2014
  co <- switch(
    gas,
    "CO2" = c(1923.6, -125.06, 4.3773, -0.085681, 0.00070284),
    "O2" = c(1745.1, -124.34, 4.8055, -0.10115, 0.00086842)
  )
  co[1] + co[2] * temp + co[3] * temp^2 + co[4] * temp^3 + co[5] * temp^4
}

# -log10(activity coeff.), gamma, function
# TK is temperature in Kelvin
# cond is specific conductance in uS/cm
# I is ionic strength in mol/L, estimated from cond unless known
act <- function(TK, cond, I = NULL) {
  if(is.null(I)) {
    I <- 1.6 * 10^-5 * cond
  }
  E <- 60954 / (TK + 116) - 68.937 #dielectric constant
  # Davies approximation for activity coefficient
  gamma <- 1.82 * 10^6 * (E * TK)^-1.5 * ((sqrt(I) / (1 + sqrt(I))) - 0.3 * I)
  return(gamma)
}

# Calculate first freshwater thermodynamic equilibrium constant
# [H+][HCO3-]/[CO2]
# at infinite dilution via Plummer and Busenberg 1982
K1 <- function(TK) {
  pK1 <- 356.3094 + 0.06091964 * TK - 21834.37 / TK - 126.8339 *
    log10(TK) + 1684915 / TK^2
  K1 <- 10^-pK1
  return(K1)
}

# Calculate second freshwater thermodynamic equilibrium constant
# [[H+][CO3--]/[HCO3-]
# at infinite dilution via Plummer and Busenberg 1982
K2 <- function(TK){
  pK2 <- 107.8871 + 0.03252849 * TK - 5151.79 / TK - 38.92561 *
    log10(TK) + 563713.9 / TK^2
  K2 <- 10^-pK2
  return(K2)
}

# Kw function, thermodynamic dissociation constant for water
Kw <- function(TK) {
  # Kw is dissociation constant for water
  pKw <- 4471 / TK + 0.01706 * TK - 6.0875
  Kw <- 10^-pKw
  return(Kw)
}

# Function for Henry's constant at given temperature and pressure (mol/kg/atm)
KH <- function(TK = 298, press = 1, alt = NULL) {
  # TK = temperature in Kelvin,
  # press = pressure in atm
  # If altitude is given, calculate pressure
  if(!is.null(alt)) { press <- 1.01325 * (1 - 2.25577E-5 * alt)^5.25588 }
  
  # using version from Plummer and Busenberg (1982)
  k0 <- 10^(108.3865 + 0.01985076 * TK - 6919.53 / TK -
              40.45154 * log10(TK) + 669365 / TK^2)
  # Correct for pressure; Weiss 1974 Marine Chemistry
  R <- 82.05736 # Gas constant, [cm3 * atm / (K * mol)]
  vbar <- 32.3 #partial molal volume of CO2 in solution (cm^3/mol)
  kh <- k0 * exp(((1 - press) * vbar) / (R * TK))
  return(kh)
}

# Apparent, or stoichiometric, solubility constant for calcite
# Mucci 1983
# Ca2++ + CO3-- = CaCO3 (s)
Ksp <- function(TK = 298.15, sal = 0) {
  # this is the -log(thermodynamic solubility constant) as a function of temp.
  # in distilled water according to Plummer and Busenberg 1982
  pKsp_0 <- -171.9065 - 0.077993 * TK + 2839.319 / TK + 71.595 * log10(TK)
  # These are how salinity affects the constant according to
  # Mucci 1983, Table 7
  b <- +(-0.77712 + 0.0028426 * TK + 178.34 / TK) * sqrt(sal)
  c <- -0.07711 * sal + 0.0041249 * sal^1.5
  log10Kspc <- pKsp_0 + b + c
  Kspc <- 10^(log10Kspc)
  return(Kspc)
}

# Density of water function kg/m3
rhow <- function(temp = 25) {
  # Martin and McCutcheon (1999)
  999.842594 + 6.793952 * 10^(-2) * temp - 9.095290 * 10^(-3) * temp^2 +
    1.001685 * 10^(-4) * temp^3 - 1.120083 * 10^(-6) * temp^4 +
    6.536335e-9 * temp^5
}

# Carbonate system function -----------------------------------------------
# Calculate carbonate system based on alkalinity and pH
# Or based on alkalinity and DIC
# TK = temperature in Kelvin
# AT = total alkalinity (mol/m3)
# pH = pH on free scale, initial guess for solving ALK/DIC system is 8
# cond = specific conductance at 25degC (uS/cm)
carb <- function(TK, AT, pH = 8, cond, TC = NULL) {
  
  # Calculate all the parameters
  AT <- AT / rhow(TK - 273.15) # mol/kg
  gamma <- act(TK, cond) # -log(monovalent act. coef)
  Kh <- KH(TK) #henry's constant uncorrected for salinity
  I <- 1.6 * 10^-5 * cond # ionic strength, cond in uS/cm
  S <- 53.974 * I #salinity from ionic strength, estimated
  aH <- 10 ^ -gamma # activity coefficient for H+
  aOH <- 10 ^ -gamma # activity coefficient for OH-
  aHCO3 <- 10 ^ -gamma # activity coefficient for HCO3-
  aCO3 <- 10 ^ (4 * -gamma) # activity coefficient for CO32-
  KWa <- Kw(TK) / (aH * aOH) # apparent dissociation coefficient
  K1a <- K1(TK) / (aH * aHCO3) # apparent dissociation coefficient
  K2a <- K2(TK) / (aH * aCO3 / aHCO3) # apparent dissociation coefficient
  KHa <- Kh + (0.023517 - 0.023656 * TK / 100 + 0.0047036 * TK / 100 * TK / 100) * S #apparent Henry constant, Weiss 1974
  
  # Calculate pH if DIC is given
  if(!is.null(TC)) {
    # Simple iterative scheme from Park 1969 in L&O
    TC <- TC / rhow(TK - 273.15) # mol/kg
    # Iterate for H and CA by repeated solution
    H <- 10 ^ (-pH)  # initial guess from arg list      
    delH <- H     
    tol <- 1.e-15 
    
    iter <- 0
    
    while (delH > tol) {     # iterate until H converges
      
      H_old <- H  # previous H
      
      # solve for carbonate alkalinity from TA
      CA <- AT  
      
      # solve quadratic for H
      a <- CA
      b <- K1a * (CA - TC)
      c <- K1a * K2a * (CA - 2 * TC)
      H <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)  
      
      # How different is new estimate from previous one?
      delH <- abs(H - H_old)
      iter <- iter + 1
      
    }
    pH <- -log10(H)
    CT <- TC
  } else {
    H <- 10^(-pH) #hydrogen ion conc.
  }
  
  OH <- KWa / H # hydroxide ion conc.
  
  # Solve the carbonate system, from Stumm and Morgan, (mol/kg)
  alpha1 <- (H * K1a) / (H^2 + K1a * H + K1a * K2a)  # HCO3 ionization fraction
  alpha2 <- (K1a * K2a) / (H^2 + K1a * H + K1a * K2a) # CO3 ionization fraction
  if(is.null(TC)) {
    CT  <- (AT - OH + H) / (alpha1 + 2 * alpha2) #total carbon, DIC
  }
  CO2 <- CT * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a)
  HCO3 <- CO2 * K1a / H
  CO3 <- HCO3 * K2a / H
  
  # Uses the apparent Henry's Law constant and converts from atm to uatm
  pCO2 <- (CO2 / KHa) * 1000000
  
  # Get concentrations in mol/m3
  CO2 <- CO2 * rhow(TK - 273.15)
  HCO3 <- HCO3 * rhow(TK - 273.15)
  CO3 <- CO3 * rhow(TK - 273.15)
  DIC <- CT * rhow(TK - 273.15)
  
  # Return the carbonate system
  data.frame(pH = pH, CO2 = CO2, pCO2 = pCO2, DIC = DIC, HCO3 = HCO3, CO3 = CO3)
}


pH_TA_CO2 <- function(TK, cond, TAi, CO2i){
  # Calculate all the parameters
  gamma = act(TK, cond) # -log(monovalent act. coef)
  aH = 10 ^ -gamma # activity coefficient for H+
  aOH = 10 ^ -gamma # activity coefficient for OH-
  aHCO3 = 10 ^ -gamma # activity coefficient for HCO3-
  aCO3 = 10 ^ (4 * -gamma) # activity coefficient for CO32-
  KWa = Kw(TK) / (aH * aOH) # apparent dissociation coefficient
  K1a = K1(TK) / (aH * aHCO3) # apparent dissociation coefficient
  K2a = K2(TK) / (aH * aCO3 / aHCO3) # apparent dissociation coefficient
  pHGuess    = 8      # this is the first guess
  pHTol      = 0.0001 # tolerance
  ln10       = log(10)
  pH         = pHGuess
  deltapH    = pHTol + pH
  while (abs(deltapH) > pHTol) {
    H         = 10^(-pH)
    HCO3      = K1a * CO2i / H
    CO3       = K1a * K2a * CO2i / (H * H)
    CAlk      = HCO3 + 2 * CO3
    OH        = KWa / H
    Residual  = TAi - CAlk - OH + H
    #  find Slope dTA/dpH
    #  (this is not exact, but keeps all important terms):
    Slope     = ln10 * (HCO3 + 4*CO3 + OH + H)
    deltapH   = Residual / Slope # this is Newton's method
    pH = pH + deltapH
  }
  return(pH)
}
