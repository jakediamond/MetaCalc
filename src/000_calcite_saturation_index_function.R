# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate the calcite saturation index 
# Date: 2023 February 2
# -------------------------------------

# Mucci 1983
# Ca2++ + CO3-- = CaCO3 (s)
# Apparent, or stoichiometric, solubility constant for calcite
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

# pK1 function 
pK1_fun <- function(TK, sal = 0){
  # following Millero et al. 2002 Deep-Sea Research
  # [H+][HCO3-]/[CO2]
  # pK1 = -8.712-9.46*0.001*sal+8.56*0.00001*sal^2+1355.1/TK+1.7976*log(TK)
  # Different calculation for freshater: American Public Health Association 2005  
  # via JAWWA 1990, 82(7) via Plummer and Busenberg 1982
  pK1 = 356.3094 + 0.06091964*TK - 21834.37 / TK - 126.8339 * 
    log10(TK) + 1684915/TK^2
}

# pK2 function
pK2_fun <- function(TK, sal = 0){
  # [H+][CO3--]/[HCO3-]
  # pK2 = 17.0001-0.01259*sal-7.9334*0.00001*sal^2+936.291/TK-1.87354*log(TK)-
  #   2.61471*sal/TK+0.07479*sal^2/TK
  # Different calculation for freshwater
  # via Plummer and Busenberg 1982
  pK2 = 107.8871 + 0.03252849*TK - 5151.79 / TK - 38.92561 * 
    log10(TK) + 563713.9/TK^2
}

# pKw function
pKw_fun <- function(TK) {
  # Kw is dissociation constant for water
  pKw = 4471 / TK + 0.01706 * TK -6.0875
}

# pfm function, fm is activity coefficient for monovalent species
pfm_fun <- function(TK, cond){
  I = 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  E = 60954/(TK+116) - 68.937 #dielectric constant
  # Davies approximation for activity coefficient
  pfm = 1.82 * 10^6 * (E*TK)^-1.5 * ((sqrt(I)/(1+sqrt(I))) - 0.3 * I)
}

# Calcite saturation function. This predicts whether water has a tendency to 
# precipitate or dissolved CaCO3. This version is the Langelier Saturation Index,
# so if it is greater than 0, precipitation is likely, and if less than 0
# dissolution is likely
SI_cal_fun <- function(temp = 25, #degC
                       pH = 7, 
                       cond = 250, #uS/cm
                       Ca = 30, #mg/L
                       HCO3 = 0.0015, #mol/L
                       sal = 0) #ppt
  {
  # American Public Health Association 2005
  # via JAWWA 1990, 82(7)
  # SI_calcite = pH - pHs; where pHs is the pH of the water when in equilibrium 
  # with calcium carbonate at the existing [Ca2++] and [HCO3-]
  # pHs = pK2 - pKs + p[Ca2++] + p[HCO3-] + 5pf_m
  # Where "p" refers to the -log10()
  
  # Some quick unit changes
  MW_Ca = 40.078 #g/mol
  Ca_mol = Ca / 1000 / MW_Ca #mol/L
  pCa = -log10(Ca_mol)
  pHCO3 = -log10(HCO3)
  TK = temp + 273.15 # temperature in kelvin
  
  # K2 is second dissociation constant for carbonic acid
  pK2 = pK2_fun(TK)
  
  # Ks is solubility product constant for CaCO3
  pKs = -log10(Ksp_fun(TK, sal))
  
  # Kw is dissociation constant for water
  pKw = pKw_fun(TK)
  
  # fm is activity coefficient for monovalent species
  pfm = pfm_fun(TK, cond)
  
  # final calculations
  pHs = pK2 - pKs + pCa + pHCO3 + 5*pfm
  SI = pH - pHs
}


# function for calcium carbonate precipitation potential (CCPP) [mg/L]
# This determines the actual amount of calcium carbonate that precipitates
# given ideal nucleation conditions and assuming equilibrium is achieved
# and can be used then to determine how much
# CO2 is released from this precipitation reaction
# Ca2++ + 2HCO3- <--> CaCO3 + CO2 + H2O
# from JAWWA 1990, 82(7) via Plummer and Busenberg 1982
# Only consideres HCO3-, CO3--, OH-, and H+ in calculating alkalinity, and 
# does not consider ion pairs. Therefore it tends to overestimate CaCO3 that can
# be precipitated and underestimates the amount that can be dissolved
CCPP_fun <- function(temp, pH, cond, alk, Ca){
  # First determine properties of water
  TK = temp + 273.15 # temperature in kelvin
  
  # monovalent activity coefficient
  pfm = pfm_fun(TK, cond)
  
  # Initial alkalinity (eq/L) from mmol/L
  Alk_i = alk / 1000

  # Initial calcium concentration
  MW_Ca = 40.078 #g/mol
  Ca_i = Ca / 1000 / MW_Ca #mol/L
  
  # initial hydrogen ion concentration (mol/L)
  H_i = 10^(pfm - pH)
  
  # get pK for equilibrium
  pK1 = pK1_fun(TK)
  pK2 = pK2_fun(TK)
  pKw = pKw_fun(TK)
  pKs = -log10(Ksp_fun(TK))
  
  # get apparent equilibrum constants
  K1 = 10^(2*pfm - pK1)
  K2 = 10^(4*pfm - pK2)
  Kw = 10^(2*pfm - pKw)
  Ks = 10^(8*pfm - pKs)
  
  # determine values for parameters, p, s, and t for initial conditions
  p = (2 * H_i + K1) / K1
  s = H_i - (Kw / H_i)
  t = (2 * K2 + H_i)/ H_i
  
  # calculate initial acidity
  Acy_i = ((Alk_i + s) / t) * p + s
  
  # Use iterative calculations to determine pH when equilibrium is established
  # Assume a pH for the equilibrated condition (pHeq) and calculate the [H+]
  # Start with pHeq = 7
  pH_start = 7
  eq_fun = function(pH_eq) {
    # Get equilbrium [H+]
    H_eq = 10^(pfm - pH_eq)
    # Get equilibrium parameters
    p_eq = (2 * H_eq + K1) / K1
    r_eq = (H_eq + 2 * K2) / K2
    s_eq = H_eq - (Kw / H_eq)
    t_eq = (2 * K2 + H_eq)/ H_eq
    # get RHS and LHS of eqn from Plummer and Busenberg 1982
    RHS = 2 * Ks * r_eq * p_eq / (t_eq * (Acy_i - s_eq)) - (t_eq * (Acy_i - s_eq)/
                                                              p_eq) + s_eq
    LHS = 2 * Ca_i - Alk_i
    # Minimization function, minimize sum of squared errors
    sum(RHS - LHS)^2
  }
  # optimize to find the equilibrium pH
  pH_eq = optim(par = pH_start, fn = eq_fun)$par
  # Get equilbrium [H+]
  H_eq = 10^(pfm - pH_eq)
  # Get equilibrium parameters
  p_eq = (2 * H_eq + K1) / K1
  r_eq = (H_eq + 2 * K2) / K2
  s_eq = H_eq - (Kw / H_eq)
  t_eq = (2 * K2 + H_eq)/ H_eq
  # This difference must be positive
  # if (Acy_i - s_eq < 0)
  #   warning("not a valid solution")
  
  # Determine equilibrium alkalinity
  Alk_eq = (t_eq / p_eq) * (Acy_i - s_eq) - s_eq
  
  # Determine CCPP (mg/L)
  CCPP = 50000 * (Alk_i - Alk_eq)
}


# Calcite precipitation kinetics 
# Plummer et al., 1978, AJS 278, 179; Appelo et al., AG 13, 257.
calc_ppt_fun <- function(){
  # Calcite mmol/cm2/s
  # -m0    3.e-3
  # -m     3.e-3

  sa_calc = 1.67e-5 # specific surface area of calcite, cm^2/mol calcite
  exp_const = 0.6 # exponent for M/M0
  
  si_cc = SI_cal_fun()
  
  if (M <= 0 & si_cc < 0){
    moles
  } 
  k1 = 10^(0.198 - 444.0 / TK )
  k2 = 10^(2.84 - 2177.0 /TK )
  if (TC <= 25){
    k3 = 10^(-5.86 - 317.0 / TK)
  } else {
    k3 = 10^(-1.1 - 1737.0 / TK )
  }
  if (M0 > 0) {
    area = sa_calc * M0 * (M / M0)^exp_const
  } else {
    area = sa_calc * M
  } 
  rate = area * (k1 * ACT("H+") + k2 * ACT("CO2") + k3 * ACT("H2O"))
  rate = rate * (1 - 10^(2/3*si_cc))
  moles = rate * 0.001 * TIME # convert from mmol to mol
  
}


# pK2_mil = 17.00012+936.291/temp-1.87354*log(temp)
# -log10(Kspa)
# log10(Oc)
# 
# # Milero 1995
# 0.023/1.6E-5
# 
# 
# 
# # following Millero 1995 Geochemica et Cosmochimica Acta,
# # which follows Roy et al. 1993 for low salinities (<5ppt)
# lnK1 = 290.9097 - 14554.21/temp - 45.0575*log(temp) +
#   (-228.39774 + 9714.6839/temp + 34.485796*log(temp)) * sqrt(sal) +
#   (54.20871 - 2310.48919/temp - 8.19515 * log(temp)) * sal +
#   (-3.969101 + 170.22169/temp + 0.603627 * log(temp)) * sal ^1.5 -
#   0.00258768 * sal^2
# 
# K1 = exp(lnK1)
# sal = 0
# temp = 293.15
# lnK2 = 207.6548 - 11843.79/temp - 33.6485*log(temp) +
#   (-167.69908 + 6551.35253/temp + 25.928788*log(temp)) * sqrt(sal) +
#   (39.75854 - 1566.13883/temp - 6.171951 * log(temp)) * sal +
#   (-2.892532 + 116.270079/temp + 0.45788501 * log(temp)) * sal ^1.5 -
#   0.00613142 * sal^2
# 
# K2 = exp(lnK2)
# pK2_2 = -log(K2)
# log(1)
