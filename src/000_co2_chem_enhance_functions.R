# temperature dependent dissociation constants for carbonic acid
diss_const_fun <- function(TK, h, sal = 0){
  # following Millero et al. 2002 Deep-Sea Research
  # [H+][HCO3-]/[CO2]
  # pK1 = -8.712-9.46*0.001*sal+8.56*0.00001*sal^2+1355.1/TK+1.7976*log(TK)
  # Different calculation for freshater: American Public Health Association 2005  
  # via JAWWA 1990, 82(7) via Plummer and Busenberg 1982
  pK1 = 356.3094 + 0.06091964*TK - 21834.37 / TK - 126.8339 * 
    log10(TK) + 1684915/TK^2
  K1 = 10^-pK1
  # [H+][CO3--]/[HCO3-]
  # pK2 = 17.0001-0.01259*sal-7.9334*0.00001*sal^2+936.291/TK-1.87354*log(TK)-
  #   2.61471*sal/TK+0.07479*sal^2/TK
  # Different calculation for freshwater
  pK2 = 107.8871 + 0.03252849*TK - 5151.79 / TK - 38.92561 * 
    log10(TK) + 563713.9/TK^2
  K2 = 10^-pK2
  tau = 1+h^2/(K1*K2+K1*h) #simplified variable from Hoover and Berkshire 1969
}

# temperature dependent rate constants for CO2 hydration/dehydration
rate_const_fun <- function(TK, h, sal = 0){
  # The combined rate constant, from Hoover and Berkshire 1969 is:
  # r = kCO2 + kOH * [OH-]
  # and [OH-] = Kw / [H+]
  # following Johnson 1982 L&O
  # CO2 + H2O = H2CO3
  lnkCO2 = 1246.98+0*sal-6.19*10^4/TK-183*log(TK)
  kCO2 = exp(lnkCO2) #1/s
  # CO2 + OH- = HCO3-
  # This empirical estimate has the dissociation constant for water built in (Kw)
  lnkOH.Kw = -930.13+0.11*sal+3.1*10^4/TK+140.9*log(TK)
  kOH.Kw = exp(lnkOH.Kw) #mol/L/s
  # Combined rate constant
  r = kCO2 + kOH.Kw / h
}

# temperature dependent CO2 diffusion coefficient (m2/s)
D_fun <- function(TK){
  # following Zeebe 2011 Geochemica et Cosmochimica Acta
  D=0.0000000146836*((TK/217.2056)-1)^1.997
}

# chemical enhancement function based on: 
# temperature (degC)
# pH
# CO2 gas exchange coefficient, KCO2 (1/d)
# depth, d (m)
chem_enh_fun <- function(temp, pH, KCO2, d) {
  TK = temp + 273.15 #temperature in Kelvin
  h = 10^-pH # activity of hydrogen ion
  k_cms = KCO2 * d * 100 / 86400 # CO2 gas exchange constant in cm/s
  # following Hoover and Berkshire 1969 and Wanninkhoff and Knox 1996
  tau = diss_const_fun(TK, h) # dimensionsless
  r = rate_const_fun(TK, h) #1/s
  D = D_fun(TK) * 1e4 # diffusivity of CO2 (cm2/s)
  Q = sqrt(r * tau / D) # 1/cm
  Z = D / k_cms # cm, boundary layer thickness
  alpha = tau / ((tau - 1) + tanh(Q*Z)/(Q*Z)) # dimensionless
}
