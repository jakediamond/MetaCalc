calcite_rxn <- function(LSI, pH, Ca, temp) {
  
  # Convert temperature to Kelvin
  T_K <- temp + 273.15
  
  # Calculate the equilibrium constant for calcite dissolution or precipitation
  # using the temperature-dependent equation from Plummer et al. (1978)
  ln_K_eq <- -171.9065 - 0.077993*T_K + 2839.319/T_K + 71.595*log(T_K)
  K_eq <- exp(ln_K_eq)
  
  I <- 1.6*10^-5 * 300
  
  # Calculate the calcium ion activity coefficient using the Davies equation
  A_Ca <- 10^(-0.5092*(2^2)*(sqrt(I)/(1+sqrt(I))) - 0.0018*(T-298.15))
  
  # Calculate the pH-dependent term in the LSI equation
  pH_term <- (10^-pH)/(10^-8.3)
  
  # Calculate the LSI
  LSI_calc <- log10((Ca*A_Ca*K_eq)/pH_term)
  
  # Calculate the amount of calcite dissolution or precipitation in mol/m3
  if (LSI > LSI_calc) {
    # Precipitation
    delta_CaCO3 <- (10^(-LSI) - 10^(-LSI_calc))*pH_term/(A_Ca*K_eq)
  } else {
    # Dissolution
    delta_CaCO3 <- (10^(-LSI_calc) - 10^(-LSI))*pH_term/(A_Ca*K_eq)
  }
  
  return(delta_CaCO3)
}
calcite_rxn(0.4, 7, 30, 25)


calcite_calc <- function(LSI, pH, Ca, SC, temp) {
  # Calculate saturation concentration of calcite
  Ksp <- 10^(-8.48 + (0.0118 * temp) + (0.0028 * (25 - temp)) - (0.000026 * temp^2))
  sat_conc <- Ksp / ((Ca * 10^-3) * exp(-0.00286 * (25 - temp)))
  
  # Calculate ion activity product
  IAP <- (10^(-pH) * (Ca * 10^-3) * 10^(SC-7))
  
  # Calculate precipitation or dissolution
  if (IAP > sat_conc) {
    return(-(IAP - sat_conc))
  } else {
    return(sat_conc - IAP)
  }
}


calcite_precipitation_dissolution <- function(LSI, pH, Ca_conc, temp) {
  
  # Convert temperature to Kelvin
  temp_K <- temp + 273.15
  
  # Calculate the saturation index (SI)
  SI <- log10(Ca_conc*10^-3.3) + (9.76 + (0.0118*temp_K) - (0.000116*temp_K^2))*pH - 13.12
  
  # Calculate the equilibrium concentration of calcium carbonate in mol/m3
  CaCO3_eq <- 10^(SI - LSI)
  
  # Determine if there is precipitation or dissolution
  if (SI > LSI) {
    result <- paste0("Calcite precipitation of ", format(CaCO3_eq, scientific = FALSE), " mol/m3")
  } else if (SI < LSI) {
    result <- paste0("Calcite dissolution of ", format(-CaCO3_eq, scientific = FALSE), " mol/m3")
  } else {
    result <- "No precipitation or dissolution"
  }
  
  return(result)
}
