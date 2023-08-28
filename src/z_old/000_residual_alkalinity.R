# For most carbonate systems
# Alk = [HCO3-] + 2[CO3--]

# An alternative form is derived from electroneturality (Sigg et al. 1994)
# [Na+] + [K+] + 2[Ca++] + 2[Mg++] = [HCO3-] + 2[CO3--] + [Cl-] + 2[SO4--] + [NO3-]

# Using the electroneutrality equation and general carbonate alkalinity eqn
# Alk = [Na+] + [K+] + 2[Ca++] + 2[Mg++] - [Cl-] - 2[SO4--] - [NO3-]

# Photosynthesis and respiration do not change alkalinity, which is conservative relatively to CO2 variations, whereas DIC can be either taken up by photosynthesis or produced by respiration.The main observable consequence is an increase of pH when photosynthesis prevails while the opposite is observed when respiration dominates

# Calcite precipitation decreases both calcium concentration, alkalinity and DIC. According to the stoichiometry the decrease in alkalinity must be twice the decrease in calcium concentration and DIC. One observable consequence is a decrease in pH

# Therefore, photosynthesis and calcite precipitation have opposite effects on pH variations.In order to verify that changes in the ionic composition are exclusively caused by calcite precipitation or photosynthesis-respiration, the concept ofresidual alkalinity (Alk') can be introduced (Al-Droubiet al., 1980; Michard, 1989).

# Alk'is defined in order to be constant if the main biogeochemical processes are calcite precipitation and photosynthesis-respiration. 

# Therefore, the residual alkalinity canbe expressed as follows [equation (10)]:
# Alk' = [Na+] + [K+] + 2[Mg+] - [Cl-] - 2[SO4--]

# Combining Alk and Alk' eqns
# Alk' = Alk - 2[Ca++] - [NO3-]

df_wq$Alk_res = df_wq$Alkalinity * 1000 - 2 * df_wq$Ca * 1000/ 40.078 - df_wq$NO3 * 1000/ 62.0049
mean(df_wq$Alk_res, na.rm = T)
sd(df_wq$Alk_res, na.rm = T)
