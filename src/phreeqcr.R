# Define the solution composition
pH <- as.numeric(df[283651, "pH"])
Alkalinity <- as.numeric(df[283651, "Alk_molkg"])
Ca <- as.numeric(df[283651, "Ca"]) / 40.078 / 1000
CO2 <- as.numeric(df[283651, "CO2_uM"]) / 1E6
HCO3 <- as.numeric(df[283651, "HCO3"]) / 1E6
temp <- as.numeric(df[283651, "temp"])
# Mg <- 0.001
# Cl <- 0.001
# Na <- 0.001
# K <- 0.001

# Define the input file text
input_text <- paste0(
  "DATABASE phreeqc.dat\n",  # contains the PWP calcite rate
  # "RATES",
  
  "SOLUTION 1 Calcite Dissolution\n",
  "  units mol/l\n",
  "  density 1.00\n",
  "  temp      ", temp, "\n",
  "  pH      ", pH, "\n",
  "  Alkalinity      ", Alkalinity, "\n",
  "  CO2      ", CO2, "\n",
  "  HCO3      ", HCO3, "\n",
  "  Ca      ", Ca, "\n",
  # "  Mg      ", Mg, "\n",
  # "  Cl      ", Cl, "\n",
  # "  Na      ", Na, "\n",
  # "  K      ", K, "\n",
  "RATES \n",
  "Calcit2\n", # simplified rate...
  "-start\n",
  "10 rate = 10^-6.91 - 10^-1.52 * (tot(", Ca, "))^2\n",  # mmol/cm2/s
  "20 save 1 * rate * time\n",                         # integrate in moles/L
  "-end\n",
  "EQUILIBRIUM_PHASES 1\n",
  "  Calcite 0\n",        #Calculates CCPP by setting solution in eq with solid
  # "REACTION_TEMPERATURE\n",
  # "  25\n",
  # "KINETICS 1\n",
  # "  Calcit2\n",
  # "  -m  7.0E-4\n",
  # "  -m0  7.0E-4\n",
  # "  -parms 5.0 0.3\n",
  # "  -step 4 in 1.0 hour\n",
  # "  -tol 1E-10\n",
  # "  -initial 0.0\n",
  # "  -surface 1.0\n",
  "SELECTED_OUTPUT\n",
  "  -file output.txt\n",
  # "  -totals     Ca\n",
  "  -time true\n",
  "  -solution true\n",
  # "  -reaction true\n",
  # "  -activities Ca+2 CO2 HCO3- CO3-2\n",
  # "  -alkalinity true\n",
  # "  -kinetic_reactants CaCO3\n",
  # "  -step false\n",
  # "  -molalities\n",
  # "  -si Calcite\n",
  # "  -kinetics true\n",
  # "  -rate Calcite\n",
  "END\n"
)

phrRunString(input_text)
phrGetSelectedOutput()


# Write the input file to disk
writeLines(input_text, file.path("data", "06_PHREEQC", "input.txt"))

# phreeqc(input = "input_file.txt", output = "output_file.txt")

phrRunFile(file.path("data", "06_PHREEQC", "input.txt"))

# phrRunString(input_text)
phrRunString(input_text)
phrGetSelectedOutput()


phrRunFile(file.path("data", "06_PHREEQC", "calc_test_rate.txt"))
x = phrGetSelectedOutput()
summary(x$n1)
?phrGetSelectedOutput
