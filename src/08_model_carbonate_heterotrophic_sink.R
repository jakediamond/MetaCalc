# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

# Load libraries
library(plotly)
library(htmltools)
library(lubridate)
library(seacarb)
library(tidyverse)

# Load data---------------------------------------------------------------
# Hourly carbonate system
df <- readRDS(file.path("data", "03_CO2", "hourly_nep_co2_data.RDS"))

# CO2 function not considering carbonate system
# co2_fun <- function(co2_mod, nep, k_co2_meta_1_d, atmosphere_co2_mmol_m3) {
#   if (is.na(co2_mod[1])) {
#     co2_mod[1] - NEP[2] - k_co2_meta_1_d[2] * (co2_mod[2] - atmosphere_co2_mmol_m3[2])}
#   else {co2_mod[1]}
# }
# 
# fill_in2 <- function(prev, new) if (is.na(new[1])) prev[1]*(1+new[2]) else new[1]
# options(pillar.sigfig = 5)
# x %>%
#   mutate(b = accumulate(pmap(list(b, growth), c)[-1], .init = b[1], fill_in2))

# Get period just before heterotrophic sink behavior
df_test <- filter(df, month == 1, year == 1994) %>%
  mutate(co2_mod = if_else(row_number() == 1, co2_mmol_m3, NA_real_))
#   mutate(
#     # Initiate model
#     co2_mod = accumulate(pmap(list(b, growth), c)[-1], .init = co2_mmol_m3[1],
#                          co2_fun),
#     # difference between model and model with pH, alk, temp
#     co2_dif= co2_mod - co2_mmol_m3
# )
# str(df_test)

for (i in 1:nrow(df_test)) {
  if(is.na(df_test$co2_mod[i])) {
    df_test$co2_mod[i] <- df_test$co2_mod[i-1] - df_test$NEP[i-1]*1000/32/df_test$depth[i-1] -
      df_test$k_co2_meta_1_d[i-1] / 24 * (df_test$co2_mod[i-1] - df_test$atmosphere_co2_mmol_m3[i-1])}
  }



ggplot(data = df_test[1:(24*30), ],
       aes(x = datetime)) +
  geom_path(aes(y = co2_mmol_m3), color = "black") +
  geom_path(aes(y = co2_mod), color = "red") +
  geom_path(aes(y = atmosphere_co2_mmol_m3), color = "blue") +
  theme_bw() +
  labs(x = "",
       y = expression(CO[2]~"("*mmol~m^{-3}*")"))

ggsave(filename = file.path("results", "2012_april_model_co2_no_carbonate.png"),
       dpi = 600,
       width = 16,
       height = 16,
       units = "cm")

# Define model parameters --------------------------------------------------------
parms <- c(

)

# Overall model function ----------------------------------------------------------
model <- function(time, state, parms){
  with(as.list(parms), {
    
    # Unpack states
    CO2 <- state[1:N]
    
    
    
    # Rates of change (g O2/m3)
    dDO = (adv_dis + storage_str + (gpp + er + reaeration) / d) * del_t
    dDO_stor = storage * del_t
    diffs = c(dDO = dDO, dDO_stor = dDO_stor)
    
    # Fluxes
    fluxes = c(GPP = gpp * del_t, ER = er * del_t)
    
    # Output of model
    return(list(diffs, fluxes))
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
DO_ini = 9 # mg/L
DO_stor_ini = 9 # mg/L
yini <- c(CO2 = rep(DO_ini, with(as.list(parms), L / dx)),
          DO_stor = rep(DO_stor_ini, with(as.list(parms), L / dx)))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 1 # time step (h)
days <- 30 # number of days to simulate
simulation_time <- 24 * days / del_t # simulation time (h)
times <- seq(0, simulation_time, by = del_t)

# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini,
              times = times,
              func = model,
              parms = parms,
              nspec = 2)

# Examine the model output ------------------------------------------------
# Reorganize the data into long-form
df <- as_tibble(out) %>%
  gather(key, value, -time) %>%
  separate(key, c("type", "reach"), "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(time_hr = time * del_t,
         dist = as.numeric(reach) * with(as.list(parms), dx)) %>%
  mutate(type_plot = recode(type,
                            `DO` = "DO~(mg~L^{-1})",
                            `DO_stor` = "DO[storage]~(mg~L^{-1})",
                            `ER.DO_stor` = "ER~(g~O^2~m^{-2}~h^{-1})",
                            `ER` = "ER~(g~O^2~m^{-2}~h^{-1})",
                            `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})",
                            `PAR` = "PAR~({`mu`}*mol~m^{-2}~s^{-1})"))


