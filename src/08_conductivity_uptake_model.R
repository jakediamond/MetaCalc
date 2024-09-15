library(deSolve)
library(patchwork)
library(scales)
library(tidyverse)

source(file.path("src", "000_functions_for_model.R"))

# Equlibrium model --------------------------------------------------------
eqmodel <- function(time, states, pars){
  with(as.list(c(states, pars)), {
    # Constants
    TK    <- temp + 273.15 # temperature in K
    I     <- 1.6 * 10^-5 * cond # ionic strength, cond in uS/cm
    S     <- 53.974 * I #salinity from ionic strength, estimated
    gamma <- act(TK, cond) # -log(monovalent act. coef)
    aHCO3 <- 10^-gamma # activity coefficient for HCO3-
    aCO3  <- 10^(4 * -gamma) # activity coefficient for CO32-
    aCa   <- 10^(4 * -gamma) # activity coefficient for CO32-
    # eq conductivities at 25C (S cm2/mol)
    hco3_eq <- 44.5
    co3_eq  <- 69.3
    ca_eq   <- 59.47 #eq conductivity at 25C (S cm2/mol)
    
    #  Get the day of the year
    doy <- 180 + floor(time/24)
    # mean PAR at the time step (umol/m2/h)
    par <- par_fun(time, day = doy, latitude = 43)
    # mean daily PAR
    meanpar <- mean(par_fun(0:23, day = doy, latitude = 43))
    
    # Gross primary production (mol CO2/kg/h)
    gpp <- (gpp_mean / 24 / 32) / 1000 * (par / meanpar)
    
    # Ecosystem respiration (mol CO2/kg/h)
    er <- (er_mean / 24 / 32) / 1000
    
    # Carbonate equilibrium
    carb_eq <- carb(TK, ALK, cond = cond, TC = DIC)
    
    # Get the equilibrium values of the carbonate system
    pH    <- carb_eq$pH
    CO2   <- carb_eq$CO2
    HCO3  <- carb_eq$HCO3
    CO3   <- carb_eq$CO3
    
    # Calcite Saturation index
    SI <- - log10((aCa * Ca) * (aCO3 * CO3) / Ksp(TK, sal = S))

    # Change in DIC due to metabolism
    # Change in DIC is 2x when due to calcite precipitation
    dDIC <- (-gpp + er) * (1-alpha) + 2 * (-gpp + er) * (alpha)
    dALK <- 2*(-gpp + er) * (alpha)
    dCa <- 0.5 * dALK

    # Other outputs 
    NEP <- gpp - er
    hco3_sc <- HCO3 * 1000 * aHCO3 * hco3_eq
    co3_sc <- CO3 * 1000 * aCO3 * co3_eq * 2
    ca_sc   <- Ca * 1000 * aCa * ca_eq * 2
    SC  <- hco3_sc + co3_sc + ca_sc
    
    # Differences
    diffs <- c(dDIC = dDIC, dALK = dALK, dCa = dCa)
    return(list(diffs, CO2 = CO2, HCO3 = HCO3, CO3 = CO3, pH = pH, NEP = NEP, SC = SC, SI = SI))
  })
}

# Get simulation times ----------------------------------------------------
del_t    <- 1            # time step (s)
days    <- 24*3       # number of days to simulate
times    <- seq(0, days, # sequence of in-model times to simulate
                by = del_t)  

# Initial conditions ------------------------------------------------------
# Start at equilibrium with atmosphere
carbeq <- seacarb::carbfull(flag = 24,
                            var1 = 400, #pCO2 [uatm]
                            var2 = 2E-3, #alkalinity
                            S = 1.6 * 10^-5 * 250 * 53.974,
                            T = 25,
                            pHscale = "F")

# Initial conditions vector
yini_eq <- c(DIC  = as.numeric(carbeq$CO2 + carbeq$CO3 + carbeq$HCO3),
             Ca = 0.5 * carbeq$HCO3,
             ALK = 2E-3)

# Get parameters ----------------------------------------------------------
pars0 <- c(
  temp = 25,
  cond = 250,
  p_atm = 1,
  CO2_atm = 400,
  z = 1,
  gpp_mean = 25,
  er_mean = 25,
  K600 = 4,
  alpha = 0 # % of DIC uptake attributed to CaCO3
)

# Run models --------------------------------------------------------------
out_eq_0 <- ode(y = yini_eq,
              times = times,
              func = eqmodel,
              parms = pars0)

out_eq_25 <- ode(y = yini_eq,
                times = times,
                func = eqmodel,
                parms = replace(pars0, "alpha", 0.25))

out_eq_50 <- ode(y = yini_eq,
                 times = times,
                 func = eqmodel,
                 parms = replace(pars0, "alpha", 0.5))

out_eq_100 <- ode(y = yini_eq,
                 times = times,
                 func = eqmodel,
                 parms = replace(pars0, "alpha", 1))
# Get into dataframe and plot--------------------------------------------------
df_fun <- function(desolve) {
  as_tibble(desolve) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(across(-c(pH, time, SC, SI), ~. * 1E6))
}

df <- list("0" = out_eq_0, "25" = out_eq_25, "50" = out_eq_50, "100" = out_eq_100) %>%
  map(., df_fun) %>%
  map_dfr(., bind_rows, .id= "calc") %>%
  mutate(calc = factor(calc, levels = c("0", "25", "50", "100")))

p_ts<- df %>% 
  group_by(calc) %>%
  slice_tail(n = 24) %>%
  ggplot(aes(x = time %% 24,
                  y = SC,
                  color = calc)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(name = expression("Calc.,"~alpha~"(%)"), 
                     values = c('#009E73', '#56B4E9', '#D55E00', "#4d4d4d")) +
  theme_classic(base_size = 10) +
  annotate(geom = "rect", xmin = 6, xmax = 18, ymin = 133, ymax = 130, fill = "gold") +
  annotate(geom = "text", x = 1, y = 207, label = "c") +
  scale_x_continuous(breaks = seq(0,24, 6)) +
  coord_cartesian(expand = FALSE) +
  theme(legend.position = c(0.82, 0.83),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = "transparent", color = "black")) +
  labs(x = "time (h)",
       y = expression(C[25]~"("*mu*S~cm^{-1}*")"))
p_ts

p_nep <- df %>% 
  group_by(calc) %>%
  slice_tail(n = 24) %>%
  mutate(delSC = (SC- lag(SC))) %>%
  ggplot(aes(x = NEP,
             y = delSC,
             color = calc))+
  geom_path(size = 1.5) +
  scale_color_manual(name = "Calc. %", values = c('#009E73', '#56B4E9', '#D55E00', "#4d4d4d")) +
  theme_classic(base_size = 10) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  ggpubr::stat_regline_equation(show.legend = FALSE, label.y.npc = 0.38, size = 3.2) +
  theme(legend.position = "none") +
  annotate(geom = "text", x = -35, y = 7.5, label = "d") +
  labs(x = expression(NEP~"("*mmol~m^{-3}~h^{-1}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-1}~h^{-1}*")"))
p_nep


p_ts_all <- ggplot(pivot_longer(df, -c(time, calc)), 
              aes(x = time / 24,
                  y = value,
                  color = calc)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_color_manual(name = "Calc. %", values = c('#009E73', '#56B4E9', '#D55E00', "#4d4d4d")) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.9, 0.15)) +
  facet_wrap(~name, scales = "free_y") + 
  labs(x = "time (d)",
       y = expression(mu*M~or~mu*M~s^{-1}~or~mu*S~cm^{-1}))
p_ts_all
ggsave(plot = p_ts_all, filename = file.path("results", "supplementary_cond_model.png"),
       dpi = 300,
       units = "cm",
       height = 12,
       width = 18.4)



# uptake conductivity effects ---------------------------------------------
# Model the carbonate system of the Loire under varying CO2 conditions
# Make a log sequence of increasing CO2 concentrations
co2 <-  1.5e-07 * 2^(0:log((800e-06/1.5e-07), 2))
dic <- seq(1.15, 2.75, 0.02)
carb(298, 1.7, cond = 250, TC = 1.15)
carb(298, 1.7, cond = 250, TC = 2.75)
carb(298, 1.7, cond = 250, pH = 10)
# activity coeff
ac_hco3 <- 0.93
ac_co3  <- 0.93^4
ac_no3  <- 0.93
# eq conductivities at 25C (S cm2/mol)
hco3_eq <- 44.5
co3_eq  <- 69.3
ca_eq   <- 59.47
no3_eq  <- 71.42

# function to get seacarb based on co2
carbfun <- function (x) {
  seacarb::carb(flag = 4, 
                var1 = x,
                var2 = 0.0017,
                T = 25,
                S = 1.5E-5 * 350 * 60,
                k1k2 = "m06",
                pHscale = "F")
}

# Get a dataframe of the results
df_sc <- #tibble(DIC = dic) %>%
  # rowwise() %>%
  # mutate(carbs = carb(298, 1.7, cond = 250, TC = DIC)) %>%
  # ungroup() %>%
tibble(CO2 = co2) |>
# Get the carbonate system based on CO2 and ALK, mol/m3 (outputs in mol/kg)
mutate(carbs = carbfun(CO2)) %>%
mutate(hco3 = carbs$HCO3 * 1000,
       co3 = carbs$CO3 * 1000,
       pH = carbs$pH,
       co2 = carbs$CO2 * 1000) %>%
  select(-carbs) %>%
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         SC = hco3_sc + co3_sc ) |> 
  # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
         dSC = SC - lag(SC))

# Plot of the carbonate system
p1 <- ggplot(data = df_sc,
             aes(x = co2 * 1000)) +
  geom_line(aes(y = hco3 * 1000), linewidth = 1.2) +
  annotate(geom = "text", x = 0.01 * 1000, y = 0.87 * 1000, label = "HCO[3]^{`-`}", parse = TRUE) +
  annotate(geom = "text", x = 0.0055 * 1000, y = 3E-2 * 1000, label = "CO[3]^{`2-`}", parse = TRUE) +
  geom_line(aes(y = co3 * 1000), linetype = "dashed", linewidth = 1.2) + 
  scale_y_log10() +
  scale_x_continuous(trans = c("log10", "reverse")) +
  annotation_logticks(sides = "bl") +
  geom_segment(aes(x = 0.038 * 1000, y = 0.002 * 1000, xend = 0.0005 * 1000, yend = 0.002 * 1000),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x = 0.0055 * 1000, y = 0.0033 * 1000, label = "decreasing~CO[2]", parse = TRUE) +
  annotate(geom = "text", x = 600, y = 2300, label = "a") +
  theme_classic(base_size = 10) +
  labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
       y = expression("("*mmol~m^{-3}*")"))
p1


# # Plot of the carbonate system
# p1_dic <- ggplot(data = df_sc,
#              aes(x = DIC * 1000)) +
#   geom_line(aes(y = hco3 * 1000), linewidth = 1.2) +
#   # annotate(geom = "text", x = 0.01 * 1000, y = 0.9 * 1000, label = "HCO[3]^{`-`}", parse = TRUE) +
#   # annotate(geom = "text", x = 0.0055 * 1000, y = 3E-2 * 1000, label = "CO[3]^{`2-`}", parse = TRUE) +
#   geom_line(aes(y = co3 * 1000), linetype = "dashed", linewidth = 1.2) + 
#   geom_line(aes(y = co2 * 1000), linetype = "dotted", linewidth = 1.2) + 
#   scale_y_log10() +
#   scale_x_continuous(trans = "reverse") +
#   annotation_logticks(sides = "bl") +
#   # geom_segment(aes(x = 0.035 * 1000, y = 0.002 * 1000, xend = 0.0005 * 1000, yend = 0.002 * 1000),
#   #              arrow = arrow(length = unit(0.5, "cm"))) +
#   # annotate(geom = "text", x = 0.0054 * 1000, y = 0.0035 * 1000, label = "decreasing~CO[2]", parse = TRUE) +
#   # annotate(geom = "text", x = 600, y = 2200, label = "a") +
#   theme_classic(base_size = 10) +
#   labs(x = expression(DIC~"("*mmol~m^{-3}*")"),
#        y = expression("("*mmol~m^{-3}*")"))
# p1_dic

# Plot of pH
ppH <- ggplot(data = df_sc,
              aes(x = co2)) +
  geom_line(aes(y = pH), linewidth = 1.2, color = "red") +
  scale_y_continuous(position = "right") +
  scale_x_continuous(trans = c("log10", "reverse")) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(color = "red"),
        axis.title.y = element_text(color = "red"),
        axis.ticks.y = element_line(color= "red"),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) +
  labs(x = "",
       y = "pH")
ppH


# # Plot of the local slope of SpC based on CO2
# p2 <- ggplot(data = df_sc,
#              aes(x = co2 * 1000)) +
#   geom_line(aes(y = -dhco3_sc), linewidth = 1.2) +
#   annotate(geom = "text", x = 0.001 * 1000, y = -6.3, label = "Delta*C[HCO[3]^{`-`}]", parse = TRUE) +
#   annotate(geom = "text", x = 0.003 * 1000, y = 5, label = "Delta*C[CO[3]^{`2-`}]", parse = TRUE) +
#   annotate(geom = "text", x = 0.002 * 1000, y = -0.2, label = "Delta*C[25]", parse = TRUE) +
#   geom_line(aes(y = -dco3_sc), linetype = "dashed", linewidth = 1.2) +
#   geom_line(aes(y = -dSC), linetype = "dotted", linewidth = 1.2) + 
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(trans = c("log10", "reverse")) +
#   annotation_logticks(sides = "b") +
#   geom_segment(aes(x = 0.1 * 1000, y = -7.5, xend = 0.003 * 1000, yend = -7.5),
#                arrow = arrow(length = unit(0.5, "cm"))) +
#   annotate(geom = "text", x = 0.02 * 1000, y = -6.8, label = "decreasing~CO[2]", parse = TRUE) +
#   annotate(geom = "text", x = 95, y = 9.5, label = "b") +
#   theme_classic() +
#   labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
#        y = expression(Delta*C[25]~"("*mu*S~cm^{-2}*")"))
# p2

# Plot of just the conductivity
p3 <- ggplot(data = df_sc,
             aes(x = co2 * 1000)) +
  geom_line(aes(y = hco3_sc), linewidth = 1.2) +
  annotate(geom = "text", x = 1, y = 60, label = "C[HCO[3]^{`-`}]", parse = TRUE) +
  annotate(geom = "text", x = 0.8, y = 15, label = "C[CO[3]^{`2-`}]", parse = TRUE) +
  annotate(geom = "text", x = 10, y = 78, label = "C[HCO[3]^{`-`}]+C[CO[3]^{`2-`}]", parse = TRUE) +
  geom_line(aes(y = co3_sc), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = SC), linetype = "dotted", linewidth = 1.2) + 
  scale_x_continuous(trans = c("log10", "reverse")) +
  annotation_logticks(sides = "b") +
  geom_segment(aes(x = 620, y = 10, xend = 10, yend = 10),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x = 110, y = 14.5, label = "decreasing~CO[2]", parse = TRUE) +
  annotate(geom = "text", x = 600, y = 80, label = "b") +
  theme_classic(base_size = 10) +
  labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
       y = expression(contribution~to~C[25]~"("*mu*S~cm^{-2}*")"))
p3



# Overall plots -----------------------------------------------------------
# Combined plot
p1pH <- (p1 + ppH) + plot_layout(design = c(area(t = 1, l = 1, b = 5, r = 5),
                                            area(t = 1, l = 1, b = 5, r = 5)))
p1pH

# Overall combined plot
p <- (p1pH | p3) / (p_ts | p_nep)
p

ggsave(plot = p,
       filename = file.path("results", "figure1_cond_changes.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12.2)
