# -------------------------------------
# Author: Jake Diamond
# Purpose: DIC removal effects on specific conductance
# Date: 2023-04-28
# -------------------------------------

# Load libraries
library(tidyverse)
library(scales)
library(patchwork)

# Model the carbonate system of the Loire under varying CO2 conditions
# High CO2 for pH 7.5 and Alk 1.7 = 9.989661e-05
# Low CO2 for pH 10 and ALK 1.7 = 3.288592e-07
# Make a log sequence of increasing CO2 concentrations
# co2 <-  3.288592e-07 * 2^(0:log((9.989661e-05/3.288592e-07), 2))
co2 <-  1.5e-07 * 2^(0:log((800e-06/1.5e-07), 2))
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

# assume NEP = 10 mmol CO2 /m3/hr = 0.01 mol CO2/m3/hr = 0.00001 mol CO2/L/hr
co2_nep <- seq(5E-6, 4E-6, -1E-7)
no3_nep <- seq(3.7/62, 3.7/62 - 1E-4, -1E-4/10)

df_nep <- tibble(co2 = co2_nep,
                 no3 = no3_nep) |>
  # Get the carbonate system based on CO2 and ALK, mol/m3 (outputs in mol/kg)
  mutate(hco3 = carbfun(co2)$HCO3 * 1000,
         co3 = carbfun(co2)$CO3 * 1000,
         pH = carbfun(co2)$pH,
         co2 = co2 * 1000) %>%
  # arrange(desc(co2)) %>%
  # mutate(no3 = 3.7/62 - (co2 - lag(co2)) / 10) %>% #C:N = 10
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         no3_sc = no3 * ac_no3 * no3_eq,
         SC = hco3_sc + co3_sc ) |> #+ no3_sc
  # # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
         dno3_sc = no3_sc - lag(no3_sc),
         dSC = SC - lag(SC))


# Get a dataframe of the results
df <- tibble(co2 = co2) |>
  # Get the carbonate system based on CO2 and ALK, mol/m3 (outputs in mol/kg)
  mutate(hco3 = carbfun(co2)$HCO3 * 1000,
         co3 = carbfun(co2)$CO3 * 1000,
         pH = carbfun(co2)$pH,
         co2 = co2 * 1000) %>%
  # arrange(desc(co2)) %>%
  # mutate(no3 = 3.7/62 - (co2 - lag(co2)) / 10) %>% #C:N = 10
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         # no3_sc = no3 * ac_no3 * no3_eq,
         SC = hco3_sc + co3_sc ) |> #+ no3_sc
  # # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
         # dno3_sc = no3_sc - lag(no3_sc),
         dSC = SC - lag(SC))

# Plot of the carbonate system
p1 <- ggplot(data = df,
       aes(x = co2 * 1000)) +
  geom_line(aes(y = hco3 * 1000), linewidth = 1.2) +
  annotate(geom = "text", x = 0.01 * 1000, y = 1.1 * 1000, label = "HCO[3]^{`-`}", parse = TRUE) +
  annotate(geom = "text", x = 0.008 * 1000, y = 3E-2 * 1000, label = "CO[3]^{`2-`}", parse = TRUE) +
  geom_line(aes(y = co3 * 1000), linetype = "dashed", linewidth = 1.2) + 
  scale_y_log10() +
  scale_x_continuous(trans = c("log10", "reverse")) +
  annotation_logticks(sides = "bl") +
  geom_segment(aes(x = 0.03 * 1000, y = 0.007 * 1000, xend = 0.0009 * 1000, yend = 0.007 * 1000),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x = 0.005 * 1000, y = 0.009 * 1000, label = "decreasing~CO[2]", parse = TRUE) +
  annotate(geom = "text", x = 120, y = 1900, label = "a") +
  theme_classic() +
  labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
       y = expression("("*mmol~m^{-3}*")"))
p1

# Plot of pH
ppH <- ggplot(data = df,
             aes(x = co2)) +
  geom_line(aes(y = pH), linewidth = 1.2, color = "red") +
  scale_y_continuous(position = "right") +
  scale_x_continuous(trans = c("log10", "reverse")) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "red"),
        axis.title.y = element_text(color = "red"),
        axis.ticks.y = element_line(color= "red"),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) +
  labs(x = "",
       y = "pH")
ppH

# Combined plot
p1pH <- (p1 + ppH) + plot_layout(design = c(area(t = 1, l = 1, b = 5, r = 5),
                                          area(t = 1, l = 1, b = 5, r = 5)))
p1pH

# Plot of the local slope of SpC based on CO2
p2 <- ggplot(data = df,
       aes(x = co2 * 1000)) +
  geom_line(aes(y = -dhco3_sc), linewidth = 1.2) +
  annotate(geom = "text", x = 0.001 * 1000, y = -6.3, label = "Delta*C[HCO[3]^{`-`}]", parse = TRUE) +
  annotate(geom = "text", x = 0.003 * 1000, y = 5, label = "Delta*C[CO[3]^{`2-`}]", parse = TRUE) +
  annotate(geom = "text", x = 0.002 * 1000, y = -0.2, label = "Delta*C[25]", parse = TRUE) +
  geom_line(aes(y = -dco3_sc), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = -dSC), linetype = "dotted", linewidth = 1.2) + 
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = c("log10", "reverse")) +
  annotation_logticks(sides = "b") +
  geom_segment(aes(x = 0.1 * 1000, y = -7.5, xend = 0.003 * 1000, yend = -7.5),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x = 0.02 * 1000, y = -6.8, label = "decreasing~CO[2]", parse = TRUE) +
  annotate(geom = "text", x = 95, y = 9.5, label = "b") +
  theme_classic() +
  labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-2}*")"))
p2

# Plot of just the conductivity
p3 <- ggplot(data = df,
             aes(x = co2 * 1000)) +
  geom_line(aes(y = hco3_sc), linewidth = 1.2) +
  annotate(geom = "text", x = 1, y = 60, label = "C[HCO[3]^{`-`}]", parse = TRUE) +
  annotate(geom = "text", x = 3, y = 5, label = "C[CO[3]^{`2-`}]", parse = TRUE) +
  annotate(geom = "text", x = 10, y = 76, label = "C[HCO[3]^{`-`}]+C[CO[3]^{`2-`}]", parse = TRUE) +
  geom_line(aes(y = co3_sc), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = SC), linetype = "dotted", linewidth = 1.2) + 
  scale_x_continuous(trans = c("log10", "reverse")) +
  annotation_logticks(sides = "b") +
  geom_segment(aes(x = 600, y = 10, xend = 10, yend = 10),
  arrow = arrow(length = unit(0.5, "cm"))) +
  annotate(geom = "text", x = 70, y = 14, label = "decreasing~CO[2]", parse = TRUE) +
  annotate(geom = "text", x = 600, y = 80, label = "b") +
  theme_classic() +
  labs(x = expression(CO[2]~"("*mmol~m^{-3}*")"),
       y = expression(contribution~to~C[25]~"("*mu*S~cm^{-2}*")"))
p3

# Overall combined plot
p <- (p1pH | p3) / p2
p

ggsave(plot = p,
       filename = file.path("results", "figureS_cond_changes_co2_uptake_mmol_v2.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12.2)

pbeta <- ggplot(data = df,
             aes(x = co2)) +
  # geom_line(aes(y = dhco3_sc), linewidth = 1.2) +
  # annotate(geom = "text", x = 0.004, y = 5, label = "HCO[3]^{`-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.004, y = -6, label = "CO[3]^{`2-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.002, y = -1, label = "C[25]", parse = TRUE) +
  # geom_line(aes(y = dco3_sc), linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = -dSC/(dco2*1000)), linetype = "dotted", linewidth = 1.2) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = expression(CO[2]~"("*mol~m^{-3}*")"),
       y = expression(beta[C[25]]~"("*mu*S~cm^{-2}*")"))
pbeta

# Look at buffer effects --------------------------------------------------

seacarb::carb(flag=15, var1=2.4E-6, var2=2.13E-6, S=35, T=15, Patm=1, 
              pHscale="F")

carbeq1 <- seacarb::carb(flag=8, var1=8.2, var2=0.00170, S=0.2, T=25, Patm=1, 
                         pHscale="F") %>%
  unlist()

names(carbeq1) <- paste0(names(carbeq1), "test")
carbeq1
buf1 <- seacarb::buffer(flag=8, var1=8.2, var2=0.00170, S=0.2, T=25, Patm=1, 
                        pHscale="F")
buf1
?seacarb::buffer
phib <- buf1$PhiB
# subtract 1E-5 mol/kg HCO3 to see the change in pH
delpH <- phib * 1E-5

# recalculate carbeq
carbeq2 <- seacarb::carb(flag=8, var1=8.2 + delpH, var2=0.00170 - 1E-5, S=0.2, T=25, Patm=1, 
                         pHscale="F")
carbeq1
carbeq2


# function to get seacarb based on co2
carbfunbuf <- function (x) {
  carbeq1 <- seacarb::carb(flag=4, 
                           var1=x, 
                           var2=2.4E-6,#0.00170, 
                           S=35,#1.5E-5 * 350 * 60, 
                           T=15,#25, 
                           pHscale="F")
  buf1 <- seacarb::buffer(flag=8, 
                          var1=8.2, 
                          var2=2.4E-6,#0.00170, 
                          S=35,#1.5E-5 * 350 * 60, 
                          T=15,#25, 
                          pHscale="F")
  phib <- buf1$PhiB
  phid <- buf1$PhiD
  # subtract 10x1E-6 mol/kg HCO3 (b) and CO2 (d) to see the change in pH
  dpH_b <- phib * -1E-5
  dpH_d <- phid * -1E-5
  
  carbeq2_b <- seacarb::carb(flag = 8, 
                var1 = carbeq1$pH + dpH_b,
                var2=2.4E-6,#0.00170, 
                S=35,#1.5E-5 * 350 * 60, 
                T=15,#25, 
                k1k2 = "m06",
                pHscale = "F") %>%
    unlist()
  # names(carbeq2_b) <- paste0(names(carbeq2_b), "_HCO3")
  
  carbeq2_d <- seacarb::carb(flag = 8, 
                             var1 = carbeq1$pH + dpH_d,
                             var2=2.4E-6,#0.00170, 
                             S=35,#1.5E-5 * 350 * 60, 
                             T=15,#25, 
                             k1k2 = "m06",
                             pHscale = "F") %>%
    unlist()
  # names(carbeq2_d) <- paste0(names(carbeq2_d), "_CO2")
  result <- data.frame(type = c(rep("HCO3", length(carbeq2_b)),
                                rep("CO2", length(carbeq2_d))),
                       variable = c(names(carbeq2_b), names(carbeq2_d)),
                       value = c(as.numeric(carbeq2_b), as.numeric(carbeq2_d))
                       )
}

# Assuming a pulse of GPP equal to 10 mmol/m2/hr = 1E-6 mol/kg/hr with depth = 1
# Get a dataframe of the results
df_buf <- tibble(co2 = co2) |>
  # Get the carbonate system based on CO2 and ALK, mol/m3 (outputs in mol/kg)
  mutate(a = carbfunbuf(co2))


x <- carbfunbuf(6.831E-7)

# old ---------------------------------------------------------------------


