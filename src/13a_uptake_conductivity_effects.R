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
# Low CO2 for pH 9.7 and ALK 1.7 = 3.288592e-07
# Make a log sequence of increasing CO2 concentrations
co2 <-  3.288592e-07 * 2^(0:log((9.989661e-05/3.288592e-07), 2))

# activity coeff
ac_hco3 <- 0.93
ac_co3 <- 0.93^4

# eq conductivities at 25C (S cm2/mol)
hco3_eq <- 44.5
co3_eq <- 69.3
ca_eq <- 59.47

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
df <- tibble(co2 = co2) |>
  # Get the carbonate system based on CO2 and ALK, mol/m3 (outputs in mol/kg)
  mutate(hco3 = carbfun(co2)$HCO3 * 1000,
         co3 = carbfun(co2)$CO3 * 1000,
         pH = carbfun(co2)$pH,
         co2 = co2 * 1000) |>
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         SC = hco3_sc + co3_sc) |>
  # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
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
# Overall combined plot
p <- p1pH | p2
p

ggsave(plot = p,
       filename = file.path("results", "figureS_cond_changes_co2_uptake_mmol.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 9.2)

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


