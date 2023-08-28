library(tidyverse)
library(scales)
library(patchwork)

?seacarb::carb
seacarb::carb(flag = 8,
              var1 = 9.7,
              var2 = 0.0017,
              T = 25,
              1.5E-5 * 350 * 60,
              k1k2 = "m06")
# High CO2 for pH 7.5 and Alk 1.7 = 9.989661e-05
# Low CO2 for pH 9.7 and ALK 1.7 = 3.288592e-07

co2 <- seq(6.37E-7, 9.989661e-05, length.out = 100)
co2 <-  3.288592e-07 * 2^(0:log((9.989661e-05/3.288592e-07), 2))
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

# activity coeff
ac_hco3 <- 0.93
ac_co3 <- 0.93^4

# eq conductivities at 25C (S cm2/mol)
hco3_eq <- 44.5
co3_eq <- 69.3
ca_eq <- 59.47

df <- tibble(co2 = co2) |>
  # Get the carbonate system based on CO2 and ALK
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

p1pH <- (p1 + ppH) + plot_layout(design = c(area(t = 1, l = 1, b = 5, r = 5),
                                          area(t = 1, l = 1, b = 5, r = 5)))
p1pH

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
  labs(x = expression(CO[2]~"("*mol~m^{-3}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-2}*")"))
p2
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


# Same for HCO3 uptake ----------------------------------------------------
?seacarb::carb
seacarb::carb(flag = 8,
              var1 = 9.7,
              var2 = 0.0017,
              T = 25,
              S = 1.5E-5 * 350 * 60,
              k1k2 = "m06")
# High HCO3 for pH 7.5 and Alk 1.7 = 0.001690458 
# Low HCO3 for pH 9.7 and ALK 1.7 = 0.0008819898

# High HCO3 for pH 7.5 and Alk 1.7 = 0.001690458 
# Low HCO3 for pH 9.7 and ALK 1.7 = 0.0008819898

hco3 <- seq(0.0008819898, 0.001690458 , length.out = 20)

dic

# function to get seacarb based on hco3
carbfun_hco3 <- function (x) {
  seacarb::carb(flag = 11, 
                var1 = x,
                var2 = 0.0017,
                T = 25,
                S = 1.5E-5 * 350 * 60,
                k1k2 = "m06",
                pHscale = "F")
}

df_hco3 <- tibble(hco3 = hco3) |>
  # Get the carbonate system based on CO2 and ALK
  mutate(co2 = carbfun_hco3(hco3)$CO2 * 1000,
         co3 = carbfun_hco3(hco3)$CO3 * 1000,
         pH = carbfun_hco3(hco3)$pH,
         dic = carbfun_hco3(hco3)$DIC * 1000,
         hco3 = hco3 * 1000) |>
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         SC = hco3_sc + co3_sc) |>
  # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
         dSC = SC - lag(SC))


p3 <- ggplot(data = df_hco3,
             aes(x = pH)) +
  geom_line(aes(y = hco3_sc), linewidth = 1.2) +
  # annotate(geom = "text", x = 8, y = 1.1, label = "HCO[3]^{`-`}", parse = TRUE) +
  # annotate(geom = "text", x = 8, y = 3E-2, label = "CO[3]^{`2-`}", parse = TRUE) +
  geom_line(aes(y = co3_sc), linetype = "dashed", linewidth = 1.2) +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic() +
  labs(x = expression(HCO[3]^{`-`}~"("*mol~m^{-3}*")"),
       y = expression("("*mol~m^{-3}*")"))
p3

p4 <- ggplot(data = df_hco3,
             aes(x = pH)) +
  geom_line(aes(y = SC), linewidth = 1.2) +
  # geom_line(aes(y = -hco3_sc), linewidth = 1.2) +
  # annotate(geom = "text", x = 0.009, y = 0.01, label = "HCO[3]^{`-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.009, y = -0.005, label = "CO[3]^{`2-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.005, y = 0.001, label = "C[25]", parse = TRUE) +
  # geom_line(aes(y = -co3_sc), linetype = "dashed", linewidth = 1.2) + 
  # geom_line(aes(y = dSC / (dhco3 * ac_hco3)), linetype = "dotted", linewidth = 1.2) +
  # scale_x_log10() +
  theme_classic() +
  labs(x = "pH", #expression(HCO[3]^{`-`}~"("*mol~m^{-3}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-2}*")"))
p4
p5 <- p3 | p4
p5
ggsave(plot = p,
       filename = file.path("results", "cond_changes_co2_uptake.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12)

# same for dic ------------------------------------------------------------
seacarb::carb(flag = 8,
              var1 = 9.7,
              var2 = 0.0017,
              T = 25,
              S = 1.5E-5 * 350 * 60,
              k1k2 = "m06")
# High DIC for pH 7.5 and Alk 1.7 = 0.001794901
# Low DIC for pH 9.7 and ALK 1.7 = 0.001258279

dic <- seq(0.001258279, 0.001794901 , length.out = 20)

# function to get seacarb based on hco3
carbfun_dic <- function (x) {
  seacarb::carb(flag = 15, 
                var1 = 0.0017,
                var2 = x,
                T = 25,
                S = 1.5E-5 * 350 * 60,
                k1k2 = "m06",
                pHscale = "F")
}

df_dic <- tibble(dic = dic) |>
  # Get the carbonate system based on CO2 and ALK
  mutate(co2 = carbfun_dic(dic)$CO2 * 1000,
         co3 = carbfun_dic(dic)$CO3 * 1000,
         pH = carbfun_dic(dic)$pH,
         hco3 = carbfun_dic(dic)$HCO3 * 1000,
         dic = dic * 1000) |>
  # Get the contribution to specific conductance of each
  mutate(hco3_sc = hco3 * ac_hco3 * hco3_eq,
         co3_sc = co3 * ac_co3 * co3_eq * 2,
         SC = hco3_sc + co3_sc) |>
  # Get the derivative of each to determine delta SC
  mutate(dco2 = co2 - lag(co2),
         dhco3_sc = hco3_sc - lag(hco3_sc),
         dco3_sc = co3_sc - lag(co3_sc),
         dSC = SC - lag(SC))


p3 <- ggplot(data = df_hco3,
             aes(x = pH)) +
  geom_line(aes(y = hco3_sc), linewidth = 1.2) +
  # annotate(geom = "text", x = 8, y = 1.1, label = "HCO[3]^{`-`}", parse = TRUE) +
  # annotate(geom = "text", x = 8, y = 3E-2, label = "CO[3]^{`2-`}", parse = TRUE) +
  geom_line(aes(y = co3_sc), linetype = "dashed", linewidth = 1.2) +
  scale_y_log10() +
  # scale_x_log10() +
  theme_classic() +
  labs(x = expression(HCO[3]^{`-`}~"("*mol~m^{-3}*")"),
       y = expression("("*mol~m^{-3}*")"))
p3

p4 <- ggplot(data = df_hco3,
             aes(x = pH)) +
  geom_line(aes(y = SC), linewidth = 1.2) +
  # geom_line(aes(y = -hco3_sc), linewidth = 1.2) +
  # annotate(geom = "text", x = 0.009, y = 0.01, label = "HCO[3]^{`-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.009, y = -0.005, label = "CO[3]^{`2-`}", parse = TRUE) +
  # annotate(geom = "text", x = 0.005, y = 0.001, label = "C[25]", parse = TRUE) +
  # geom_line(aes(y = -co3_sc), linetype = "dashed", linewidth = 1.2) + 
  # geom_line(aes(y = dSC / (dhco3 * ac_hco3)), linetype = "dotted", linewidth = 1.2) +
  # scale_x_log10() +
  theme_classic() +
  labs(x = "pH", #expression(HCO[3]^{`-`}~"("*mol~m^{-3}*")"),
       y = expression(Delta*C[25]~"("*mu*S~cm^{-2}*")"))
p4
p5 <- p3 | p4
ggsave(plot = p,
       filename = file.path("results", "cond_changes_co2_uptake.png"),
       dpi = 300,
       units = "cm",
       width = 18.4,
       height = 12)


fm_fun <- function(temp, I){
  TK = temp + 273.15 #temp in Kelvin
  # I = 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  E = 60954/(TK+116) - 68.937 #dielectric constant
  # Davies estimate
  pfm = 1.82 * 10^6 * (E*TK)^-1.5 * ((sqrt(I)/(1+sqrt(I))) - 0.3 * I)
  fm = 10^-pfm
}

