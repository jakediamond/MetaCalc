# -------------------------------------
# Author: Jake Diamond
# Purpose: Plot the chemical enhancement factor
# Date: 6 March 2023
# -------------------------------------

# Load libraries
library(tidyverse)

df <- readRDS(file.path("data", "hourly_data_final.RDS"))

ggplot(data = drop_na(df, trophlux),
       aes(x = enh)) +
  stat_density(aes(y = after_stat(scaled))) +
  theme_classic() +
  facet_grid(troph~sourcesink) + 
  scale_x_continuous(limits = c(1, 2)) +
  labs(x = expression(alpha[enh]~"("*`-`*")"),
       y = "scaled density")

ggsave(file.path("results", "alpha_enh_scaled_density.png"),
       dpi = 300,
       units = "cm",
       width = 9.2,
       height = 9.2)

df |>
  filter(sourcesink == "sink") |>
  summarize(q = quantile(enh, na.rm = T),
            m = mean(enh, na.rm = T),
            sd = sd(enh, na.rm = T))
