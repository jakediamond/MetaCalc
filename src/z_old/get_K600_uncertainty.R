library(streamMetabolizer)
library(tidyverse)
df <- readRDS("Data/Loire_DO/metab_extremelyconstrainedK_gppconstrained_all_discharge_bins_fullmodel")
mm <- df %>%
  mutate(met = map(mm, get_params)) %>%
  select(-mm) %>%
  unnest(met) %>%
  select(date, K600.daily, K600.daily.sd)

mm2 <- df %>%
  mutate(met = map(mm, get_fit)) %>%
  select(-mm)

x = unnest(mm2, met)

y=pluck(x, 2, 1) %>%
  bind_rows(pluck(x, 2, 8)) %>%
  bind_rows(pluck(x, 2, 15)) %>%
  bind_rows(pluck(x, 2, 22)) %>%
  bind_rows(pluck(x, 2, 29)) %>%
  select(date, starts_with("K600")) %>%
  select(-contains("pred"))

write_csv2(y, "dampierre_K600_estimates.csv")
