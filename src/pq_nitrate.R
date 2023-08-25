# Load water chemistry data
df_wq <- readxl::read_xlsx(file.path("data", "00_water chemistry", 
                                     "WQ_Loire_Moyenne.xlsx"))

# Quick clean
df_wq <- dplyr::mutate(df_wq,
                       dplyr::across(everything(), ~ifelse(. == -1, NA_real_, .))) %>%
  filter(site_no %in% c("4046000",
                        "4046800",
                        "4048000")) %>% 
  mutate(date = ymd(paste(year, month, day)))
colnames(df_wq)

df_n <- df_wq %>%
  select(site_no, year, month, date, NH4, NO2, NO3) %>%
  group_by(site_no, year, month, date) %>%
  mutate(NH4 = NH4 / 18 *1000,
         NO2 = NO2 / 46 * 1000,
         NO3 = NO3 / 62 * 1000) |>
  mutate(DIN = rowSums(across(where(is.numeric)), na.rm = T)) %>%
  mutate(nit_rat = NO3 / DIN,
         CN = if_else(year < 2005, 8, 20)) %>%
  mutate(PQ_n = nit_rat * 2 / CN)

ggplot(data = df_n,
       aes(x = date,
           y = nit_rat)) +
  geom_point()

mean(df_n$nit_rat, na.rm =T)

c = df_n |>
  ungroup() |>
  group_by(year) |>
  summarize(rat = mean(nit_rat, na.rm = T))
