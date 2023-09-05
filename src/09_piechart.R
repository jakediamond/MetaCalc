# -------------------------------------
# Author: Jake Diamond
# Purpose: Get a pie chart of temporal occurence trophlux states
# Date: 3 March 2023
# -------------------------------------
# load data
df <- readRDS(file.path("data", "daily_trophlux.RDS"))

# get trophlux summary
df_sum <- df |>
  group_by(trophlux) |>
  summarize(n = n()) |>
  ungroup() |>
  mutate(tot = sum(n),
         per = n/tot*100) |>
  arrange(desc(per))
  

# Create Data
data <- data.frame(
  group=c("4", "3", "2", "1"),
  value=df_sum$per
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=-3.14/1.4) +
  theme_void() + 
  theme(legend.position="none") +
  # geom_text(aes(y = ypos, label = round(value)), color = "white", size=3.4) +
  scale_fill_manual(values = c("#0072B2", "black", "#E69F00",  "#009E73"))
ggsave(filename = file.path("results", "piechart.png"),
       dpi = 300,
       units = "cm",
       height = 2.54,
       width = 2.54)
