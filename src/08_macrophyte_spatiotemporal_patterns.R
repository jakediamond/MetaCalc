# -------------------------------------
# Author: Jake Diamond
# Purpose: look at macrophyte cover over time
# Date: 15 November 2022
# -------------------------------------
# load libraries
library(lubridate)
library(sf)
# library(plotly)
library(tidyverse)
library(units)
library(scales)
library(tmap)
library(tidytable)

# Load data ---------------------------------------------------------------
# Data directory
data_dir <- file.path("data", "01_EDF", "Shapefiles_surfaces")

# File paths
files <- fs::dir_ls(data_dir, regexp = "\\.shp$", recurse = TRUE)

# Load all shapefile data
df <- files %>% 
  map_dfr(st_read, .id = "source", quiet = TRUE, type = 6) %>%
  st_as_sf()

# Get date info
df <- df %>%
  dplyr::mutate(source = str_extract(source, "([^/]+$)"),
                source = str_remove(source, ".shp")) %>%
  dplyr::mutate(dateinfo = if_else(str_detect(source, "surface"),
                                   str_extract(source, "[^_]+"),
                                   str_extract(source, "(?<=_)[^_]+(?=$)"))) %>%
  dplyr::mutate(dateinfo = str_replace(dateinfo, "aout", "08"),
                dateinfo = str_replace(dateinfo, "juillet", "07"),
                dateinfo = str_replace(dateinfo, "region", "0817"),
                dateinfo = str_replace(dateinfo, "2015", "0715")) %>%
  dplyr::mutate(date = parse_date_time(dateinfo, c("%Y%m%d", "%m%y")))

# Get site names and macrophyte type
df_site <- df %>%
  dplyr::mutate(site = str_extract(source, "ampierre|hinon|aurent|ivaux|ville|ienne")) %>%
  dplyr::mutate(site = case_when.(site == "ampierre" ~ "dampierre",
                           site == "hinon" ~ "chinon",
                           site == "aurent" ~ "stlaurent",
                           site %in% c("ivaux", "ienne") ~ "civaux",
                           site == "ville" ~ "belleville")) %>%
  dplyr::mutate(vegtype = if_else(is.na(Nom_espece), Nom_Espece, Nom_espece),
                vegtype = case_when.(str_detect(vegtype, "fini") ~ "undefined",
                                     is.na(vegtype) & str_detect(source, "fini") ~ "undefined",
                                     is.na(vegtype) & str_detect(source, "otamo") ~ "Potamogeton nodosus",
                                     is.na(vegtype) & str_detect(source, "enon") ~ "Ranunculus penicillatus",
                                     is.na(vegtype) & str_detect(source, "ussie") ~ "Ludwigia grandiflora",
                                     !is.na(vegtype) ~ vegtype,
                                     TRUE ~ "undefined"))

# Get areas of all polygons and the sites
df_area <- df_site %>%
  select(site, vegtype, date) %>%
  dplyr::mutate(area = st_area(df)) %>%
  group_by(site, vegtype, date) %>%
  summarize(sumveg = sum(area))

# Get rid of units and refactor 
df_area <- drop_units(df_area) %>%
  mutate(site_f = factor(site),
         year = year(date))
df_area$site_f = fct_relevel(df_area$site_f, "belleville", "dampierre", "chinon", 
                     "civaux", "stlaurent")
# All data
df_sumarea <- ungroup(df_area) %>%
  mutate(year = year(date)) %>%
  group_by(site_f, year) %>%
  summarize(sumveg = sum(sumveg))

# plot over time
ggplot(data = df_area,
       aes(x = year, y = sumveg)) +
  stat_summary(aes(color = vegtype)) + 
  stat_summary(aes(color = vegtype), geom = "line") +
  facet_grid(cols = vars(site_f)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # scale_x_date(date_labels = "%y") +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#009E73", "black")) +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  labs(x = "",
       y = expression("macrophyte area ("*m^2*")"),
       title = "surfaces d’herbiers")
ggsave(file.path("results", "macrophytes", "timeseries_shapefiles.png"),
       dpi = 600,
       width = 20,
       height = 16,
       units = "cm")

# Plot spatial macrophytes over time
ggplot() + 
  geom_sf(data = df_site) + 
  facet_wrap(year(date)~site, ncol = 5) +
  theme_classic()

df_tmap <- df_site %>%
  dplyr::mutate(year = year(date),
                site_f = factor(site))
df_tmap$site_f = fct_relevel(df_tmap$site_f, "belleville", "dampierre", "chinon", 
                             "civaux", "stlaurent")
tmap_mode("plot")
tm_shape(df_tmap) +
  tm_polygons(col = "green") +
  # tm_fill(col = "vegtype", title = "Species") +
  tm_facets(by=c("site_f", "year")) +
  tm_compass() +
  tm_scale_bar()

x = dplyr::filter(df_tmap, site == "dampierre", year == 2016,
                   surface > 0)
y = st_centroid(x)

tm_shape(dplyr::filter(df_tmap, !(site == "dampierre" & year == 2016 & surface == 0)) %>%
           st_make_valid(.)) +
  tm_polygons(col = "green") +
  tm_fill(col = "vegtype", title = "Species") +
  tm_facets(by=c("site_f", "year")) +
    tm_layout(panel.label.size = 3.8)
tmap_save(
  filename = file.path("results", "macrophytes", "facets_shapefiles.png"),
  # device = "png",
  width = 20,
  height = 16,
  units = "cm",
  dpi = 600
)


# df_tmap2 <- df_tmap %>%
#   dplyr::group_by(site, year, vegtype) %>%
#   st_cast(st_union(st_make_valid(.$geometry)), "POLYGON")

df_tmap2 <- terra::aggregate(st_make_valid(dplyr::filter(df_tmap, surface > 0)), 
                             by = "site_f", 
                             dissolve = TRUE)

tm_shape(filter(df_tmap2, site == "chinon")) +
  # tm_polygons(col = "green") +
  tm_fill(col = "vegtype", title = "Species") +
  tm_facets(by="year") +
  tm_compass() +
  tm_scale_bar()

str(df_tmap2)
