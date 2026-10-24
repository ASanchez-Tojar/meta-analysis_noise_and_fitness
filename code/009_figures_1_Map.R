################################################################################
# Authors: 
#
# Jorge Pereira (University of Coimbra): original code
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): code ajusting

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the models run in scripts 005 to 008, and plot some
# figures for the main text of the meta-analysis:

# Noise as a global change driver: a meta-analysis of its fitness consequences 
# for birds

# by Alfredo Sánchez-Tójar, Nicholas P. Moran, Lena de Framond, Alberto Comin, 
# Henrik Brumm

################################################################################
# Packages needed
################################################################################

# install.packages("pacman")
pacman::p_load(readr,
               dplyr,
               ggplot2,
               ggplot2,
               patchwork,
               sf,
               rnaturalearth)

# cleaning environment
rm(list=ls())

################################################################################
# MAPPING THE EVIDENCE
################################################################################

# Study_ID is a code, not a number, so we read it as text
#es_raw <- read_csv("map_dataset_noise_and_birds.csv", col_types = cols(Study_ID = col_character()))
# es_raw <- read_csv("map_dataset_noise_and_birds.csv", col_types = cols(Study_ID = col_character()))

# Load cleaned dataset
es_raw <- read.csv("data/04_processed/noisemeta_datasheet_processed_3.csv",
                   header=T)

glimpse(es_raw)
summary(es_raw$yi.SMD.H.signed)

# The original names are long, so we give shorter ones
es <- es_raw |>
  rename(yi      = yi.SMD.H.signed,
         lat     = Latitude_N,
         lon     = Longitude_E,
         species = species.updated)

# How much evidence do we have?
n_effects <- nrow(es)
n_studies <- n_distinct(es$Study_ID)
n_species <- n_distinct(es$species)
n_sites   <- nrow(distinct(es, lat, lon))

n_effects
n_studies
n_species
n_sites

################################################################################
# Turn the table into a spatial object.

# coords = which columns are x and y (longitude first!)
# crs = 4326 means plain longitude/latitude (WGS84, the GPS system)
# remove = FALSE keeps lon and lat as normal columns too
es_sf <- st_as_sf(es, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Country polygons from Natural Earth (scale 50 = medium detail)
world <- ne_countries(scale = 50, returnclass = "sf")

# Many effect sizes share the same coordinates, so points are on top of each
# other. We make a table with one row per site (unique lat/lon)
sites <- es |>
  group_by(lat, lon) |>
  summarise(n_es      = n(),                    # effect sizes at the site
            n_studies = n_distinct(Study_ID),   # studies at the site
            mean_yi   = mean(yi),
            .groups   = "drop")

# Give each site a number, we will use it later to join tables
sites$site_id <- 1:nrow(sites)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Now with a proper projection. Equal Earth keeps the areas right, so Africa
# is not smaller than Greenland. EPSG:8857 is the code
crs_world <- "EPSG:8857"

################################################################################
# What does a coordinate represent?

# Step 1: coordinates as text
lat_text <- as.character(sites$lat)
lon_text <- as.character(sites$lon)

# Step 2: keep only what comes after the dot and count the characters.
# sub(".*\\.", "", x) means "delete everything up to the dot".
# If there is no dot (e.g. "30") there are 0 decimals
sites$dec_lat <- ifelse(grepl(".", lat_text, fixed = TRUE),
                        nchar(sub(".*\\.", "", lat_text)), 0)
sites$dec_lon <- ifelse(grepl(".", lon_text, fixed = TRUE),
                        nchar(sub(".*\\.", "", lon_text)), 0)

# Step 3: the less precise of the two, and the approximate error in km
sites$decimals     <- pmin(sites$dec_lat, sites$dec_lon)
sites$precision_km <- 111 / 10^sites$decimals

# Sites ordered from the least to the most precise
sites |>
  select(lat, lon, decimals, precision_km) |>
  arrange(decimals)

table(sites$decimals)


# 4. Countries and continents

# We only need the name and the continent of each country
world_simple <- world |>
  select(country = name, continent)

# Shorter names
world_simple$country[world_simple$country == "United States of America"] <- "USA"
world_simple$country[world_simple$country == "United Kingdom"] <- "UK"

# In which country is each site?
# st_nearest_feature gives the closest country. We use it instead of a normal
# "point inside polygon" join because coastal sites can fall in the sea
nearest <- st_nearest_feature(sites_sf, world_simple)

sites$country   <- world_simple$country[nearest]
sites$continent <- world_simple$continent[nearest]

# Natural Earth puts the whole of Russia in Europe. Our Russian site is in
# Siberia (east of the Urals, ~60 E), so we move it to Asia by hand
sites$continent[sites$country == "Russia" & sites$lon > 60] <- "Asia"

# Update the spatial version with the new columns
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Bring site_id, country and continent back to every effect size
es <- es |>
  left_join(select(sites, lat, lon, site_id, country, continent),
            by = c("lat", "lon"))


# 4a. How many studies per country and continent?
studies_country <- es |>
  group_by(country, continent) |>
  summarise(n_studies = n_distinct(Study_ID),
            n_es      = n(),
            .groups   = "drop") |>
  arrange(desc(n_studies))

studies_country

studies_continent <- studies_country %>%
  group_by(continent) %>%
  summarise(
    n_studies = sum(n_studies),
    .groups = "drop"
  )

studies_continent

# 5. Map the effect sizes: one dot per study

# The study is our unit of replication, so we want one dot per study.
# A few studies worked at more than one site. We keep one dot per study and
# site, because averaging their coordinates would put a dot where nobody worked
studies_sites <- es |>
  group_by(Study_ID, lat, lon, continent) |>
  summarise(n_es    = n(),
            mean_yi = mean(yi),
            .groups = "drop")

studies_sites_sf <- st_as_sf(studies_sites, coords = c("lon", "lat"),
                             crs = 4326, remove = FALSE)

summary(studies_sites_sf$n_es)

# We fix the colour and size limits here
lim_yi <- max(abs(studies_sites$mean_yi))
max_es <- max(studies_sites$n_es)

# Several studies were done at the same site and their dots overlap.
# st_jitter moves each dot a little (up to 60 km) so we can see all of them.

theme_set(theme_bw(base_size = 13))   # setting the font size

set.seed(1)

studies_world <- studies_sites_sf |>
  st_transform(crs_world) |>
  st_jitter(amount = 60000)          # metres, because Equal Earth is in metres


# Colours: red = negative, blue = positive, white = around zero
map_studies <- ggplot() +
  geom_sf(data = world, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = studies_world, aes(size = n_es, fill = mean_yi),
          shape = 21, colour = "grey40", stroke = 0.3, alpha = 1) +
  scale_fill_distiller(direction = -1,
                       limits = c(-lim_yi, lim_yi),
                       name = "Mean effect\nper study\n(SMD.H)") +
  scale_size_area(max_size = 8, limits = c(1, max_es),
                  name = "Effect sizes\nper study") +
  coord_sf(crs = crs_world)

map_studies_2 <- map_studies + theme(legend.position = c(0.955,0.69))   # same clean look for all plots

map_studies_2

# PNG, 600 dpi
ggsave("figures/Figure_1.png", plot = map_studies_2,
       width = 360, height = 180, units = "mm",
       dpi = 600, bg = "white")
