rm(list=ls())

# load packages
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# field site locations
field.sites <- data.frame(
  name = c("Site 3", "Site 2", "Site 1"),
  lat  = c(52.532, 52.413, 52.989),
  long = c(-0.192, -0.221, -0.912),
  ynudge = c(-0.35, 0.2, 0.25),
  xnudge = c(.55, .9, .5)
)

# Convert to sf object
sites_sf <- st_as_sf(field.sites, coords = c("long", "lat"), crs = 4326)

# Get UK boundary
uk <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")

# Plot
ggplot() +
  geom_sf(data = uk, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = sites_sf, size = 4, color="red", shape=18) +
  geom_sf_text(data = sites_sf, aes(label = name, nudge_y=ynudge, nudge_x=xnudge), size=7) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  coord_sf(xlim = c(-8, 2), ylim = c(49, 61)) +
  theme_void() +
  labs(x = "Longitude", y = "Latitude")

ggsave(filename="figs/15_sitemap.pdf", width=6, height=12, units="in") 
