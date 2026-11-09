##adding ecoregion

library(sf)
library(dplyr)
library(here)


##1. Loading shapefile
ecoregions <- st_read(here("Ecoregions", "ecoregions.shp"))
ecoregions_clean <- ecoregions |>
  select(ECOREGION, REGION_NAM)

##2 matching crs
points_sf <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)

# 3. Ensure both layers are in the same crs
if (st_crs(points_sf) != st_crs(ecoregions_clean)) {
  points_sf <- st_transform(points_sf, st_crs(ecoregions_clean))
}

# 4. Perform the spatial join
joined_points <- st_join(points_sf, ecoregions_clean, join = st_intersects)

# 5. Convert back to a regular data frame for your GLM
dat <- as.data.frame(joined_points) %>% select(-geometry)

# 6. Explicitly tell R it's a categorical factor variable
dat$ecoregion <- as.factor(dat$ecoregion)





