##adding ecoregion

library(sf)
library(dplyr)
library(here)

# 1. Convert your RSF points data frame into a spatial object
# (Replace X and Y with your actual column names. Set crs to match your points, e.g., 4326 for WGS84)
points_sf <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)

# 2. Read in your ecoregion shapefile
ecoregions <- st_read(here("Ecoregions", "ecoregions.shp")

# 3. Ensure both layers are in the EXACT same Coordinate Reference System (CRS)
if (st_crs(points_sf) != st_crs(ecoregions)) {
  points_sf <- st_transform(points_sf, st_crs(ecoregions))
}

# 4. Perform the spatial join (extracts ecoregion attributes to points)
joined_points <- st_join(points_sf, ecoregions, join = st_intersects)

# 5. Convert back to a regular data frame for your GLM
dat <- as.data.frame(joined_points) %>% select(-geometry)

# 6. Explicitly tell R it's a categorical factor variable
dat$ecoregion <- as.factor(dat$ecoregion)



