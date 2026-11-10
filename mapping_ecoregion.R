##adding ecoregion

library(sf)
library(dplyr)
library(here)
library(terra)


##1. Loading shapefile
ecoregions <- st_read(here("Ecoregions", "ecoregions.shp"))
ecoregions_clean <- ecoregions |>
  select(ECOREGION, REGION_NAM)

##2 matching crs
points_sf <- st_as_sf(used_available_covariates, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# 3. Ensure both layers are in the same crs
if (st_crs(points_sf) != st_crs(ecoregions_clean)) {
  points_sf <- st_transform(points_sf, st_crs(ecoregions_clean))
}

# 4. Perform the spatial join
joined_points <- st_join(points_sf, ecoregions_clean, join = st_intersects)

# 5. Convert back to a regular data frame for your GLM
used_available_covariates <- as.data.frame(joined_points) %>% select(-geometry)

# 6. Explicitly tell R it's a categorical factor variable
used_available_covariates$REGION_NAM <- as.factor(used_available_covariates$REGION_NAM)

#7. determining what ecoregion beavers most occur in based on unobserved points
##setting the infrastructure columns to 0 which is the average and removes the observer effects
newd <- used_available_covariates
newd$dist_developed <- 0
newd$dist_roads     <- 0

used_available_covariates$pred <- predict(m_global, newdata = newd, type = "response") ##asking what would the predictive relative 
                                                                 ##likelihood be if the distance to roads and developed were zero at every location


used_available_covariates |>
  filter(used == 0) |>     ##filtering to unobserved               
  group_by(REGION_NAM) |>  ##grouping by region name
  summarise(n = n(),       
            mean_pred   = mean(pred),
            median_pred = median(pred)) |> ##calculating mean and median predicted occurrence likelihood
  arrange(desc(mean_pred)) ##comparing across regions

##plotting
used_available_covariates |> filter(used == 0) |>
  ggplot(aes(reorder(REGION_NAM, pred, median), pred)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = NULL, y = "Predicted relative likelihood of occurrence")



