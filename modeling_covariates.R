## modelling

##installing packages
install.packages("glmmTMB")
install.packages("performance")
install.packages("bbmle")
install.packages("ggeffects")
install.packages("DHARMa")
install.packages("geosphere")
##library
library(here)
library(dplyr)
library(glmmTMB)
library(performance)
library(bbmle)
library(ggeffects)
library(ggplot2)
library(DHARMa)
library(geosphere)

##loading csvs##
used_covariates <- read.csv(here("covariates", "beaver_RSF_covs_used.csv"))

used_available_covariates <- read.csv(here("covariates", "beaver_RSF_covariates.csv"))

covs <- c("dist_wetland","dist_water","dist_fen", "dist_bog", "dist_marsh",
          "dist_swamp", "dist_deciduous","dist_coniferous","dist_developed", "dist_roads")

# log-transform, then standardize, and SAVE the means/sds
used_available_covariates[covs] <- lapply(used_available_covariates[covs], log1p)
sc <- lapply(used_available_covariates[covs], \(x) c(mean = mean(x), sd = sd(x)))
used_available_covariates[covs] <- lapply(covs, \(v) (used_available_covariates[[v]] - sc[[v]][["mean"]]) / sc[[v]]["sd"])
saveRDS(sc, "scaling.rds")

##checking for collinearity
cor(used_available_covariates[covs])

##1. global model
m_global <- glmmTMB(used ~ dist_fen + dist_bog + dist_swamp + dist_marsh + dist_water +
                      dist_deciduous + dist_coniferous + dist_developed + dist_roads + 
                     (1|REGION_NAM) ,
                    family = binomial, data = used_available_covariates)

summary(m_global)
performance::check_collinearity(m_global)

##performing moran's i test to test for spatial autocorrelation after including ecoregion as a random effect
set.seed(6)

sample_dat <- used_available_covariates |>
       slice_sample(n = 10000)
 
m_sample <- glmmTMB(
         used ~ dist_fen + dist_bog + dist_swamp + dist_marsh +
              dist_water + dist_deciduous + dist_coniferous +
               dist_developed + dist_roads + (1 | REGION_NAM),
         family = binomial,
         data = sample_dat
     )
simulationOutput <- simulateResiduals(
  m_sample,
  n = 250
)
    
 dist_matrix <- distm(
      cbind(used_available_covariates$lon, used_available_covariates$lat),
      fun = distHaversine
    )
 testSpatialAutocorrelation(
      simulationOutput,
      x = used_available_covariates$lon,
      y = used_available_covariates$lat,
      distMat = dist_matrix
    )
   
##global model without developed and roads
m_global_noinf <- glmmTMB(used ~ dist_wetland + dist_water +
                            dist_deciduous * dist_coniferous + (1|REGION_NAM),
                          family = binomial, data = used_available_covariates)
summary(m_global_noinf)

bbmle::AICtab(m_global, m_global_noinf)



##2. infrastructure model
i_model <- glmmTMB(used ~ dist_developed + dist_roads + (1|REGION_NAM),
                   family = binomial, data = used_available_covariates)
summary(i_model)


##3. deciduous forest model
d_model <- glmmTMB(used ~ dist_deciduous +  (1|REGION_NAM),
                   family = binomial, data = used_available_covariates)
summary(d_model)


##4. dist to water and wetland
w_model <- glmmTMB(used ~ dist_water + dist_wetland + (1|REGION_NAM),
                   family = binomial, data = used_available_covariates)
summary(w_model)


##5. deciduous forest * water
dfw_model <- glmmTMB(used ~ dist_water*dist_deciduous + (1|REGION_NAM),
                     family = binomial, data = used_available_covariates)
summary(dfw_model)

bbmle::AICtab(d_model, dfw_model, i_model, w_model, m_global)


##distance to wetland
wetland_m <- glmmTMB(used ~ dist_wetland + (1|REGION_NAM),
                     family = binomial, data = used_available_covariates)

##distance to wetland types
wetlandtypes_m <- glmmTMB(used ~ dist_fen + dist_bog + dist_marsh + dist_swamp + (1|REGION_NAM),
                          family = binomial, data = used_available_covariates)

bbmle::AICtab(wetland_m, wetlandtypes_m, weights = TRUE)

##distance to each wetland type 
m_swamp <- glmmTMB(used ~ dist_swamp + (1|REGION_NAM), binomial, data = used_available_covariates)
m_marsh <- glmmTMB(used ~ dist_marsh + (1|REGION_NAM), binomial, data = used_available_covariates)
m_fen   <- glmmTMB(used ~ dist_fen + (1|REGION_NAM), binomial, data = used_available_covariates)
m_bog   <- glmmTMB(used ~ dist_bog + (1|REGION_NAM), binomial, data = used_available_covariates)
m_wet <- glmmTMB(used ~ dist_wetland + (1|REGION_NAM), binomial, data = used_available_covariates)
m_null  <- glmmTMB(used ~ 1 + (1|REGION_NAM), binomial, data = used_available_covariates)
water_model <- glmmTMB(used ~ dist_water + (1|REGION_NAM), binomial, data = used_available_covariates)

AICtab(m_null, m_wet, m_swamp, m_marsh, m_fen, m_bog, water_model, weights = TRUE)

##testing whether water or wetland is overall more important
AICtab(water_model, wetland_m, weights = TRUE)

##global model without wetland 
m_globalwat <- glmmTMB(used ~ dist_water + dist_deciduous + dist_coniferous + dist_developed
                       + dist_roads + (1|REGION_NAM),
                       family = binomial, data = used_available_covariates)
AICtab(m_global, m_globalwat)


#plot showing the relationship between distance to deciduous forest and distance to water 
predictions <- predict_response(
  dfw_model, 
  terms = c("dist_deciduous [all]", "dist_water [-1, 0, 2]")
)

# Plot the curves
plot(predictions) +
  labs(
    title = "Beaver Selection: Deciduous Forest Effect by Distance to Water",
    x = "Distance to Deciduous Forest (Log-Scaled)",
    y = "Probability of Selection",
    color = "Distance to Water"
  ) +
  theme_minimal()

