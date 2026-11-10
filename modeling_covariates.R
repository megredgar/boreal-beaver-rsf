## modelling

##installing packages
install.packages("glmmTMB")
install.packages("performance")
install.packages("bbmle")
install.packages("ggeffects")
##library
library(here)
library(dplyr)
library(glmmTMB)
library(performance)
library(bbmle)
library(ggeffects)
library(ggplot2)

##loading csvs##
used_covariates <- read.csv(here("covariates", "beaver_RSF_covs_used.csv"))

used_available_covariates <- read.csv(here("covariates", "beaver_RSF_covariates.csv"))

dat  <- bind_rows(used_covariates, used_available_covariates)

covs <- c("dist_wetland","dist_water","dist_fen", "dist_bog", "dist_marsh",
          "dist_swamp", "dist_deciduous","dist_coniferous","dist_developed", "dist_roads")

# log-transform, then standardize, and SAVE the means/sds
dat[covs] <- lapply(dat[covs], log1p)
sc <- lapply(dat[covs], \(x) c(mean = mean(x), sd = sd(x)))
dat[covs] <- lapply(covs, \(v) (dat[[v]] - sc[[v]][["mean"]]) / sc[[v]]["sd"])
saveRDS(sc, "scaling.rds")

##checking for collinearity
cor(dat[covs])

##1. global model
m_global <- glmmTMB(used ~ dist_wetland + dist_water +
                      dist_deciduous + dist_coniferous + dist_developed + dist_roads,
                    family = binomial, data = dat)

summary(m_global)
performance::check_collinearity(m_global)

##global model without developed and roads
m_global_noinf <- glmmTMB(used ~ dist_wetland + dist_water +
                            dist_deciduous * dist_coniferous,
                          family = binomial, data = dat)
summary(m_global_noinf)

bbmle::AICtab(m_global, m_global_noinf)



##2. infrastructure model
i_model <- glmmTMB(used ~ dist_developed + dist_roads,
                   family = binomial, data = dat)
summary(i_model)


##3. deciduous forest model
d_model <- glmmTMB(used ~ dist_deciduous,
                   family = binomial, data = dat)
summary(d_model)


##4. dist to water and wetland
w_model <- glmmTMB(used ~ dist_water + dist_wetland,
                   family = binomial, data = dat)
summary(w_model)


##5. deciduous forest * water
dfw_model <- glmmTMB(used ~ dist_water*dist_deciduous,
                     family = binomial, data = dat)
summary(dfw_model)

bbmle::AICtab(d_model, dfw_model, i_model, w_model, m_global)


##distance to wetland
wetland_m <- glmmTMB(used ~ dist_wetland,
                     family = binomial, data = dat)

##distance to wetland types
wetlandtypes_m <- glmmTMB(used ~ dist_fen + dist_bog + dist_marsh + dist_swamp,
                          family = binomial, data = dat)

bbmle::AICtab(wetland_m, wetlandtypes_m, weights = TRUE)

##distance to each wetland type 
m_swamp <- glmmTMB(used ~ dist_swamp , binomial, data = dat)
m_marsh <- glmmTMB(used ~ dist_marsh , binomial, data = dat)
m_fen   <- glmmTMB(used ~ dist_fen, binomial, data = dat)
m_bog   <- glmmTMB(used ~ dist_bog, binomial, data = dat)
m_wet <- glmmTMB(used ~ dist_wetland, binomial, data = dat)
m_null  <- glmmTMB(used ~ 1, binomial, data = dat)
water_model <- glmmTMB(used ~ dist_water, binomial, data = dat)

bbmle::AICtab(m_null, m_wet, m_swamp, m_marsh, m_fen, m_bog, water_model, weights = TRUE)

##testing whether water or wetland is overall more important
water_model <- glmmTMB(used ~ dist_water,
                   family = binomial, data = dat)

bbmle::AICtab(water_model, wetland_m, weights = TRUE)

##global model without wetland 
m_globalwat <- glmmTMB(used ~ dist_water + dist_deciduous + dist_coniferous + dist_developed
                       + dist_roads,
                       family = binomial, data = dat)
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

