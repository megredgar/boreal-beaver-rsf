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

##loading csvs##
used_covariates <- read.csv(here("covariates", "beaver_RSF_covs_used.csv"))

available_covariates <- read.csv(here("covariates", "beaver_RSF_covariates.csv"))

dat  <- bind_rows(used_covariates, available_covariates)

covs <- c("dist_wetland","dist_water","dist_fen", "dist_bog", "dist_marsh",
          "dist_swamp", "dist_deciduous","dist_coniferous","dist_developed", "dist_roads")

# log-transform, then standardize, and SAVE the means/sds
dat[covs] <- lapply(dat[covs], log1p)
sc <- lapply(dat[covs], \(x) c(mean = mean(x), sd = sd(x)))
dat[covs] <- lapply(covs, \(v) (dat[[v]] - sc[[v]]["mean"]) / sc[[v]]["sd"])
saveRDS(sc, "scaling.rds")

## checking for collinearity
cor(dat[vars])

##global model
m_global <- glmmTMB(used ~ dist_wetland + dist_water +
                      dist_deciduous + dist_coniferous + dist_developed + dist_roads,
                    family = binomial, data = dat)

##global model without developed and roads
m_global_noinf <- glmmTMB(used ~ dist_wetland + dist_water +
                            dist_deciduous * dist_coniferous,
                          family = binomial, data = dat)
summary(m_global_noinf)
bbmle::AICtab(m_global, m_global_noinf)

summary(m_global)
performance::check_collinearity(m_global)

##infrastructure model
i_model <- glmmTMB(used ~ dist_developed + dist_roads,
                   family = binomial, data = dat)
summary(i_model)

##deciduous forest model
d_model <- glmmTMB(used ~ dist_deciduous,
                   family = binomial, data = dat)
summary(d_model)

##dist to water and wetland##
w_model <- glmmTMB(used ~ dist_water + dist_wetland,
                   family = binomial, data = dat)
summary(w_model)

##deciduous forest + water
dw_model <- glmmTMB(used ~ dist_water + dist_deciduous,
                    family = binomial, data = dat)
summary(dw_model)

##deciduous forest * water
dfw_model <- glmmTMB(used ~ dist_water + dist_wetland*dist_deciduous,
                     family = binomial, data = dat)
summary(dfw_model)

bbmle::AICtab(d_model, dfw_model, dw_model, i_model, w_model, m_global)

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

bbmle::AICtab(m_null, m_wet, m_swamp, m_marsh, m_fen, m_bog, weights = TRUE)
