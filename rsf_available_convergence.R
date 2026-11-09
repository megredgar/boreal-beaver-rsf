# How many available points? Convergence check on the GEE export.
# Fits the global model at increasing N (random draws from the 50k pool) and
# plots each coefficient vs N. Pick the smallest N where every slope has
# flattened and the draw-to-draw spread is small relative to its SE.

library(dplyr)
library(ggplot2)

covs <- c("dist_water", "dist_wetland", "dist_deciduous",
          "dist_coniferous", "dist_developed", "dist_roads")

used  <- filter(beaver_RSF_covariates, used == 1)
avail <- filter(beaver_RSF_covariates, used == 0)

N_grid <- c(1000, 2500, 5000, 10000, 20000, 35000, 50000, 60000)
n_rep  <- 10   # independent draws per N (except N = full pool)

fit_once <- function(N, rep) {
  a  <- slice_sample(avail, n = min(N, nrow(avail)))
  dd <- bind_rows(used, a)
  # standardize on the fitted data (training mean/SD), as in the proposal
  dd[covs] <- scale(dd[covs])
  # large weight on available points -> slopes approximate the IPP
  # (Warton & Shepherd 2010; Fithian & Hastie 2013); intercept is not interpretable
  dd$w <- ifelse(dd$used == 1, 1, 1000)
  m <- glm(reformulate(covs, "used"), family = binomial, data = dd, weights = w)
  data.frame(N = N, rep = rep, term = names(coef(m))[-1], est = coef(m)[-1])
}

set.seed(2026)
res <- bind_rows(lapply(N_grid, function(N)
  bind_rows(lapply(seq_len(ifelse(N >= nrow(avail), 1, n_rep)),
                   function(r) fit_once(N, r)))))

ggplot(res, aes(N, est)) +
  geom_point(alpha = 0.4) +
  stat_summary(fun = mean, geom = "line") +
  facet_wrap(~term, scales = "free_y") +
  scale_x_log10() +
  labs(x = "Number of available points", y = "Standardized coefficient") +
  theme_bw()
