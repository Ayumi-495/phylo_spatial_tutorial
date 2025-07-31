pacman::p_load(brms, MCMCglmm, metafor, metadat, tidyverse, data.table, crayon, here, ape,
               dplyr, tidyr, purrr, stringr, readr, readxl, lubridate, magrittr, janitor, 
               flextable, kableExtra, geosphere, rotl, data.table, openxlsx)


dat_coetzee <- read.xlsx(here("data", "examples", "Coetzee_2014.xlsx"), sheet = 2)
dat_coetzee$const <- 1 # add constant for spatial models

coords <- cbind(dat_coetzee$long, dat_coetzee$lat)
ref <- coords[1, ]

# dat_coetzee$x_km <- distGeo(cbind(ref[1], dat_coetzee$lat), coords) / 1000
# dat_coetzee$x_km[dat_coetzee$long < ref[1]] <- -dat_coetzee$x_km[dat_coetzee$long < ref[1]]
# 
# dat_coetzee$y_km <- distGeo(cbind(dat_coetzee$long, ref[2]), coords) / 1000
# dat_coetzee$y_km[dat_coetzee$lat < ref[2]] <- -dat_coetzee$y_km[dat_coetzee$lat < ref[2]]
# 
# coords_km <- cbind(dat_coetzee$x_km, dat_coetzee$y_km)
# dist_matrix_euclid <- as.matrix(dist(coords_km))

###
# Euclidean distance matrix
###

# project to a planar coordinate system and convert to kilometers
dat_sf <- st_as_sf(dat_coetzee, coords = c("long", "lat"), crs = 4326)
dat_sf_proj <- st_transform(dat_sf, crs = 3857)

coords_m <- st_coordinates(dat_sf_proj)
dat_coetzee$x_km <- coords_m[,1] / 1000
dat_coetzee$y_km <- coords_m[,2] / 1000


coords_km <- cbind(dat_coetzee$x_km, dat_coetzee$y_km)
dist_matrix_euclid <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid) <- dat_coetzee$datapt_id
colnames(dist_matrix_euclid) <- dat_coetzee$datapt_id

###
# spherical
###
coords_ll <- cbind(dat_coetzee$long, dat_coetzee$lat)

dist_matrix_spherical_m <- distm(coords_ll, fun = distGeo)
dist_matrix_spherical_km <- dist_matrix_spherical_m / 1000

rownames(dist_matrix_spherical_km) <- dat_coetzee$datapt_id
colnames(dist_matrix_spherical_km) <- dat_coetzee$datapt_id
 
# Euclidean
EXP_m1 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|datapt_id,
                   ~ 1|study_id, 
                   ~ datapt_id|const
                 ), 
                 struct = "SPEXP", 
                 data = dat_coetzee,
                 dist = list(datapt_id = dist_matrix_euclid)
)
   
summary(EXP_m1)

# scale to [0,1] to match metefor internal
dist_matrix_scaled <- dist_matrix_euclid / max(dist_matrix_euclid)

# rerun model
EXP_m1_scaled <- rma.mv(yi, vi, 
                        random = list(
                          ~ 1|datapt_id,
                          ~ 1|study_id, 
                          ~ datapt_id|const
                        ), 
                        struct = "SPEXP", 
                        data = dat_coetzee,
                        dist = list(datapt_id = dist_matrix_scaled)
)

summary(EXP_m1_scaled)

# spherical
EXP_m2 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|datapt_id,
                   ~ 1|study_id, 
                   ~ datapt_id|const
                 ), 
                 struct = "SPEXP", 
                 data = dat_coetzee,
                 dist = list(datapt_id = dist_matrix_spherical_km)
)

summary(EXP_m2)

# use default (Euclid)
EXP_m3 <- rma.mv(yi, vi, 
                   random = list(
                     ~ 1|datapt_id,
                     ~ 1|study_id, 
                     ~ lat + long | const
                   ), 
                   struct = "SPEXP", 
                   data = dat_coetzee
                   )
summary(EXP_m3)

EXP_m3.1 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|datapt_id,
                   ~ 1|study_id, 
                   ~ long + lat | const
                 ), 
                 struct = "SPEXP", 
                 data = dat_coetzee
)
summary(EXP_m3.1)


EXP_m3_2 <- rma.mv(
  yi, vi,
  random = list(
    ~ 1|datapt_id,
    ~ 1|study_id,
    ~ x_km + y_km | const
  ),
  struct = "SPEXP",
  data = dat_coetzee
)

summary(EXP_m3_2)

# use dist = "gcd" for the great-circle distance (WGS84 ellipsoid method)
EXP_m4 <- rma.mv(yi, vi, 
                   random = list(
                     ~ 1|datapt_id,
                     ~ 1|study_id, 
                     ~ lat + long | const
                   ), 
                   struct = "SPEXP", 
                   dist = "gcd",
                   data = dat_coetzee
                 )

summary(EXP_m4)

# brms
fit_1 <- bf(yi | se(sqrt(vi)) ~ 1 + 
              (1 | datapt_id) + 
              (1 | study_id) + 
              gp(y_km, x_km, 
                 cov = "exponential",
                 scale = FALSE))

# prior <- get_prior(formula = fit_1,
#                    data = dat_coetzee,
#                    family = gaussian()
# )

# prior <- c(
#   prior(normal(0, 1), class = "Intercept"),
#   prior(exponential(1), class = "sdgp"),
#   prior(lognormal(log(1000), 0.5), class = "lscale")

  prior <- c(
    prior(exponential(1), class = "sdgp", coef = "gpy_kmx_km"),
    prior(lognormal(log(1000), 0.5), class = "lscale", coef = "gpy_kmx_km")
  )

m_exp <- brm(formula = fit_1,
             data = dat_coetzee,
             family = gaussian(),
             # prior = prior,
             iter = 8000,
             warmup = 6000,
             chain = 2, 
             thin = 1
)

summary(m_exp)
Warning messages:
  1: In pgamma(1/q, shape, rate = rate, lower.tail = !lower.tail, log.p = log.p) :
  NaNs produced
2: In pgamma(1/q, shape, rate = rate, lower.tail = !lower.tail, log.p = log.p) :
  NaNs produced
3: There were 61 divergent transitions after warmup. See
https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
to find out why this is a problem and how to eliminate them. 
4: Examine the pairs() plot to diagnose sampling problems


> summary(m_exp)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: yi | se(sqrt(vi)) ~ 1 + (1 | datapt_id) + (1 | study_id) + gp(y_km, x_km, cov = "exponential", scale = FALSE) 
# Data: dat_coetzee (Number of observations: 1484) 
# Draws: 2 chains, each with iter = 8000; warmup = 6000; thin = 1;
# total post-warmup draws = 4000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gpy_kmx_km)       0.15      0.12     0.01     0.43 1.00      459      600
# lscale(gpy_kmx_km)  3850.51  10567.30   377.47 19974.97 1.00     1509     1720
# 
# Multilevel Hyperparameters:
#   ~datapt_id (Number of levels: 1484) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.25      0.03     1.19     1.32 1.00      882     1819
# 
# ~study_id (Number of levels: 127) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.86      0.09     0.69     1.03 1.01      685     1365
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.45      0.12     0.22     0.69 1.00      819      664
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
#   There were 61 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
