pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape, sf,
               purrr, stringr, readr, lubridate, magrittr, janitor, openxlsx,
               geosphere, rotl)


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


EXP_m5 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|datapt_id,
                   ~ 1|study_id, 
                   ~ x_km + y_km | const
                 ), 
                 struct = "SPEXP", 
                 dist = "manhattan",
                 data = dat_coetzee
)

summary(EXP_m5)

# brms
# fit_1 <- bf(yi | se(sqrt(vi)) ~ 1 + 
#               (1 | datapt_id) + 
#               (1 | study_id) + 
#               gp(y_km, x_km, 
#                  cov = "exponential",
#                  scale = FALSE))
# 
# prior <- get_prior(formula = fit_1,
#                    data = dat_coetzee,
#                    family = gaussian()
# )
# 
# hist(dat_coetzee$x_km)
# 
# prior <- c(
#   prior(normal(0, 5), class = "Intercept"),
#   prior(student_t(3, 0, 2.5), class = "sd", group = "datapt_id"),
#   prior(student_t(3, 0, 2.5), class = "sd", group = "study_id"),
#   prior(exponential(0.001), class = "lscale"),  
#   prior(exponential(1), class = "sdgp")
# )


max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)
options(mc.cores = num_chains) 

fit_1 <- bf(
  yi | se(sqrt(vi)) ~ 1 + 
    (1 | datapt_id) + 
    (1 | study_id) + 
    gp(y_km, x_km, cov = "exponential", scale = FALSE)
  )

prior <- c(
    prior(normal(0, 5), class = "Intercept"),
    prior(student_t(3, 0, 2.5), class = "sd", group = "datapt_id"),
    prior(student_t(3, 0, 2.5), class = "sd", group = "study_id"),
    prior(lognormal(9, 1), class = "lscale"),   
    prior(student_t(3, 0, 0.5), class = "sdgp")
  )

m_exp <- brm(formula = fit_1,
             data = dat_coetzee,
             family = gaussian(),
            # prior = prior,
             iter = 11000,
             warmup = 9000,
             chains = num_chains,
             backend = "cmdstanr",
             threads = threading(threads_per_chain),
             control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(m_exp)


# Gaussian process model with square-exponential covariance function
fit_2 <- bf(
  yi | se(sqrt(vi)) ~ 1 + 
    (1 | datapt_id) + 
    (1 | study_id) + 
    gp(y_km, x_km, cov = "exp_quad", scale = FALSE)
)

m_gau <- brm(formula = fit_2,
             data = dat_coetzee,
             family = gaussian(),
             # prior = prior,
             iter = 11000,
             warmup = 9000,
             chains = num_chains,
             backend = "cmdstanr",
             threads = threading(threads_per_chain),
             control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(m_gau)

GAU_1 <-  rma.mv(yi, vi, 
                   random = list(
                     ~ 1|datapt_id,
                     ~ 1|study_id, 
                     ~ datapt_id|const
                   ), 
                   struct = "SPGAU", 
                   data = dat_coetzee,
                   dist = list(datapt_id = dist_matrix_euclid)
  )

summary(GAU_1)


GAU_2 <- rma.mv(
  yi, vi,
  random = list(
    ~ 1|datapt_id,
    ~ 1|study_id,
    ~ x_km + y_km | const
  ),
  struct = "SPGAU",
  data = dat_coetzee
)

summary(GAU_2)
