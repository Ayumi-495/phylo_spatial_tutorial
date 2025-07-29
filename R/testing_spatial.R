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

204766.3216 
22097415.7402
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
             prior = prior,
             iter = 2000,
             warmup = 1000,
             chain = 2, 
             thin = 1
)

