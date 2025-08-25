

#' @test

# for comparison the results between metafor and brms, to use scaled latitude and longtude for metafor
# and scale = FALSE in brms…


dat_coetzee <- read.xlsx(here("data", "examples", "Coetzee_2014.xlsx"), sheet = 2)
dat_coetzee$const <- 1 # add constant for spatial models

# make distance matrix
coords <- cbind(dat_coetzee$long, dat_coetzee$lat)
dist_matrix <- distm(coords, fun = distHaversine)
dist_matrix[1:5,1:5]

rownames(dist_matrix) <- dat_coetzee$datapt_id
colnames(dist_matrix) <- dat_coetzee$datapt_id


###
# GP model
###

# metafor
### non-scaled
GAU_maximum <- rma.mv(yi, vi, 
                   random = list(
                     ~ 1|datapt_id,
                     ~ 1|study_id, 
                     ~ lat + long | const
                   ), 
                   struct="SPGAU", 
                   dist = "maximum",
                   data = dat_coetzee)

summary(GAU_maximum)

GAU_manhattan <- rma.mv(yi, vi, 
                      random = list(
                        ~ 1|datapt_id,
                        ~ 1|study_id, 
                        ~ lat + long | const
                      ), 
                      struct="SPGAU", 
                      dist = "manhattan",
                      data = dat_coetzee)

summary(GAU_manhattan)


### scaled

dat_coetzee$long_s <- scale(dat_coetzee$long, center = TRUE, scale = TRUE)
dat_coetzee$lat_s <- scale(dat_coetzee$lat, center = TRUE, scale = TRUE)


EXP_mod1.1 <- rma.mv(yi, vi, 
                     random = list(
                       ~ 1|datapt_id,
                       ~ 1|study_id, 
                       ~ lat_s + long_s | const # changed
                     ), 
                     struct="SPGAU", 
                     data = dat_coetzee,
                     # dist=list(dist_matrix),
                     # control=list(rho.init = 10)
)

summary(EXP_mod1.1)

## brms

### scaled
fit_2 <- bf(yi | se(sqrt(vi)) ~ 1 + 
              (1 | datapt_id) + 
              (1 | study_id) + 
              gp(long, lat, cov = "exp_quad"))

prior2 <- get_prior(formula = fit_2,
                    data = dat_coetzee,
                    family = gaussian()
)

m_gp <- brm(formula = fit_2,
            data = dat_coetzee,
            family = gaussian(),
            prior = prior2,
            iter = 15000,
            warmup = 13000,
            chain = 2, 
            thin = 1
)
summary(m_gp)

### non-scaled
fit_1.1 <- bf(yi | se(sqrt(vi)) ~ 1 + 
                (1 | datapt_id) + 
                (1 | study_id) + 
                gp(long, lat, 
                   cov = "exp_quad",
                   scale = FALSE) # added
) 

prior1.1 <- get_prior(formula = fit_2.1,
                      data = dat_coetzee,
                      family = gaussian()
)

m_exp1.1 <- brm(formula = fit_2.1,
                data = dat_coetzee,
                family = gaussian(),
                prior = prior2.1,
                iter = 8000,
                warmup = 6000,
                chain = 2, 
                thin = 1
)

summary(m_exp1.1)


library(sf)
pts <- st_as_sf(dat_coetzee, coords = c("long", "lat"), crs = 4326)     # WGS-84
pts <- st_transform(pts, 3857)                                  # Web-Merc (m)
dat_coetzee$x_km <- st_coordinates(pts)[,1] / 1000                      # km
dat_coetzee$y_km <- st_coordinates(pts)[,2] / 1000

fit <- brm( yi | se(sqrt(vi)) ~ 1 + (1|datapt_id) + (1|study_id) +
              gp(x_km, y_km, cov = "exp_quad", scale = FALSE),
            data = dat_coetzee,
            prior = c(prior(student_t(3,0,5), class = "sdgp"),
                      prior(gamma(2,0.002),   class = "lscale")),
            iter = 6000,
            warmup = 4000,
            thin = 1)
summary(fit)