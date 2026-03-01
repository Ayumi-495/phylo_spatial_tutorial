dat_scholer2020 <- read.csv("data/Scholer_2020/Scholer_2020.csv")

head(dat_scholer2020)
# logit_survival -> yi, 1/se -> vi^2

dat_scholer2020 <- dat_scholer2020 |> 
  mutate(vi = se^2)

dat_scholer2020$const <- 1 # add constant for spatial models
dat_scholer2020$effect_id <- seq_len(nrow(dat_scholer2020))

coords2 <- cbind(dat_scholer2020$long, dat_scholer2020$lat)

# project to a planar coordinate system and convert to kilometers
dat_sf2 <- st_as_sf(dat_scholer2020, coords = c("long", "lat"), crs = 4326)
dat_sf2_proj <- st_transform(dat_sf2, crs = 3857)

coords_m <- st_coordinates(dat_sf2_proj)
dat_scholer2020$x_km <- coords_m[,1] / 1000
dat_scholer2020$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_scholer2020$x_km, dat_scholer2020$y_km)
dist_matrix_euclid2 <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid2) <- dat_scholer2020$effect_id
colnames(dist_matrix_euclid2) <- dat_scholer2020$effect_id

dat_scholer2020 <- dat_scholer2020 |>
  group_by(lat, long) |>
  mutate(site_id = cur_group_id()) |>
  ungroup()
names(dat_scholer2020)
dat_scholer2020 <- dat_scholer2020 |>
  mutate(site_id = as.factor(site_id),
         effect_id = as.factor(effect_id))
nrow(dat_scholer2020)
# [1] 949

cat_vars <- dat_scholer2020 |>
  dplyr::select(where( ~ is.factor(.x) || is.character(.x)))

n_levels <- cat_vars |> 
  dplyr::summarise(across(
    everything(),
    ~ dplyr::n_distinct(.x)
  )) |> 
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_levels"
  )

n_levels
# A tibble: 15 × 2
# variable     n_levels
# <chr>           <int>
# 1 ref               205
# 2 order              15
# 3 species           636
# 4 common            636
# 5 tip_label         625
# 6 migration           2
# 7 sitename          207
# 8 country            46
# 9 hemisphere          2
# 10 passerine           2
# 11 realm               2
# 12 study_type          4
# 13 auk_used            2
# 14 geomean_used        9
# 15 effect_id         949
# 16 site_id           454

ma_sp1_exp0_mf <- rma.mv(logit_survival, vi, 
                         random = list(
                           ~ 1|effect_id, 
                           ~ 1|ref,
                           ~ 1|site_id,                           
                           ~ effect_id|const
                         ), 
                         struct = "SPEXP", 
                         data = dat_scholer2020,
                         dist = list(site = dist_matrix_euclid2),
                         test = "t",
                         method = "REML",
                         sparse = TRUE,
                         verbose = TRUE)
summary(ma_sp1_exp0_mf)


sp_eg1_exp_mf <- rma.mv(logit_survival, vi,
                                    random = list(
                                      ~ 1|effect_id, 
                                      ~ 1|ref,
                                      # ~ 1|site_id, 
                                      ~ x_km + y_km |const
                                    ), 
                                    struct = "SPEXP", 
                                    data = dat_scholer2020,
                                    sparse = TRUE,
                                    verbose = TRUE,
                                    method = "REML",
                                    test = "t")
summary(sp_eg1_exp_mf)


# site id version ---------------------------------------------------------
#' [this is better as this dataset includes multiple effect size per location]
site_coords <- dat_scholer2020 |>
  distinct(site_id, x_km, y_km) |>
  arrange(site_id)

coords_site <- as.matrix(site_coords[, c("x_km", "y_km")])
dist_site <- as.matrix(dist(coords_site))

rownames(dist_site) <- as.character(site_coords$site_id)
colnames(dist_site) <- as.character(site_coords$site_id)
dist_site

dat_Roger$site_id <- as.character(dat_Roger$site_id)

## metafor ------------------------------------------
ma_sp_site <- rma.mv(
  yi = logit_survival,
  V  = vi,                      
  random = list(
    ~ 1 | effect_id,          
    ~ 1 | ref,                   
    # ~ 1 | site_id,               
    ~ site_id | const           
  ),
  struct = "SPEXP",
  data = dat_scholer2020,
  dist = list(site = dist_site), 
  method = "REML",
  test = "t",
  sparse = TRUE,
  verbose = TRUE
)

ma_sp_site


sp_eg1_exp_mf <- rma.mv(logit_survival, vi,
                        random = list(
                          ~ 1|effect_id, 
                          ~ 1|ref,
                          # ~ 1|site_id, 
                          ~ x_km + y_km | const
                        ), 
                        struct = "SPEXP", 
                        data = dat_scholer2020,
                        sparse = TRUE,
                        verbose = TRUE,
                        method = "REML",
                        test = "t")
summary(sp_eg1_exp_mf)

# >>> Gaussian kernel version
ma_sp_site_gau <- rma.mv(
  yi = logit_survival,
  V  = vi,                      
  random = list(
    ~ 1 | effect_id,          
    ~ 1 | ref,                   
    # ~ 1 | site_id,               
    ~ site_id | const           
  ),
  struct = "SPGAU",
  data = dat_scholer2020,
  dist = list(site = dist_site), 
  method = "REML",
  test = "t",
  sparse = TRUE,
  verbose = TRUE
)

summary(ma_sp_site_gau)

## brms ------------------------------------------
### exponential ver
dat_scholer2020$effect_id <- as.character(dat_scholer2020$effect_id)
vcv <- diag(dat_scholer2020$vi)

rownames(vcv) <- colnames(vcv) <- dat_scholer2020$effect_id

fit_gp_brms <- bf(logit_survival ~ 1 +
                    (1|ref) + # this is site-level random effect
                    (1|gr(effect_id, cov = vcv)) + # this is m (sampling error)
                    gp(x_km, y_km, 
                       cov = "exponential",
                       scale = FALSE)
)

# generate default priors
prior <- default_prior(fit_gp_brms, 
                       data = dat_scholer2020, 
                       data2 = list(vcv = vcv),
                       family = gaussian())
prior$prior[5] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior

sp_exp_brms <- brm(
  formula = fit_gp_brms,
  data = dat_scholer2020,
  data2 = list(vcv = vcv),
  prior = prior,
  iter = 2000, 
  warmup = 1000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
summary(sp_exp_brms)


### gaussian ver
ma_brms_gau <- bf(logit_survival | se(sqrt(vi)) ~ 1 +
                      (1 | ref) +
                      (1| effect_id) + 
                      gp(x_km, y_km, scale = FALSE))

prior_brms <- get_prior(
  formula = ma_brms_gau,
  data = dat_scholer2020,
  family = gaussian()
)

fit_gp <- brm(
  formula = ma_brms_gau,
  prior = prior_brms,
  data = dat_scholer2020,
  family = gaussian(),
  iter = 2000, 
  warmup = 1000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)


dat_scholer2020$effect_id <- as.character(dat_scholer2020$effect_id)
vcv <- diag(dat_scholer2020$vi)

rownames(vcv) <- colnames(vcv) <- dat_scholer2020$effect_id

fit_gp_brms <- bf(logit_survival ~ 1 +
                    (1|ref) + # this is site-level random effect
                    (1|gr(effect_id, cov = vcv)) + # this is m (sampling error)
                    gp(x_km, y_km, 
                       cov = "exp_quad",
                       scale = FALSE)
)

# generate default priors
prior <- default_prior(fit_gp_brms, 
                       data = dat_scholer2020, 
                       data2 = list(vcv = vcv),
                       family = gaussian())
prior$prior[5] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior

sp_exp_quad_brms <- brm(
  formula = fit_gp_brms,
  data = dat_scholer2020,
  data2 = list(vcv = vcv),
  prior = prior,
  iter = 2000, 
  warmup = 1000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
summary(sp_exp_quad_brms)
