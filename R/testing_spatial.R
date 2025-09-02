pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape, sf,
               purrr, stringr, readr, lubridate, magrittr, janitor, openxlsx,
               geosphere, rotl)


# for running brms 
max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)
options(mc.cores = num_chains) 


# example 1 Coetzee_2014----
dat_coetzee <- read.xlsx(here("data", "examples", "Coetzee_2014.xlsx"), sheet = 2)
dat_coetzee$const <- 1 # add constant for spatial models

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

## linear ----
head(dat_coetzee)
### metafor ----
metafor_1 <- rma.mv(
  yi, vi,
  random = list(
    ~ 1|datapt_id,
    ~ 1|study_id,
    ~ x_km + y_km | const
  ),
  struct = "SPLIN",
  data = dat_coetzee,
  # sparse = TRUE,
  verbose = TRUE
)
# when I sparse = TRUE, I got error…
# Error: *(<dgeMatrix>, <AsIs>) is not yet implemented; ask maintainer("Matrix") to implement the missing method
# In addition: Warning message:
# In class(x) <- unique.default(c("AsIs", oldClass(x))) :
#   Setting class(x) to multiple strings ("AsIs", "lsyMatrix", ...); result will no longer be an S4 object

summary(metafor_1)
# Multivariate Meta-Analysis Model (k = 1484; method: REML)

#     logLik    Deviance         AIC         BIC        AICc   
# -2927.4735   5854.9470   5864.9470   5891.4561   5864.9876   

# Variance Components:

#             estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  1.5639  1.2506   1484     no  datapt_id 
# sigma^2.2  0.7205  0.8488    127     no   study_id 

# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 100)

#                estim    sqrt  fixed 
# tau^2         0.0000  0.0004     no 
# rho        8091.9219             no 

# Test for Heterogeneity:
# Q(df = 1483) = 10075.9124, p-val < .0001

# Model Results:

# estimate      se    zval    pval   ci.lb   ci.ub      
#   0.4386  0.0956  4.5880  <.0001  0.2513  0.6260  *** 


### glmmTMB ----
library(glmmTMB)
dat_coetzee$datapt_id <- factor(dat_coetzee$datapt_id)

D <- dist_matrix_euclid[sort(rownames(dist_matrix_euclid)), sort(rownames(dist_matrix_euclid))]
rownames(D) <- colnames(D) <- dat_coetzee$datapt_id
D[1:5, 1:5]
# ASK: D is raw distance metrix (diagnal = 0, other entries = large positive values), it does not meet the requirements foe the matrix passed to propto()? (On my understanding, it must be a cor (or cov) matrix? i.e diagonal = 1 and positive definite matrix?)

VCV <- diag(dat_coetzee$vi, nrow = nrow(dat_coetzee)) 
rownames(VCV)<- colnames(VCV)<- dat_coetzee$datapt_id
VCV[1:5, 1:5]
# equalto(): the component that fixes each data point's known sampling variance
# propto(): the component that estimates spatial autocorrelation as "sigma^2 * correlation matrix" derived geographic distance

# okay - let's make D correlation matrix… using exponential kernel
alpha <- 0.005
R <- exp(-alpha * D)
R[1:50, 1:50]
diag(R) <- 1

tmb_1 <- glmmTMB(
  yi ~ 1 +
    equalto(0 + datapt_id | const, VCV) + 
    (1 | study_id) +
    propto(0 + datapt_id | const, R),     
  data = dat_coetzee,
  REML = TRUE
)
# Error in chol.default(cr) : the leading minor of order 5 is not positive

## exponential ----
metafor_2 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|datapt_id,
                   ~ 1|study_id, 
        ~ x_km + y_km | const
                 ), 
                 struct = "SPEXP", 
                 data = dat_coetzee,
                 sparse = TRUE,
                 verbose = TRUE
                 )
   
summary(metafor_2)

# Multivariate Meta-Analysis Model (k = 1484; method: REML)

#     logLik    Deviance         AIC         BIC        AICc   
# -2927.4735   5854.9470   5864.9470   5891.4561   5864.9876   

# Variance Components:

#             estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  1.5639  1.2506   1484     no  datapt_id 
# sigma^2.2  0.7205  0.8488    127     no   study_id 

# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 100)

#                    estim    sqrt  fixed 
# tau^2             0.0005  0.0227     no 
# rho        22091692.2263             no 

# Test for Heterogeneity:
# Q(df = 1483) = 10075.9124, p-val < .0001

# Model Results:

# estimate      se    zval    pval   ci.lb   ci.ub      
#   0.4386  0.0983  4.4639  <.0001  0.2460  0.6312  *** 
a
### brms ----
fit_1 <- bf(
  yi | se(sqrt(vi)) ~ 1 + 
    (1 | datapt_id) + 
    (1 | study_id) + 
    gp(y_km, x_km, cov = "exponential", scale = FALSE)
  )

# prior <- c(
#     prior(normal(0, 5), class = "Intercept"),
#     prior(student_t(3, 0, 2.5), class = "sd", group = "datapt_id"),
#     prior(student_t(3, 0, 2.5), class = "sd", group = "study_id"),
#     prior(lognormal(9, 1), class = "lscale"),   
#     prior(student_t(3, 0, 0.5), class = "sdgp")
#   )

brms_1 <- brm(formula = fit_1,
             data = america_dat_coetzee,
             family = gaussian(),
            # prior = prior,
             iter = 11000,
             warmup = 9000,
             chains = num_chains,
             backend = "cmdstanr",
             threads = threading(threads_per_chain),
             control = list(adapt_delta = 0.98, max_treedepth = 15)
)

summary(brms_1)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: yi | se(sqrt(vi)) ~ 1 + (1 | datapt_id) + (1 | study_id) + gp(y_km, x_km, cov = "exponential", scale = FALSE) 
#    Data: america_dat_coetzee (Number of observations: 490) 
#   Draws: 2 chains, each with iter = 11000; warmup = 9000; thin = 1;
#          total post-warmup draws = 4000

# Gaussian Process Hyperparameters:
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gpy_kmx_km)       0.59      0.31     0.05     1.30 1.00      340      578
# lscale(gpy_kmx_km)  1338.81   3279.22   154.76  6950.50 1.00     1259     1153

# Multilevel Hyperparameters:
# ~datapt_id (Number of levels: 490) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.75      0.07     1.62     1.90 1.00      734     1430

# ~study_id (Number of levels: 30) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.43      0.26     0.02     0.98 1.01      216      812

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.26      0.33    -0.33     0.96 1.00      764      940

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Gaussian process model with square-exponential covariance function

## squared exponential ----
### metafor ----
# metafor_3 <-  rma.mv(yi, vi, 
#                  random = list(
#                    ~ 1|datapt_id,
#                    ~ 1|study_id, 
#                    ~ datapt_id|const
#                  ), 
#                  struct = "SPGAU", 
#                  data = dat_coetzee,
#                  dist = list(datapt_id = dist_matrix_euclid)
# )

# summary(metafor_3)

metafor_3 <- rma.mv(
  yi, vi,
  random = list(
    ~ 1|datapt_id,
    ~ 1|study_id,
    ~ x_km + y_km | const
  ),
  struct = "SPGAU",
  data = dat_coetzee
)

summary(metafor_3)

### brms ----
fit_2 <- bf(
  yi | se(sqrt(vi)) ~ 1 + 
    (1 | datapt_id) + 
    (1 | study_id) + 
    gp(y_km, x_km, cov = "exp_quad", scale = FALSE)
)

brms_2 <- brm(formula = fit_2,
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

summary(brms_2)

## Northern and South America ----

table(dat_coetzee$continent)
america_dat_coetzee  <-  dat_coetzee[dat_coetzee$continent %in% c("NAM", "SA"), ]
table(america_dat_coetzee$continent)

# project to a planar coordinate system and convert to kilometers
dat_sf <- st_as_sf(america_dat_coetzee, coords = c("long", "lat"), crs = 4326)
dat_sf_proj <- st_transform(dat_sf, crs = 3857)

coords_m <- st_coordinates(dat_sf_proj)
america_dat_coetzee$x_km <- coords_m[,1] / 1000
america_dat_coetzee$y_km <- coords_m[,2] / 1000


coords_km <- cbind(america_dat_coetzee$x_km, america_dat_coetzee$y_km)
dist_matrix_euclid <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid) <- america_dat_coetzee$datapt_id
colnames(dist_matrix_euclid) <- america_dat_coetzee$datapt_id

### exponential ----
#### metafor ----
# EXP_m1 <- rma.mv(yi, vi, 
#                  random = list(
#                    ~ 1|datapt_id,
#                    ~ 1|study_id, 
#                    ~ datapt_id|const
#                  ), 
#                  struct = "SPEXP", 
#                  data = america_dat_coetzee,
#                  dist = list(datapt_id = dist_matrix_euclid),
#                  sparse = TRUE,
#                  verbose = TRUE
#                  )
   
# summary(EXP_m1)

metafor_US1 <- rma.mv(
  yi, vi,
  random = list(
    ~ 1|datapt_id,
    ~ 1|study_id,
    ~ x_km + y_km | const
  ),
  struct = "SPEXP",
  data = america_dat_coetzee,
  sparse = TRUE,
  verbose = TRUE
)

summary(metafor_US1)
# Multivariate Meta-Analysis Model (k = 490; method: REML)

#     logLik    Deviance         AIC         BIC        AICc   
# -1060.6555   2121.3110   2131.3110   2152.2728   2131.4352   

# Variance Components:

#             estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  3.0633  1.7502    490     no  datapt_id 
# sigma^2.2  0.0000  0.0003     30     no   study_id 

# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 23)

#               estim    sqrt  fixed 
# tau^2        0.4331  0.6581     no 
# rho        627.8325             no 

# Test for Heterogeneity:
# Q(df = 489) = 4848.5903, p-val < .0001

# Model Results:

# estimate      se    zval    pval    ci.lb   ci.ub    
#   0.1943  0.2303  0.8438  0.3988  -0.2570  0.6457    

profile(metafor_US1, rho=1, xlim=c(0,800), steps=100)
confint(metafor_US1)
#           estimate  ci.lb  ci.ub 
# sigma^2.1   3.0633 2.6138 3.6017 
# sigma.1     1.7502 1.6167 1.8978 

#           estimate  ci.lb  ci.ub 
# sigma^2.2   0.0000 0.0000 0.7838 
# sigma.2     0.0003 0.0000 0.8853 

#       estimate  ci.lb    ci.ub 
# tau^2   0.4331 0.0000 >43.3121 
# tau     0.6581 0.0000  >6.5812 

#     estimate    ci.lb      ci.ub 
# rho 627.8325 <62.7833 >6278.3253 

#### glmmTMB ----
library(glmmTMB)
head(america_dat_coetzee, 10)
america_dat_coetzee$datapt_id <- factor(america_dat_coetzee$datapt_id)

D <- dist_matrix_euclid[sort(rownames(dist_matrix_euclid)), sort(rownames(dist_matrix_euclid))]
rownames(D) <- colnames(D) <- america_dat_coetzee$datapt_id
D[1:5, 1:5]
# ASK: D is raw distance metrix (diagnal = 0, other entries = large positive values), it does not meet the requirements foe the matrix passed to propto()? (On my understanding, it must be a cor (or cov) matrix? i.e diagonal = 1 and positive definite matrix?)

VCV <- diag(america_dat_coetzee$vi, nrow = nrow(america_dat_coetzee)) 
rownames(VCV)<- colnames(VCV)<- america_dat_coetzee$datapt_id
VCV[1:5, 1:5]
# equalto(): the component that fixes each data point's known sampling variance

america_dat_coetzee$pos <- numFactor(america_dat_coetzee$x_km, america_dat_coetzee$y_km)

tmb_1 <- glmmTMB(
  yi ~ 1 +
    equalto(0 + datapt_id | const, VCV) + 
    (1 | study_id) +
    exp(pos + 0 | const),     
  data = america_dat_coetzee,
  REML = TRUE
)

head(confint(tmb_1), 10)

#### brms ----
fit_3 <- bf(
  yi | se(sqrt(vi)) ~ 1 + 
    (1 | datapt_id) + 
    (1 | study_id) + 
    gp(y_km, x_km, cov = "exp_quad", scale = FALSE)
)

brms_US1 <- brm(formula = fit_3,
             data = america_dat_coetzee,
             family = gaussian(),
             # prior = prior,
             iter = 11000,
             warmup = 9000,
             chains = num_chains,
             backend = "cmdstanr",
             threads = threading(threads_per_chain),
             control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(brms_US1)


#######################################
#  example 2 ----

dat_romano <- read.csv(here("data", "examples", "Romano_2023.csv"))

dat_romano$yi <- dat_romano$Slope
dat_romano$vi <- 1 / dat_romano$N_years # the author said they used 1/N_years for the variance

dat_romano$const <- 1 # add constant for spatial models
dat_romano$effect_id <- seq_len(nrow(dat_romano))

names(dat_romano)
summary(dat_romano)
coords <- cbind(dat_romano$Longitude, dat_romano$Latitude)
dat_romano <- dat_romano %>%
  filter(!is.na(Longitude))

# project to a planar coordinate system and convert to kilometers
dat_sf <- st_as_sf(dat_romano, coords = c("Longitude", "Latitude"), crs = 4326)
dat_sf_proj <- st_transform(dat_sf, crs = 3857)

coords_m <- st_coordinates(dat_sf_proj)
dat_romano$x_km <- coords_m[,1] / 1000
dat_romano$y_km <- coords_m[,2] / 1000


coords_km <- cbind(dat_romano$x_km, dat_romano$y_km)
dist_matrix_euclid <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid) <- dat_romano$effect_id
colnames(dist_matrix_euclid) <- dat_romano$effect_id

EXP_eg2_1 <- rma.mv(yi, vi, 
                 random = list(
                   ~ 1|effect_id,
                   ~ 1|Study, 
                   ~ effect_id|const
                 ), 
                 struct = "SPEXP", 
                 data = dat_romano,
                 dist = list(effect_id = dist_matrix_euclid),
                 sparse = TRUE,
                 verbose = TRUE
)

EXP_eg2_2 <- rma.mv(yi, vi, 
                    random = list(
                      ~ 1|effect_id,
                      ~ 1|Study, 
                      ~ x_km + y_km |const
                    ), 
                    struct = "SPEXP", 
                    data = dat_romano,
                    sparse = TRUE,
                    verbose = TRUE
)



summary(EXP_eg2_2)

# Multivariate Meta-Analysis Model (k = 5585; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc   
# -3369.3302   6738.6604   6748.6604   6781.7987   6748.6711   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  0.1242  0.3524   5585     no  effect_id 
# sigma^2.2  0.0882  0.2970    405     no      Study 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 648)
# 
# estim    sqrt  fixed 
# tau^2        0.0551  0.2347     no 
# rho        266.8900             no 
# 
# Test for Heterogeneity:
#   Q(df = 5584) = 30373.8064, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval    ci.lb    ci.ub      
# -0.1573  0.0304  -5.1815  <.0001  -0.2168  -0.0978  *** 
orchaRd::i2_ml(EXP_eg2_2)

fit_brms <- bf(
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