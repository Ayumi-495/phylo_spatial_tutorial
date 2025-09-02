pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, sf, openxlsx, purrr, stringr, readr, lubridate, magrittr, janitor, glmmTMB)

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

## metafor ----
system.time(
  metafor_eg1 <- rma.mv(
  yi, vi,
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
)
#    user  system elapsed 
# 357.534  14.692 376.312 
summary(metafor_eg1)

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

# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## glmmTMB ----
dat_coetzee$datapt_id <- factor(dat_coetzee$datapt_id)
VCV <- diag(dat_coetzee$vi, nrow = nrow(dat_coetzee)) 
rownames(VCV)<- colnames(VCV)<- dat_coetzee$datapt_id
VCV[1:5, 1:5]

dat_coetzee$pos <- numFactor(dat_coetzee$x_km, dat_coetzee$y_km)

system.time(
  tmb_1 <- glmmTMB(
    yi ~ 1 +
      equalto(0 + datapt_id | const, VCV) + 
      (1 | study_id) +
      exp(pos + 0 | const),     
    data = dat_coetzee,
    REML = TRUE
  ))
#    user  system elapsed 
#  12.250   0.394  12.665 

head(confint(tmb_1), 10)
sigma(tmb_1)^2
summary(tmb_1)
tmb_1_varcor <- VarCorr(tmb_1)$cond
# head(tmb_1_varcor, 10)

maire2019 <- dat.maire2019$dat
dmat <- dat.maire2019$dmat
dmat[1:5,1:5]


# example 2 Romano_2023----
dat_romano <- read.csv(here("data", "examples", "Romano_2023.csv"))

dat_romano$yi <- dat_romano$Slope
dat_romano$vi <- 1 / dat_romano$N_years # the author said they used 1/N_years for the variance

dat_romano$const <- 1 # add constant for spatial models
dat_romano$effect_id <- seq_len(nrow(dat_romano))

names(dat_romano)
summary(dat_romano)
coords2 <- cbind(dat_romano$Longitude, dat_romano$Latitude)
dat_romano <- dat_romano %>%
  filter(!is.na(Longitude))

# project to a planar coordinate system and convert to kilometers
dat_sf2 <- st_as_sf(dat_romano, coords = c("Longitude", "Latitude"), crs = 4326)
dat_sf2_proj <- st_transform(dat_sf2, crs = 3857)

coords_m <- st_coordinates(dat_sf2_proj)
dat_romano$x_km <- coords_m[,1] / 1000
dat_romano$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_romano$x_km, dat_romano$y_km)
dist_matrix_euclid2 <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid2) <- dat_romano$effect_id
colnames(dist_matrix_euclid2) <- dat_romano$effect_id

## metafor ----
system.time(EXP_eg2 <- rma.mv(yi, vi, 
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
)
#      user    system   elapsed 
#  8346.963   217.651 10501.001
summary(EXP_eg2)
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
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## glmmTMB ----
dat_romano$effect_id <- factor(dat_romano$effect_id)
VCV <- diag(dat_romano$vi, nrow = nrow(dat_romano)) 
rownames(VCV)<- colnames(VCV)<- dat_romano$effect_id
VCV[1:5, 1:5]

dat_romano$pos <- numFactor(dat_romano$x_km, dat_romano$y_km)

tmb_2 <- glmmTMB(
  yi ~ 1 +
    equalto(0 + effect_id | const, VCV) + 
    (1 | Study) +
    exp(pos + 0 | const),     
  data = dat_romano,
  REML = TRUE
)


head(confint(tmb_2), 10)
# 2.5 %      97.5 %
#   (Intercept)                                              -0.2178940 -0.09667794
# Std.Dev.(Intercept)|Study                                 0.2480820  0.35549281
# Std.Dev.pos(18378.8479299695,-12600.9031463533)|const.1   0.1806699  0.30485608
# Std.Dev.pos(-7010.90153016037,-10329.1461608005)|const.1  0.1806699  0.30485608
# Std.Dev.pos(15584.7287110583,-9923.52887926705)|const.1   0.1806699  0.30485608
# Std.Dev.pos(-7135.57935984884,-9542.826396182)|const.1    0.1806699  0.30485608
# Std.Dev.pos(-6730.37641336132,-8997.99883386902)|const.1  0.1806699  0.30485608
# Std.Dev.pos(-6483.24714380025,-8882.8934725054)|const.1   0.1806699  0.30485608
# Std.Dev.pos(17693.1198666829,-7276.93514117913)|const.1   0.1806699  0.30485608
# Std.Dev.pos(4233.48023486819,-7170.15629399995)|const.1   0.1806699  0.30485608
# Estimate
# (Intercept)                                              -0.1572860
# Std.Dev.(Intercept)|Study                                 0.2969703
# Std.Dev.pos(18378.8479299695,-12600.9031463533)|const.1   0.2346877
# Std.Dev.pos(-7010.90153016037,-10329.1461608005)|const.1  0.2346877
# Std.Dev.pos(15584.7287110583,-9923.52887926705)|const.1   0.2346877
# Std.Dev.pos(-7135.57935984884,-9542.826396182)|const.1    0.2346877
# Std.Dev.pos(-6730.37641336132,-8997.99883386902)|const.1  0.2346877
# Std.Dev.pos(-6483.24714380025,-8882.8934725054)|const.1   0.2346877
# Std.Dev.pos(17693.1198666829,-7276.93514117913)|const.1   0.2346877
# Std.Dev.pos(4233.48023486819,-7170.15629399995)|const.1   0.2346877
# > 
sigma(tmb_2)^2
summary(tmb_2)
tmb_2_varcor <- VarCorr(tmb_2)$cond
# head(tmb_2_varcor, 10)

tmb_2$fit$par
