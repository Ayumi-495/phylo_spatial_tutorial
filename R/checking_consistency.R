pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape, patchwork,
               dplyr, tidyr,MCMCglmm, bayesplot, tidybayes, posterior, orchaRd,purrr, stringr, readr,
               lubridate, magrittr, janitor, glmmTMB, rotl, sf)

dat_Roger <- read.csv(here("data", "examples", "Roger_etal_2024", "Roger_etal_2024.csv")) |>
  dplyr::filter(!is.na(longitude)) %>%
  dplyr::mutate(
    const = 1,
    effect_id = seq_len(dplyr::n())
  )

dat_Roger$effect_id_chr <- as.character(dat_Roger$effect_id)

lev_num <- as.character(sort(as.integer(unique(dat_Roger$effect_id_chr))))

dat_Roger$effect_id <- factor(dat_Roger$effect_id_chr, levels = lev_num)
head(levels(dat_Roger$effect_id), 20) 

coords2 <- cbind(dat_Roger$longitude, dat_Roger$latitude)


# project to a planar coordinate system and convert to kilometers
dat_sf2 <- st_as_sf(dat_Roger, coords = c("longitude", "latitude"), crs = 4326)
dat_sf2_proj <- st_transform(dat_sf2, crs = 3857)

coords_m <- st_coordinates(dat_sf2_proj)
dat_Roger$x_km <- coords_m[,1] / 1000
dat_Roger$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_Roger$x_km, dat_Roger$y_km)
dist_matrix_euclid2 <- as.matrix(dist(coords_km))

id <- as.character(dat_Roger$effect_id)
rownames(dist_matrix_euclid2) <- colnames(dist_matrix_euclid2) <- id

stopifnot(identical(rownames(dist_matrix_euclid2), id))

### metafor ----

system.time(sp_eg1_exp_mf <- rma.mv(d_Hedges, var_Hedges, 
                                    random = list(
                                      ~ 1|effect_id,
                                      ~ x_km + y_km |const
                                    ), 
                                    struct = "SPEXP", 
                                    data = dat_Roger,
                                    sparse = TRUE,
                                    verbose = TRUE,
                                    method = "REML",
                                    test = "t"
                                    )
            )

summary(sp_eg1_exp_mf)
# Multivariate Meta-Analysis Model (k = 2361; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc   
# -4034.4852   8068.9704   8076.9704   8100.0360   8076.9874   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed     factor 
# sigma^2    0.7935  0.8908   2361     no  effect_id 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 383)
# 
# estim    sqrt  fixed 
# tau^2      1.2330  1.1104     no 
# rho        0.1744             no 
# 
# Test for Heterogeneity:
#   Q(df = 2360) = 21273.5149, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     tval    df    pval    ci.lb    ci.ub      
# -0.3349  0.0644  -5.2006  2360  <.0001  -0.4612  -0.2086  *** 

confint(sp_eg1_exp_mf)
#         estimate  ci.lb  ci.ub 
# sigma^2   0.7935 0.7257 0.8678 
# sigma     0.8908 0.8519 0.9316 

#       estimate  ci.lb  ci.ub 
# tau^2   1.2330 1.0147 1.4976 
# tau     1.1104 1.0073 1.2238 

#     estimate   ci.lb  ci.ub 
# rho   0.1744 <0.0174 0.8572 
sp_eg1_exp_mf <- readRDS(here("Rdata", "tutorial_v2", "sp_eg1_exp_mf.rds"))
I2_total <- (sp_eg1_exp_mf$sigma2[1] + sp_eg1_exp_mf$tau2 ) / (sp_eg1_exp_mf$sigma2[1] + sp_eg1_exp_mf$tau2 + mean(dat_Roger$var_Hedges)) * 100
I2_effect_id <- sp_eg1_exp_mf$sigma2[1] / (sp_eg1_exp_mf$sigma2[1] + sp_eg1_exp_mf$tau2 + mean(dat_Roger$var_Hedges)) * 100
I2_spatial <- sp_eg1_exp_mf$tau2 / (sp_eg1_exp_mf$sigma2[1] + sp_eg1_exp_mf$tau2 + mean(dat_Roger$var_Hedges)) * 100
### brms ----

vcv <- diag(dat_Roger$var_Hedges)
rownames(vcv) <- colnames(vcv) <- as.character(dat_Roger$effect_id)

stopifnot(setequal(rownames(vcv), id_levels))

vcv <- vcv[id_levels, id_levels]
stopifnot(identical(rownames(vcv), id_levels))

fit_sp_eg1_exp_brms <- bf(d_Hedges ~ 1 +
                            (1|gr(effect_id, cov = vcv)) + # this is m (sampling error)
                            gp(x_km, y_km, 
                               cov = "exponential",
                               scale = FALSE)
)

# generate default priors
prior <- default_prior(fit_sp_eg1_exp_brms, 
                       data = dat_Roger, 
                       data2 = list(vcv = vcv),
                       family = gaussian())
prior$prior[5] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior
# fitting model
sp_eg1_exp_brms <- brm(
  formula = fit_sp_eg1_exp_brms,
  data = dat_Roger,
  data2 = list(vcv=vcv),
  chains = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(sp_eg1_exp_brms)
# Family: gaussian 
# Links: mu = identity 
# Formula: d_Hedges ~ 1 + (1 | gr(effect_id, cov = vcv)) + gp(x_km, y_km, cov = "exponential", scale = FALSE) 
# Data: dat_Roger (Number of observations: 2361) 
# Draws: 2 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 6000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gpx_kmy_km)       1.14      0.07     1.00     1.29 1.00     1720     3259
# lscale(gpx_kmy_km)   127.84     26.58    85.01   186.52 1.00     2336     3293
# 
# Multilevel Hyperparameters:
#   ~effect_id (Number of levels: 2361) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.00      0.00     1.00     1.00   NA       NA       NA
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -0.32      0.09    -0.49    -0.14 1.00     1691     2849
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.97      0.02     0.92     1.01 1.00     2786     4157

saveRDS(sp_eg1_exp_mf, here("Rdata", "tutorial_v2", "sp_eg1_exp_brms.rds"))


### glmmTMB ----
id_levels <- levels(dat_Roger$effect_id)
dat_Roger$pos <- numFactor(dat_Roger$x_km, dat_Roger$y_km)
VCV <- diag(dat_Roger$var_Hedges) 
rownames(VCV) <- colnames(VCV) <- as.character(dat_Roger$effect_id)

stopifnot(setequal(rownames(VCV), id_levels))

VCV <- VCV[id_levels, id_levels]

stopifnot(identical(rownames(VCV), id_levels))
stopifnot(identical(colnames(VCV), id_levels))

system.time(
  sp_eg1_exp_tmb <- glmmTMB(d_Hedges ~ 1 
                            + equalto(0+effect_id|const, VCV)
                            # + (1|site_id)
                            + exp(pos+0|const),
                            data = dat_Roger, 
                            REML=TRUE)
)

summary(sp_eg1_exp_tmb)　

#  Family: gaussian  ( identity )
# Formula:          
# d_Hedges ~ 1 + equalto(0 + effect_id | const, VCV) + exp(pos +      0 | const)
# Data: dat_Roger

#       AIC       BIC    logLik -2*log(L)  df.resid 
#    8084.7    8107.8   -4038.4    8076.7      2357 

# Random effects:

# Conditional model:
#  Groups   Name                                     Variance  Std.Dev. Corr     
#  const    effect_id1                                0.202208 0.44968           
#           effect_id2                                0.166676 0.40826  0.00     
#           effect_id3                                0.225770 0.47515  0.00 0.00
#           effect_id4                                0.197994 0.44497  0.00 0.00
#           effect_id5                                3.676968 1.91754  0.00 0.00
#           effect_id6                                0.315616 0.56180  0.00 0.00
#           effect_id7                                0.305460 0.55268  0.00 0.00
#           effect_id8                                0.313426 0.55984  0.00 0.00
#           effect_id9                                0.413540 0.64307  0.00 0.00
#           effect_id10                               0.325726 0.57072  0.00 0.00
#           ...                                      ...       ...      ...  ... 
#  const.1  pos(16262.6644099893,-5204.51951303405)   1.232997 1.11040           
#           pos(16280.4755285163,-5116.14629455727)   1.232997 1.11040  0.00     
#           pos(-7965.96086752978,-5019.4733863405)   1.232997 1.11040  0.00 0.00
#           pos(-7959.34359171906,-4693.06364429579)  1.232997 1.11040  0.00 0.00
#           pos(-7960.45678662699,-4691.63535918108)  1.232997 1.11040  0.00 0.00
#           pos(-7096.61753807119,-4685.92422115392)  1.232997 1.11040  0.00 0.00
#           pos(-7992.73943895704,-4607.717759142)    1.232997 1.11040  0.00 0.00
#           pos(-7903.68384632242,-4579.4258128701)   1.232997 1.11040  0.00 0.00
#           pos(16237.0609271069,-4578.0132445546)    1.232997 1.11040  0.00 0.00
#           pos(16304.594788628,-4572.36489569299)    1.232997 1.11040  0.00 0.00
#           ...                                      ...       ...      ...  ... 
#  Residual                                           0.793492 0.89078           
                                       
                                       
                                       
                                       
#  0.00                                  
#  0.00 0.00                             
#  0.00 0.00 0.00                        
#  0.00 0.00 0.00 0.00                   
#  0.00 0.00 0.00 0.00 0.00              
#  0.00 0.00 0.00 0.00 0.00 0.00         
#  0.00 0.00 0.00 0.00 0.00 0.00 0.00    
#  ...  ...  ...  ...  ...  ...  ...  ...
                                       
                                       
                                       
#  0.00                                  
#  0.00 0.00                             
#  0.00 0.00 0.00                        
#  0.00 0.00 0.00 0.00                   
#  0.00 0.00 0.00 0.00 0.00              
#  0.00 0.00 0.00 0.00 0.00 0.00         
#  0.00 0.00 0.00 0.00 0.00 0.00 0.00    
#  ...  ...  ...  ...  ...  ...  ...  ...
                                       
# Number of obs: 2361, groups:  const, 1

# Dispersion estimate for gaussian family (sigma^2): 0.793 

# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.33493    0.06444  -5.197 2.02e-07 ***

 head(confint(sp_eg1_exp_tmb), 10)
# 2.5 %     97.5 %   Estimate
# (Intercept)                                              -0.461234 -0.2086251 -0.3349296
# Std.Dev.pos(16262.6644099893,-5204.51951303405)|const.1   1.007168  1.2242221  1.1104039
# Std.Dev.pos(16280.4755285163,-5116.14629455727)|const.1   1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7965.96086752978,-5019.4733863405)|const.1   1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7959.34359171906,-4693.06364429579)|const.1  1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7960.45678662699,-4691.63535918108)|const.1  1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7096.61753807119,-4685.92422115392)|const.1  1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7992.73943895704,-4607.717759142)|const.1    1.007168  1.2242221  1.1104039
# Std.Dev.pos(-7903.68384632242,-4579.4258128701)|const.1   1.007168  1.2242221  1.1104039
# Std.Dev.pos(16237.0609271069,-4578.0132445546)|const.1    1.007168  1.2242221  1.1104039

sigma(sp_eg1_exp_tmb) 
# [1] 0.8907814
sp_eg1_exp_tmb_varcor <- VarCorr(sp_eg1_exp_tmb)$cond
exp(sp_eg1_exp_tmb$fit$par)
# betadisp     theta     theta 
# 0.8907814 1.1104039 0.1744060 

saveRDS(sp_eg1_exp_tmb, here("Rdata", "tutorial_v2", "sp_eg1_exp_tmb.rds"))

diagnose(sp_eg1_exp_tmb)
# Unusually large Z-statistics (|x|>5):
#   
#   disp~(Intercept) 
# -5.070227 
# 
# Large Z-statistics (estimate/std err) suggest a *possible* failure of the Wald approximation - often also
# associated with parameters that are at or near the edge of their range (e.g. random-effects standard deviations  approaching 0).  
# (Alternately, they may simply represent very well-estimated parameters; intercepts of non-centered models may fall in this category.) While the Wald p-values and standard errors listed in summary() may be unreliable, profile confidence intervals (see ?confint.glmmTMB) and likelihood ratio test p-values derived by comparing models (e.g. ?drop1) are probably still OK.  (Note that the LRT is conservative when the null value is on the boundary, e.g. a variance or zero-inflation value of 0 (Self and Liang 1987; Stram and Lee 1994; Goldman and Whelan 2000); in simple cases these p-values are approximately twice as large as they should be.)

library(DHARMa)
res = simulateResiduals(sp_eg1_exp_tmb)
plot(res, rank = T)

testUniformity(res)
testDispersion(res)
testOutliers(res)

# 2) conditional: fitted values
plotResiduals(res, form = predict(sp_eg1_exp_tmb))

# 3) conditional: sampling variance
plotResiduals(res, form = dat_Roger$var_Hedges)

# 4) group-wise: study ID 
plotResiduals(res, form = dat_Roger$study_id)  

# 5) spatial autocorrelation
testSpatialAutocorrelation(res, x = dat_Roger$x_km, y = dat_Roger$y_km)