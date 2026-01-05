pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, sf, openxlsx, purrr, stringr, readr, lubridate, magrittr, janitor, glmmTMB)

# example 1 Coetzee_2014 (do not use for MS) ----
dat_coetzee <- read.xlsx(here("data", "examples", "Coetzee_2014.xlsx"), sheet = 2)
dat_coetzee$const <- 1 # add constant for spatial models

# project to a planar coordinate system and convert to kilometres
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

prof_tau2 <- profile(metafor_eg1, cline = TRUE)

pr_tau <- profile(metafor_eg1, tau2 = 1, xlim = c(1e-5, 5e-3), steps = 100)
plot(pr_tau, main = expression("Profile log-likelihood for " * tau^2))
confint(pr_tau)  

pr_rho <- profile(metafor_eg1, rho = 1, xlim = c(0, 8e9), steps = 100)
plot(pr_rho, main = expression("Profile log-likelihood for " * rho))
confint(pr_rho)

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


# example 2 Romano_2023 (do not use for MS as the brms cannot converge )----
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

exp_eg2 <- readRDS(here("Rdata", "EXP_eg2_metafor.rds"))
confint(exp_eg2)

# add moderator
dat_romano$Nest_type <- factor(dat_romano$Nest_type)

system.time(EXP_eg2_1 <- rma.mv(yi, vi, 
                              mods = ~ Nest_type - 1,
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

saveRDS(EXP_eg2_1, here("Rdata", "EXP_eg2.1_metafor.rds"))

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

tmb_2$fit$par[[4]]
tmb_2 <- readRDS(here("Rdata", "EXP_eg2_glmmTMB.rds"))

## brms ---- 
# 1) Ensure the session can handle UTF-8
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

# 2) Convert all character columns to UTF-8
char_to_utf8 <- function(x) if (is.character(x)) iconv(x, from = "", to = "UTF-8") else x
dat_romano[] <- lapply(dat_romano, char_to_utf8)

# 3) Create CLEAN factor IDs for every grouping factor used in your formula
#    (replace names here with the ones you actually use, e.g., Study, crop_species, effect_id, etc.)
library(stringi)
clean_fac <- function(x) {
  x <- as.character(x)
  x <- stri_trans_general(x, "Latin-ASCII")             # drop diacritics
  x <- gsub("[^A-Za-z0-9]+", "_", x)                    # non-alnum -> underscore
  x <- gsub("^_+|_+$", "", x)                           # trim
  factor(make.unique(x, sep = "_"))                     # ensure uniqueness
}

dat_romano$study_id     <- clean_fac(dat_romano$Study)
dat_romano$species_phylogeny <- clean_fac(dat_romano$Species_phylogeny)
dat_romano$effect_id    <- clean_fac(dat_romano$effect_id)  # if effect_id is character
# If effect_id was numeric already, just do: dat_romano$effect_id <- factor(dat_romano$effect_id)

# Optional: keep a key so you can map back later
key_study <- unique(data.frame(study_id = dat_romano$study_id, Study = dat_romano$Study))

max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)
options(mc.cores = num_chains) 

fit_eg2 <- bf(yi | se(sqrt(vi)) ~ 1 + 
                (1 | study_id) + 
                (1 | effect_id) + 
                gp(x_km, y_km, 
                   cov = "exponential",
                   scale = FALSE))

prior <- get_prior(formula = fit_eg2,
                   data = dat_romano,
                   family = gaussian()
)

m_exp_eg2 <- brm(formula = fit_eg2,
             data = dat_romano,
             family = gaussian(),
             prior = prior,
             iter = 10000,
             warmup = 7000,
             thin = 1,
             chains = num_chains,
             backend = "cmdstanr",
             threads = threading(threads_per_chain), 
             control = list(adapt_delta = 0.95, max_treedepth = 12)
)

summary(m_exp_eg2)
saveRDS(m_exp_eg2, here("Rdata", "EXP_eg2_brms.rds"))

## summary of results ----
metafor2 <- data.frame(
  model  = "metafor",
  logLik = as.numeric(logLik(exp_eg2)),  
  est = as.numeric(exp_eg2$b[1]),
  se = as.numeric(exp_eg2$se[1]),
  sigma2.u = as.numeric(exp_eg2$sigma2[2]),
  sigma2.m = as.numeric(exp_eg2$sigma2[1]),
  tau2 = as.numeric(exp_eg2$tau2[1]),   
  rho = as.numeric(exp_eg2$rho[1])     
  )

tmb2 <- data.frame(
  model   = "glmmTMB",
  logLik  = as.numeric(logLik(tmb_2)[1]),
  est = unlist(fixef(tmb_2))[[1]],
  se  = sqrt(vcov(tmb_2)$cond[1,1]),                
  sigma2.u = as.numeric(tmb_2_varcor$Study)[1],
  sigma2.m = sigma(tmb_2)^2,
  tau2 = as.numeric(tmb_2_varcor$const.1[1]),
  rho = as.numeric(exp(tmb_2$fit$par[[4]]))
  )

output <- rbind(metafor2, tmb2)
knitr::kable(output)


tmb_2_1 <- glmmTMB(
  yi ~ Nest_type - 1 +
    equalto(0 + effect_id | const, VCV) + 
    (1 | Study) +
    exp(pos + 0 | const),     
  data = dat_romano,
  REML = TRUE
)



## check the dataset for simulated dataset
#  distribution of the other dataset and ideas for the parameter values

summary(dat_romano)
ggplot(dat_romano, aes(x = x_km)) + 
  geom_histogram(binwidth = 500) + theme_classic()

ggplot(dat_romano, aes(x = y_km)) + 
  geom_histogram(binwidth = 500) + theme_classic()

# example 3 Maire_2019 ----
# In the original dataset, spatial correlations among sampling stations were modelled using station-level coordinates along rivers, but precise latitude-longitude data for the individual stations were not available - and only the study sites could be georeferenced. So we used site-level coordinates as proxies. While this approximation allowed us to compare the implementation of spatial meta-analysis models across brms, glmmTMB, and metafor, the results are not directly comparable to those reported in Maire et al. (2019), where station-level spatial random effects were used…

dat <- dat.maire2019$dat
dmat <- dat.maire2019$dmat
dmat[1:35,1:35]

sites <- tribble(
  ~site_id, ~site_name,   ~lat,    ~lon,
  "site1",  "Belleville",  47.5,    2.65,
  "site2",  "Bugey",       45.8,    5.27,
  "site3",  "Chinon",      47.2,    0.18,
  "site4",  "Chooz",       50.1,    4.80,
  "site5",  "Civaux",      46.4,    0.65,
  "site6",  "Cruas",       44.6,    4.76,
  "site7",  "Dampierre",   47.7,    2.51,
  "site8",  "Nogent",      48.5,    3.52,
  "site9",  "St Alban",    45.3,    4.75,
  "site10", "St Laurent",  47.7,    2.53,
  "site11", "Tricastin",   44.3,    4.73
)


dat_maire2019 <- dat %>%
  left_join(sites, by = c("site" = "site_id")) |> 
  mutate(id = row_number())

dat_maire2019$const <- 1 # add constant for spatial models

# project to a planar coordinate system and convert to kilometers
dat_sf <- st_as_sf(dat_maire2019, coords = c("lon", "lat"), crs = 4326)
dat_sf_proj <- st_transform(dat_sf, crs = 3857)

coords_m <- st_coordinates(dat_sf_proj)
dat_maire2019$x_km <- coords_m[,1] / 1000
dat_maire2019$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_maire2019$x_km, dat_maire2019$y_km)
dist_matrix_euclid <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid) <- dat_maire2019$id
colnames(dist_matrix_euclid) <- dat_maire2019$id

head(dat_maire2019)

## metafor ----
metafor_eg3_exp <- rma.mv(s1, vars1, 
                          random = list(
                            ~ 1|id,
                            ~ 1|site,
                            ~ x_km + y_km | const), 
                  struct = "SPEXP", 
                  data = dat_maire2019,
                  sparse = TRUE,
                  verbose = TRUE
                  )

summary(metafor_eg3_exp)
# Multivariate Meta-Analysis Model (k = 35; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -199.0469   398.0939   408.0939   415.7257   410.2367   
# 
# Variance Components:
#   
#   estim     sqrt  nlvls  fixed  factor 
# sigma^2.1   819.5017  28.6269     35     no      id 
# sigma^2.2  1653.5229  40.6635     11     no    site 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 11)
# 
# estim     sqrt  fixed 
# tau^2      8548.0233  92.4555     no 
# rho          97.3090              no 
# 
# Test for Heterogeneity:
#   Q(df = 34) = 191.4837, p-val < .0001
# 
# Model Results:
#   
#   estimate       se    zval    pval     ci.lb     ci.ub      
# 190.8255  40.9187  4.6635  <.0001  110.6263  271.0246  ***

confint(metafor_eg3_exp)

# estimate  ci.lb     ci.ub 
# sigma^2.1 819.5017 0.0000 4240.3900 
# sigma.1    28.6269 0.0000   65.1183 
# 
# estimate  ci.lb      ci.ub 
# sigma^2.2 1653.5229 0.0000 21698.5229 
# sigma.2     40.6635 0.0000   147.3042 
# 
# estimate  ci.lb        ci.ub 
# tau^2 8548.0233 0.0000 >854802.3263 
# tau     92.4555 0.0000    >924.5552 
# 
# estimate   ci.lb     ci.ub 
# rho  97.3090 <9.7309 >973.0899 

metafor_eg3_gau <- rma.mv(s1, vars1, 
                          random = list(
                            ~ 1|id,
                            ~ 1|site,
                            ~ x_km + y_km | const), 
                          struct = "SPGAU", 
                          data = dat_maire2019,
                          sparse = TRUE,
                          verbose = TRUE
)

summary(metafor_eg3_gau)
# Multivariate Meta-Analysis Model (k = 35; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -199.1047   398.2094   408.2094   415.8412   410.3522   
# 
# Variance Components:
#   
#   estim     sqrt  nlvls  fixed  factor 
# sigma^2.1   805.9835  28.3898     35     no      id 
# sigma^2.2  3557.3247  59.6433     11     no    site 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 11)
# 
# estim     sqrt  fixed 
# tau^2      6637.4779  81.4707     no 
# rho         143.9412              no 
# 
# Test for Heterogeneity:
#   Q(df = 34) = 191.4837, p-val < .0001
# 
# Model Results:
#   
#   estimate       se    zval    pval     ci.lb     ci.ub      
# 186.2795  40.6944  4.5775  <.0001  106.5199  266.0392  *** 

confint(metafor_eg3_gau)
## glmmTMB ----

dat_maire2019$id <- factor(dat_maire2019$id)
VCV <- diag(dat_maire2019$vars1, nrow = nrow(dat_maire2019)) 
rownames(VCV)<- colnames(VCV)<- dat_maire2019$id
VCV[1:5, 1:5]

dat_maire2019$pos <- numFactor(dat_maire2019$x_km, dat_maire2019$y_km)

tmb_eg3_exp <- glmmTMB(
  s1 ~ 1 +
    equalto(0 + id | const, VCV) + 
    (1 | site) +
    exp(pos + 0 | const),     
  data = dat_maire2019,
  REML = TRUE
  )


head(confint(tmb_eg3_exp), 10)
#                                                           2.5 %   97.5 %  Estimate
# (Intercept)                                            110.019957 271.6310 190.82547
# Std.Dev.(Intercept)|site                                 3.368291 490.9096  40.66358
# Std.Dev.pos(526.541191452184,5511.98585679308)|const.1  41.433901 206.3050  92.45552 <- tau

sigma(tmb_eg3_exp) # for id
# 28.62694

tmb_3_varcor <- VarCorr(tmb_eg3_exp)$cond
exp(tmb_eg3_exp$fit$par[[4]]) # for spatial range parameter rho
# 97.30914

tmb_eg3_gau <- glmmTMB(
  s1 ~ 1 +
    equalto(0 + id | const, VCV) + 
    (1 | site) +
    gau(pos + 0 | const),     
  data = dat_maire2019,
  REML = TRUE
)
head(confint(tmb_eg3_gau), 10)
exp(tmb_eg3_gau$fit$par[[4]])

## brms ----

fit_eg3 <- bf(s1 | se(sqrt(vars1)) ~ 1 + 
              (1 | site) + 
              (1 | id) + 
              gp(x_km, y_km, 
                 cov = "exponential",
                 scale = FALSE))

prior <- get_prior(formula = fit_eg3,
                   data = dat_maire2019,
                   family = gaussian()
)

m_exp <- brm(formula = fit_eg3,
             data = dat_maire2019,
             family = gaussian(),
             prior = prior,
             iter = 15000,
             warmup = 12000,
             chain = 2, 
             thin = 1,
             control = list(adapt_delta = 0.98, max_treedepth = 15)
)

summary(m_exp)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: s1 | se(sqrt(vars1)) ~ 1 + (1 | site) + (1 | id) + gp(x_km, y_km, cov = "exponential", scale = FALSE) 
# Data: dat_maire2019 (Number of observations: 35) 
# Draws: 2 chains, each with iter = 15000; warmup = 12000; thin = 1;
# total post-warmup draws = 6000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gpx_kmy_km)      83.40     42.99     5.89   172.73 1.00     1749     1961
# lscale(gpx_kmy_km)   105.02    322.12    11.91   495.18 1.00     5619     3178
# 
# Multilevel Hyperparameters:
#   ~id (Number of levels: 35) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    29.60     17.25     2.11    67.63 1.00     1644     2343
# 
# ~site (Number of levels: 11) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    63.28     39.59     4.28   151.35 1.00     1467     3078
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept   192.58     42.74   107.67   278.37 1.00     3713     3596
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# pp_check(m_exp, ndraws = 200)

saveRDS(m_exp, here("Rdata", "EXP_eg3_brms.rds"))

fit_eg3_gau <- bf(s1 | se(sqrt(vars1)) ~ 1 + 
                (1 | site) + 
                (1 | id) + 
                gp(x_km, y_km, 
                   cov = "exp_quad",
                   scale = FALSE))

prior <- get_prior(formula = fit_eg3_gau,
                   data = dat_maire2019,
                   family = gaussian()
)

m_gau <- brm(formula = fit_eg3_gau,
             data = dat_maire2019,
             family = gaussian(),
             prior = prior,
             iter = 15000,
             warmup = 12000,
             chain = 2, 
             thin = 1,
             control = list(adapt_delta = 0.98, max_treedepth = 15)
)

summary(m_gau)
saveRDS(m_gau, here("Rdata", "GAU_eg3_brms.rds"))

## summary of results ----
metafor3 <- data.frame(model = "metafor", 
                      logLik = logLik(metafor_eg3_exp),
                      est = metafor_eg3_exp$b[[1]], 
                      se = metafor_eg3_exp$se[[1]], 
                      sigma2.u = metafor_eg3_exp$sigma2[2], ## among study variance
                      sigma2.m = metafor_eg3_exp$sigma2[1], ## within study variance
                      tau2 = metafor_eg3_exp$tau2,
                      rho = metafor_eg3_exp$rho
                      )


tmb3 <- data.frame(model = "glmmTMB",
                     logLik = logLik(tmb_eg3_exp)[1],
                     est = unlist(fixef(tmb_eg3_exp))[[1]], #overall mean
                     se = as.numeric(sqrt(vcov(tmb_eg3_exp)[[1]])), #overall mean SE
                     sigma2.u = (tmb_3_varcor$site[1]), ##among study variance estimate
                     sigma2.m = sigma(tmb_eg3_exp)^2, ##within study variance estimate
                     tau2 = (tmb_3_varcor$const.1[1]), ##spatial variance estimate
                     rho = exp(tmb_eg3_exp$fit$par[[4]]) ## rho
)

posterior_summary(m_exp)
summ <- posterior_summary(
  m_exp,
  variable = c(
    "b_Intercept",
    "sd_id__Intercept",
    "sd_site__Intercept",
    "sdgp_gpx_kmy_km",
    "lscale_gpx_kmy_km"
  )
)

brms3 <- data.frame(
  model = "brms",
  logLik = NA_real_, 
  est = summ["b_Intercept","Estimate"],
  se  = summ["b_Intercept","Est.Error"],
  sigma2.u = (summ["sd_site__Intercept","Estimate"])^2,
  sigma2.m = (summ["sd_id__Intercept","Estimate"])^2,
  tau2 = (summ["sdgp_gpx_kmy_km","Estimate"])^2,
  rho = summ["lscale_gpx_kmy_km","Estimate"]
)

output_eg3 <- rbind(metafor3, tmb3, brms3)

knitr::kable(output_eg3)
# |model   |    logLik|      est|       se| sigma2.u| sigma2.m|     tau2|       rho|
#   |:-------|---------:|--------:|--------:|--------:|--------:|--------:|---------:|
#   |metafor | -199.0469| 190.8255| 40.91869| 1653.523| 819.5017| 8548.023|  97.30899|
#   |glmmTMB | -200.8246| 190.8255| 41.22806| 1653.527| 819.5016| 8548.023|  97.30914|
#   |brms    |        NA| 192.5790| 42.74387| 4004.832| 876.2053| 6955.341| 105.02247|


# example 4 Roger et al. 2024 ----

dat_Roger <- read.csv(here("data", "examples", "Roger_etal_2024", "Roger_etal_2024.csv"))
head(dat_Roger)


dat_Roger$const <- 1 # add constant for spatial models
dat_Roger$effect_id <- seq_len(nrow(dat_Roger))

dat_Roger <- dat_Roger %>%
  filter(!is.na(longitude))
coords2 <- cbind(dat_Roger$longitude, dat_Roger$latitude)


# project to a planar coordinate system and convert to kilometers
dat_sf2 <- st_as_sf(dat_Roger, coords = c("longitude", "latitude"), crs = 4326)
dat_sf2_proj <- st_transform(dat_sf2, crs = 3857)

coords_m <- st_coordinates(dat_sf2_proj)
dat_Roger$x_km <- coords_m[,1] / 1000
dat_Roger$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_Roger$x_km, dat_Roger$y_km)
dist_matrix_euclid2 <- as.matrix(dist(coords_km))

# row and column names
rownames(dist_matrix_euclid2) <- dat_Roger$effect_id
colnames(dist_matrix_euclid2) <- dat_Roger$effect_id

dat_Roger <- dat_Roger %>%
  group_by(latitude, longitude) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup()
names(dat_Roger)


## metafor ----
system.time(EXP_eg4 <- rma.mv(d_Hedges, var_Hedges, 
                              random = list(
                                ~ 1|effect_id,
                                ~ 1|site_id,
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
#     user   system  elapsed 
#   969.641   24.550 1018.565 

summary(EXP_eg4)
# logLik    Deviance         AIC         BIC        AICc   
# -4031.5291   8063.0583   8073.0583   8101.8903   8073.0837   

# Variance Components:

#             estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  0.7940  0.8911   2361     no  effect_id 
# sigma^2.2  0.9288  0.9637    383     no    site_id 

# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 383)

#               estim    sqrt  fixed 
# tau^2        0.2989  0.5467     no 
# rho        188.1869             no 

# Test for Heterogeneity:
# Q(df = 2360) = 21273.5149, p-val < .0001

# Model Results:

# estimate      se     zval    pval    ci.lb    ci.ub      
#  -0.3347  0.0765  -4.3727  <.0001  -0.4846  -0.1847  *** 

saveRDS(EXP_eg4, here("Rdata","EXP_eg4_metafor_site_id"))

EXP_eg4_metafor <- readRDS(here("Rdata", "EXP_eg4_metafor_site_id.rds"))

confint(EXP_eg4)
# estimate  ci.lb  ci.ub 
# sigma^2.1   0.7940 0.7262 0.8684 
# sigma.1     0.8911 0.8522 0.9319 
# 
# estimate  ci.lb  ci.ub 
# sigma^2.2   0.9288 0.6287 1.2893 
# sigma.2     0.9637 0.7929 1.1355 
# 
# estimate ci.lb ci.ub 
# tau^2   0.2989    NA    NA 
# tau     0.5467    NA    NA 
# 
# estimate   ci.lb      ci.ub 
# rho 188.1869 37.0656 >1881.8688 

## glmmTMB ----
dat_Roger$effect_id <- factor(dat_Roger$effect_id)
dat_Roger$site_id <- factor(dat_Roger$site_id)

VCV <- diag(dat_Roger$var_Hedges, nrow = nrow(dat_Roger)) 
rownames(VCV)<- colnames(VCV)<- dat_Roger$effect_id
VCV[1:5, 1:5]

dat_Roger$pos <- numFactor(dat_Roger$x_km, dat_Roger$y_km)

system.time(
  tmb_4 <- glmmTMB(d_Hedges ~ 1 
                   + equalto(0+effect_id|const, VCV)
                   + (1|site_id)
                   + exp(pos+0|const),
                   data = dat_Roger, 
                   REML=TRUE)
)

head(confint(tmb_4), 10)
#                                                               2.5 %     97.5 %   Estimate
# (Intercept)                                              -0.46105870 -0.2085821 -0.3348204
# Std.Dev.(Intercept)|site_id                               0.05190986  6.8087100  0.5945075
# Std.Dev.pos(16262.6644099893,-5204.51951303405)|const.1   0.34802159  2.5227325  0.9369981
# Std.Dev.pos(16280.4755285163,-5116.14629455727)|const.1   0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7965.96086752978,-5019.4733863405)|const.1   0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7959.34359171906,-4693.06364429579)|const.1  0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7960.45678662699,-4691.63535918108)|const.1  0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7096.61753807119,-4685.92422115392)|const.1  0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7992.73943895704,-4607.717759142)|const.1    0.34802159  2.5227325  0.9369981
# Std.Dev.pos(-7903.68384632242,-4579.4258128701)|const.1   0.34802159  2.5227325  0.9369981

sigma(tmb_4) 
# 0.8907115
tmb_4_varcor <- VarCorr(tmb_4)$cond
exp(tmb_4$fit$par[[2]])
# 1.110404

## brms ----

fit_eg4 <- bf(d_Hedges | se(sqrt(var_Hedges)) ~ 1 + 
                 (1 | study_id) + 
                (1 | effect_id) + 
                gp(x_km, y_km, 
                   cov = "exponential",
                   scale = FALSE))

prior <- get_prior(formula = fit_eg4,
                   data = dat_Roger,
                   family = gaussian()
)

system.time(
  m_exp4_brms <- brm(formula = fit_eg4,
             data = dat_Roger,
             family = gaussian(),
             prior = prior,
             iter = 6000,
             warmup = 3000,
             chain = 2, 
             thin = 1,
             control = list(adapt_delta = 0.98, max_treedepth = 15)
)
)
#   user   system  elapsed 
# 23276.560   588.641 76759.311

summary(m_exp4_brms)
# Family: gaussian 
# Links: mu = identity 
# Formula: d_Hedges | se(sqrt(var_Hedges)) ~ 1 + (1 | study_id) + (1 | effect_id) + gp(x_km, y_km, cov = "exponential", scale = FALSE) 
# Data: dat_Roger (Number of observations: 2361) 
# Draws: 2 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 6000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gpx_kmy_km)       0.22      0.14     0.01     0.52 1.00      270      734
# lscale(gpx_kmy_km)  6184.77  81608.27   419.46 23837.67 1.00      506     1264
# 
# Multilevel Hyperparameters:
#   ~effect_id (Number of levels: 2361) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.88      0.02     0.84     0.92 1.00     1223     2935
# 
# ~study_id (Number of levels: 393) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.09      0.06     0.98     1.21 1.00      947     1414
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -0.33      0.11    -0.53    -0.09 1.00     1119     1057
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


saveRDS(m_exp4_brms, here("Rdata", "EXP_eg4_brms.rds"))
m_exp4_brms <- readRDS(here("Rdata", "EXP_eg4_brms.rds"))

### different way to specify vcv matrix ----
# dat_Roger$effect_id <- seq_len(nrow(dat_Roger))
dat_Roger$effect_id <- as.character(dat_Roger$effect_id)
vcv <- diag(dat_Roger$var_Hedges)

rownames(vcv) <- colnames(vcv) <- dat_Roger$effect_id

fit_eg4_v1 <- bf(d_Hedges ~ 1 +
                 (1|study_id) + # this is u_l (the between-study effect)
                 (1|gr(effect_id, cov = vcv)) + # this is m (sampling error)
                gp(x_km, y_km, 
                   cov = "exponential",
                   scale = FALSE)
                )

# generate default priors
prior <- default_prior(fit_eg4_v1, 
                           data=dat_Roger, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior$prior[5] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior

# fitting model
fit_exp4_brms_v1 <- brm(
  formula = fit_eg4_v1,
  data = dat_Roger,
  data2 = list(vcv=vcv),
  chains = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_exp4_brms_v1)


fit_eg4_1 <- bf(d_Hedges | se(sqrt(var_Hedges)) ~ response - 1 + 
                # (1 | site) + 
                (1 | effect_id) + 
                gp(x_km, y_km, 
                   cov = "exponential",
                   scale = FALSE))

prior <- get_prior(formula = fit_eg4_1,
                   data = dat_Roger,
                   family = gaussian()
)

system.time(
  m_exp4_1_brms <- brm(formula = fit_eg4_1,
                     data = dat_Roger,
                     family = gaussian(),
                     prior = prior,
                     iter = 2000,
                     warmup = 1000,
                     chain = 2, 
                     thin = 1,
                     control = list(adapt_delta = 0.98, max_treedepth = 15)
  )
)
# saveRDS(m_exp4_1_brms, here("Rdata", "EXP_eg4_1_brms.rds"))
summary(m_exp4_1_brms)
### summarise ----
metafor4 <- data.frame(model = "metafor", 
                       logLik = logLik(EXP_eg4_metafor),
                       est = EXP_eg4_metafor$b[[1]], 
                       se = EXP_eg4_metafor$se[[1]], 
                       # sigma2.u = EXP_eg4$sigma2[2], ## among study variance
                       sigma2.m = EXP_eg4_metafor$sigma2[1], ## within study variance
                       tau2 = EXP_eg4_metafor$tau2,
                       rho = EXP_eg4_metafor$rho
)



tmb4 <- data.frame(model = "glmmTMB",
                   logLik = logLik(tmb_4)[1],
                   est = unlist(fixef(tmb_4))[[1]], #overall mean
                   se = as.numeric(sqrt(vcov(tmb_4)[[1]])), #overall mean SE
                   sigma2.u = (tmb_4_varcor$site_id[1]), ## among site variance estimate
                   sigma2.m = sigma(tmb_4)^2, ## within study variance estimate -> effect id
                   tau2 = (tmb_4_varcor$const.1[1]), ##spatial variance estimate
                   rho = exp(tmb_4$fit$par[]) ## rho
)



posterior_summary(m_exp4_brms)
summ <- posterior_summary(
  m_exp4_brms,
  variable = c(
    "b_Intercept",
    "sd_effect_id__Intercept",
    "sdgp_gpx_kmy_km",
    "lscale_gpx_kmy_km"
  )
)

brms4 <- data.frame(
  model = "brms",
  logLik = NA_real_, 
  est = summ["b_Intercept","Estimate"],
  se  = summ["b_Intercept","Est.Error"],
  sigma2.m = (summ["sd_effect_id__Intercept","Estimate"])^2,
  tau2 = (summ["sdgp_gpx_kmy_km","Estimate"])^2,
  rho = summ["lscale_gpx_kmy_km","Estimate"]
)

output_eg4 <- rbind(metafor4, tmb4, brms4)

knitr::kable(output_eg4)
