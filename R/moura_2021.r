pacman::p_load(brms, metafor, metadat, glmmTMB, tidyverse, data.table, crayon, here, ape,
               purrr, stringr, readr, lubridate, magrittr, janitor, rotl)

# data  ----
dat_moura2021 <- dat.moura2021$dat
tree <- dat.moura2021$tree

# calculate r-to-z transformed correlations and corresponding sampling variances
dat_moura2021 <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat_moura2021)

# turn tree into an ultrametric one
tree <- compute.brlen(tree)

# compute phylogenetic correlation matrix
A <- vcv(tree, corr = TRUE)

# make copy of the species.id variable
dat_moura2021$species.id.phy <- dat_moura2021$species.id

# BM model ----
## metafor ----
system.time(
    BM_metafor <- rma.mv(yi, vi,
   random = list(~ 1 | study.id, 
   ~ 1 | effect.size.id, 
   ~ 1 | species.id, 
   ~ 1 | species.id.phy),
   R = list(species.id.phy=A), 
   data = dat_moura2021,
   verbose = TRUE,
   sparse = TRUE)
)

summary(BM_metafor)
# Multivariate Meta-Analysis Model (k = 1828; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -167.6726   335.3453   345.3453   372.8974   345.3782   

# Variance Components:

#             estim    sqrt  nlvls  fixed          factor    R 
# sigma^2.1  0.0192  0.1384    457     no        study.id   no 
# sigma^2.2  0.0145  0.1202   1828     no  effect.size.id   no 
# sigma^2.3  0.0557  0.2359    341     no      species.id   no 
# sigma^2.4  0.0512  0.2263    341     no  species.id.phy  yes 

# Test for Heterogeneity:
# Q(df = 1827) = 10743.8076, p-val < .0001

# Model Results:

# estimate      se    zval    pval   ci.lb   ci.ub     
#   0.3682  0.1300  2.8311  0.0046  0.1133  0.6230  ** 

# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

confint(BM_metafor)

#           estimate  ci.lb  ci.ub 
# sigma^2.1   0.0192 0.0108 0.0325 
# sigma.1     0.1384 0.1038 0.1802 

#           estimate  ci.lb  ci.ub 
# sigma^2.2   0.0145 0.0121 0.0172 
# sigma.2     0.1202 0.1099 0.1311 

#           estimate  ci.lb  ci.ub 
# sigma^2.3   0.0557 0.0334 0.0788 
# sigma.3     0.2359 0.1827 0.2807 

#           estimate  ci.lb  ci.ub 
# sigma^2.4   0.0512 0.0179 0.1792 
# sigma.4     0.2263 0.1336 0.4233 

orchaRd::i2_ml(BM_metafor)

phylo_heritability <- BM_metafor$sigma2[4]^2 / (BM_metafor$sigma2[4]^2 + BM_metafor$sigma2[3]^2)

### add one explanatory variable ----
system.time(
    BM_metafor1 <- rma.mv(yi, vi,
    mod = ~ 1 + spatially.pooled ,
   random = list(~ 1 | study.id, 
   ~ 1 | effect.size.id, 
   ~ 1 | species.id, 
   ~ 1 | species.id.phy),
   R = list(species.id.phy = A), 
   data = dat_moura2021,
   verbose = TRUE,
   sparse = TRUE)
)
#    user  system elapsed 
# 402.725  11.616 426.134 
View(dat_moura2021)
summary(BM_metafor1)
# saveRDS(BM_metafor1, here("Rdata", "moura2021_BM_metafor1.rds"))
## glmmTMB ----
A <- A[sort(rownames(A)), sort(rownames(A))]

dat_moura2021$g <- 1 

VCV <- diag(dat_moura2021$vi, nrow = nrow(dat_moura2021)) 
rownames(VCV)<- colnames(VCV)<- dat_moura2021$effect.size.id

dat_moura2021$effect.size.id <- as.factor(dat_moura2021$effect.size.id)
head(dat_moura2021)

system.time(
    BM_tmb1 <- glmmTMB(yi ~ 1 + spatially.pooled + equalto(0 + effect.size.id|g, VCV) +
                         (1|study.id) + 
                         (1|species.id.phy) +
                         propto(0 + species.id.phy|g, A),
                       data = dat_moura2021,
                       REML = T)
)
#  user  system elapsed 
#  23.254   0.502  23.86
# saveRDS(BM_tmb1, here("Rdata", "moura2021_BM_tmb1.rds"))

head(confint(BM_tmb1), 10)
#                                                                    2.5 %
# (Intercept)                                                             0.10679928
# spatially.pooledyes                                                    -0.01118452
# Std.Dev.(Intercept)|study.id                                            0.10562538
# Std.Dev.(Intercept)|species.id.phy                                     0.19247823
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1             0.13029865
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1                0.13029865
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1              0.13029865
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                     0.13029865
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1  0.13029865
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                      0.13029865
#                                                                           97.5 %
# (Intercept)                                                            0.6204510
# spatially.pooledyes                                                    0.0983696
# Std.Dev.(Intercept)|study.id                                           0.1817868
# Std.Dev.(Intercept)|species.id.phy                                     0.2874944
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1            0.3990331
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1               0.3990331
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1             0.3990331
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                    0.3990331
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1 0.3990331
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                     0.3990331
#                                                                          Estimate
# (Intercept)                                                            0.36362512
# spatially.pooledyes                                                    0.04359254
# Std.Dev.(Intercept)|study.id                                           0.13856875
# Std.Dev.(Intercept)|species.id.phy                                     0.23523692
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1            0.22802077
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1               0.22802077
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1             0.22802077
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                    0.22802077
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1 0.22802077
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                     0.22802077

sigma(BM_tmb1)^2
# 0.01441594
# summary(BM_tmb)
BM_tmb_varcor <- VarCorr(BM_tmb1)$cond
head(BM_tmb_varcor, 10)

## brms ----

fit_brms2 <- bf(yi | se(sqrt(vi)) ~ 1 + spatially.pooled + 
                  (1 | study.id) + 
                  (1 | effect.size.id) + 
                  (1 | species.id) + 
                  (1 | gr(species.id.phy, cov = A))
)

prior_brms2 <- get_prior(
  formula = fit_brms2,
  data = dat_moura2021,
  data2 = list(A = A),
  family = gaussian()
)

max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)
options(mc.cores = num_chains) 

system.time(
    m_brms2 <- brm(
  formula = fit_brms2,
  family = gaussian(),
  data = dat_moura2021,
  data2 = list(A = A),
  prior = prior_brms2,
  iter = 10000, 
  warmup = 7000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
)
summary(m_brms2)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: yi | se(sqrt(vi)) ~ 1 + spatially.pooled + (1 | study.id) + (1 | effect.size.id) + (1 | species.id) + (1 | gr(species.id.phy, cov = A)) 
#    Data: dat_moura2021 (Number of observations: 1828) 
#   Draws: 2 chains, each with iter = 10000; warmup = 7000; thin = 1;
#          total post-warmup draws = 6000

# Multilevel Hyperparameters:
# ~effect.size.id (Number of levels: 1828) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.12      0.01     0.11     0.13 1.00     1494     2868

# ~species.id (Number of levels: 341) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.23      0.03     0.17     0.28 1.00      535      624

# ~species.id.phy (Number of levels: 341) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.27      0.09     0.15     0.49 1.00      611      467

# ~study.id (Number of levels: 457) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.14      0.02     0.10     0.19 1.00      454      794

# Regression Coefficients:
#                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept               0.36      0.16     0.04     0.68 1.00     2050     2076
# spatially.pooledyes     0.04      0.03    -0.01     0.10 1.00     1800     2714

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

## summarise ----  

# BM_metafor <- readRDS(here("Rdata", "moura2021_BM_metafor.rds"))
# m_brms1 <- readRDS(here("Rdata", "moura2021_BM_brms.rds"))
# BM_tmb <- readRDS(here("Rdata", "moura2021_BM_tmb.rds"))

# est = unlist(fixef(BM_tmb))[[1]] #overal mean
# se = as.numeric(sqrt(vcov(BM_tmb)[[1]])) #overall mean SE
# sigma2.u = BM_tmb_varcor$study.id[1] ##among study variance estimate
# sigma2.non.phylo = BM_tmb_varcor$species.id.phy[1] ##non-phylo variance estimate
# sigma2.phylo = BM_tmb_varcor$g.1[1] ##phylo variance component
# sigma2.m = sigma(BM_tmb) ##within study variance estimate
# sigma2.total = sum(sigma2.u+sigma2.non.phylo+sigma2.phylo+sigma2.m)

metafor_1 <- data.frame(model = "BM_metafor", 
                      logLik = logLik(BM_metafor),
                      est = BM_metafor$b[[1]], 
                      se = BM_metafor$se[[1]], 
                      sigma2.u = sqrt(BM_metafor$sigma2[1]), ## among study variance
                      sigma2.m = sqrt(BM_metafor$sigma2[2]), ## within study variance
                      sigma2.non.phylo = sqrt(BM_metafor$sigma2[3]),
                      sigma2.phylo = sqrt(BM_metafor$sigma2[4])
                      )

glmmTMB_1 <- data.frame(model = "BM_tmb",
                      logLik = logLik(BM_tmb)[1],
                      est = unlist(fixef(BM_tmb))[[1]], #overall mean
                      se = as.numeric(sqrt(vcov(BM_tmb)[[1]])), #overall mean SE
                      sigma2.u = sqrt(BM_tmb_varcor$study.id[1]), ##among study variance estimate
                      sigma2.m = sigma(BM_tmb), ##within study variance estimate
                      sigma2.non.phylo = sqrt(BM_tmb_varcor$species.id[1]), ##non-phylo variance estimate
                      sigma2.phylo = sqrt(BM_tmb_varcor$g.1[1]) ##phylo variance component
                      )

summ <- posterior_summary(
  m_brms1,
  variable = c(
    "b_Intercept",
    "sd_effect.size.id__Intercept",
    "sd_study.id__Intercept",
    "sd_species.id__Intercept",
    "sd_species.id.phy__Intercept"
  )
)

brms_1 <- data.frame(
  model = "BM_brms",
  logLik = NA_real_, 
  est = summ["b_Intercept","Estimate"],
  se  = summ["b_Intercept","Est.Error"],
  sigma2.u        = summ["sd_study.id__Intercept","Estimate"],
  sigma2.m        = summ["sd_effect.size.id__Intercept","Estimate"],
  sigma2.non.phylo= summ["sd_species.id__Intercept","Estimate"],
  sigma2.phylo    = summ["sd_species.id.phy__Intercept","Estimate"]
)

output <- rbind(metafor_1, glmmTMB_1, brms_1)
knitr::kable(output)     
# |model      |    logLik|       est|        se|  sigma2.u|  sigma2.m| sigma2.non.phylo| sigma2.phylo|
# |:----------|---------:|---------:|---------:|---------:|---------:|----------------:|------------:|
# |BM_metafor | -167.6726| 0.3681659| 0.1300448| 0.1384141| 0.1202087|        0.2359274|    0.2263264|
# |BM_tmb     | -171.4281| 0.3681658| 0.1300559| 0.1384142| 0.1202088|        0.2359271|    0.2263262|
# |BM_brms    |        NA| 0.3655350| 0.1652585| 0.1416805| 0.1205363|        0.2288716|    0.2713950|

# OU model ----
## metafor ----
### obtain the rho for get alpha

dat_moura2021$const <- 1
I  <- matrix(1, nrow = length(tree$tip.label), ncol = length(tree$tip.label))
D <- I-A

D[1:5, 1:5]

names(dat_moura2021)
system.time(
    spatial_metafor <- rma.mv(yi, vi,
   random = list(~ 1 | study.id, 
   ~ 1 | effect.size.id, 
   ~ 1 | species.id, 
   ~ species.id.phy | const),
   dist = list(D), 
   struct = "SPEXP",
   control = list(rho.init = 1),
   data = dat_moura2021,
   verbose = TRUE,
   sparse = TRUE
   )
   )
#     user   system  elapsed 
# 1155.914   36.850 1214.621 

summary(spatial_metafor)
# Multivariate Meta-Analysis Model (k = 1828; method: REML)

#    logLik   Deviance       AIC        BIC       AICc   
# -160.3783   320.7567   332.7567   365.8192   332.8028   

# Variance Components:

#             estim    sqrt  nlvls  fixed          factor 
# sigma^2.1  0.0157  0.1254    457     no        study.id 
# sigma^2.2  0.0144  0.1201   1828     no  effect.size.id 
# sigma^2.3  0.0000  0.0000    341     no      species.id 

# outer factor: const           (nlvls = 1)
# inner term:   ~species.id.phy (nlvls = 341)

#             estim    sqrt  fixed 
# tau^2      0.1029  0.3208     no 
# rho        0.0182             no 

# Test for Heterogeneity:
# Q(df = 1827) = 10743.8076, p-val < .0001

# Model Results:

# estimate      se    zval    pval   ci.lb   ci.ub      
#   0.3514  0.0355  9.9090  <.0001  0.2819  0.4209  *** 

# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


rho <- spatial_metafor$rho
alpha <- 1/rho
A_OU <- exp(-alpha*D)

system.time(
    OU_metafor <- rma.mv(yi, vi,
   random = list(~ 1 | study.id, 
   ~ 1 | effect.size.id, 
   ~ 1 | species.id, 
   ~ 1 | species.id.phy),
   R = list(species.id.phy=A_OU), 
   data = dat_moura2021,
   verbose = TRUE,
   sparse = TRUE)
)
#    user   system  elapsed 
# 2756.764   53.458 2857.989 
summary(OU_metafor)
# Multivariate Meta-Analysis Model (k = 1828; method: REML)

#    logLik   Deviance        AIC        BIC       AICc   
# -160.3783   320.7567   330.7567   358.3088   330.7896   

# Variance Components:

#             estim    sqrt  nlvls  fixed          factor    R 
# sigma^2.1  0.0157  0.1254    457     no        study.id   no 
# sigma^2.2  0.0144  0.1201   1828     no  effect.size.id   no 
# sigma^2.3  0.0000  0.0000    341     no      species.id   no 
# sigma^2.4  0.1029  0.3208    341     no  species.id.phy  yes 

# Test for Heterogeneity:
# Q(df = 1827) = 10743.8076, p-val < .0001

# Model Results:

# estimate      se    zval    pval   ci.lb   ci.ub      
#   0.3514  0.0355  9.9090  <.0001  0.2819  0.4209  *** 

saveRDS(spatial_metafor, file = here("Rdata", "OU_spatial_metafor.rds"))

BM_metafor$sigma2[3] + BM_metafor$sigma2[1] + BM_metafor$sigma2[2]  + BM_metafor$sigma2[4]
OU_metafor$sigma2[3] + OU_metafor$sigma2[1] + OU_metafor$sigma2[2]  + OU_metafor$sigma2[4]
