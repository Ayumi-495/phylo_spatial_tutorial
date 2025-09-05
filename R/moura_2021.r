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
A <- vcv(tree, corr=TRUE)

# make copy of the species.id variable
dat_moura2021$species.id.phy <- dat_moura2021$species.id

# BM model ----
## metafor ----
BM_metafor <- rma.mv(yi, vi,
   random = list(~ 1 | study.id, 
   ~ 1 | effect.size.id, 
   ~ 1 | species.id, 
   ~ 1 | species.id.phy),
   R = list(species.id.phy=A), 
   data = dat_moura2021,
   verbose = TRUE,
   sparse = TRUE)

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

orchaRd::i2_ml(BM_metafor)

phylo_heritability <- BM_metafor$sigma2[4]^2 / (BM_metafor$sigma2[4]^2 + BM_metafor$sigma2[3]^2)


## glmmTMB ----
A <- A[sort(rownames(A)), sort(rownames(A))]

dat_moura2021$g <- 1 

VCV <- diag(dat_moura2021$vi, nrow = nrow(dat_moura2021)) 
rownames(VCV)<- colnames(VCV)<- dat_moura2021$effect.size.id

dat_moura2021$effect.size.id <- as.factor(dat_moura2021$effect.size.id)
head(dat_moura2021)

BM_tmb <- glmmTMB(yi ~ 1 + equalto(0 + effect.size.id|g, VCV) +
                         (1|study.id) + 
                         (1|species.id.phy) +
                         propto(0 + species.id.phy|g, A),
                       data = dat_moura2021,
                       REML = T)
# saveRDS(BM_tmb, here("Rdata", "moura2021_BM_tmb.rds"))

head(confint(BM_tmb), 10)

# (Intercept)                                                            0.1132610
# Std.Dev.(Intercept)|study.id                                           0.1056462
# Std.Dev.(Intercept)|species.id.phy                                     0.1932854
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1            0.1293612
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1               0.1293612
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1             0.1293612
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                    0.1293612
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1 0.1293612
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                     0.1293612
# Std.Dev.species.id.phyAgalychnis_callidryas_ott9483|g.1                0.1293612
#                                                                           97.5 %
# (Intercept)                                                            0.6230707
# Std.Dev.(Intercept)|study.id                                           0.1813458
# Std.Dev.(Intercept)|species.id.phy                                     0.2879763
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1            0.3959730
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1               0.3959730
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1             0.3959730
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                    0.3959730
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1 0.3959730
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                     0.3959730
# Std.Dev.species.id.phyAgalychnis_callidryas_ott9483|g.1                0.3959730
#                                                                         Estimate
# (Intercept)                                                            0.3681658
# Std.Dev.(Intercept)|study.id                                           0.1384142
# Std.Dev.(Intercept)|species.id.phy                                     0.2359271
# Std.Dev.species.id.phyAcanthurus_leucosternon_ott388125|g.1            0.2263262
# Std.Dev.species.id.phyAcanthurus_nigricans_ott467313|g.1               0.2263262
# Std.Dev.species.id.phyAcanthurus_nigrofuscus_ott605289|g.1             0.2263262
# Std.Dev.species.id.phyAchatina_fulica_ott997087|g.1                    0.2263262
# Std.Dev.species.id.phyAegithalos_glaucogularis_vinaceus_ott5560982|g.1 0.2263262
# Std.Dev.species.id.phyAethia_pusilla_ott855484|g.1                     0.2263262
# Std.Dev.species.id.phyAgalychnis_callidryas_ott9483|g.1                0.2263262

sigma(BM_tmb)^2
# 0.01445014
# summary(BM_tmb)
BM_tmb_varcor <- VarCorr(BM_tmb)$cond
head(BM_tmb_varcor, 10)

## brms ----

fit_brms1 <- bf(yi | se(sqrt(vi)) ~ 1 + 
                  (1 | study.id) + 
                  (1 | effect.size.id) + 
                  (1 | species.id) + 
                  (1 | gr(species.id.phy, cov = A))
)

prior_brms1 <- get_prior(
  formula = fit_brms1,
  data = dat_moura2021,
  data2 = list(A = A),
  family = gaussian()
)

max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)
options(mc.cores = num_chains) 

m_brms1 <- brm(
  formula = fit_brms1,
  family = gaussian(),
  data = dat_moura2021,
  data2 = list(A = A),
  prior = prior_brms1,
  iter = 10000, 
  warmup = 7000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
# saveRDS(m_brms1, here("Rdata", "moura2021_BM_brms.rds"))
summary(m_brms1)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: yi | se(sqrt(vi)) ~ 1 + (1 | study.id) + (1 | effect.size.id) + (1 | species.id) + (1 | gr(species.id.phy, cov = A)) 
#    Data: dat_moura2021 (Number of observations: 1828) 
#   Draws: 2 chains, each with iter = 10000; warmup = 7000; thin = 1;
#          total post-warmup draws = 6000

# Multilevel Hyperparameters:
# ~effect.size.id (Number of levels: 1828) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.12      0.01     0.11     0.13 1.00     1702     3408

# ~species.id (Number of levels: 341) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.23      0.03     0.17     0.28 1.00      636      912

# ~species.id.phy (Number of levels: 341) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.28      0.09     0.15     0.51 1.00      827      784

# ~study.id (Number of levels: 457) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.14      0.02     0.10     0.18 1.00      588     1121

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.37      0.16     0.03     0.68 1.00     2163     2207

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA

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
|model      |    logLik|       est|        se|  sigma2.u|  sigma2.m| sigma2.non.phylo| sigma2.phylo|
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
