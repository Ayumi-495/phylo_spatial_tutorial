pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape, sf,
                purrr, stringr, readr, lubridate, magrittr, janitor, 
              　geosphere, rotl)
rm(list = ls())
dat_scholer <- read.csv(here("data", "Scholer_2020", "Scholer_2020.csv"))
dat_scholer$effect_id <- seq_len(nrow(dat_scholer))
dat_scholer$non_phylo <- dat_scholer$tip_label

tr_scholer <- read.nexus(here("data", "Scholer_2020", "Scholer.nex"))
tree <- tr_scholer[[1]]

# phylogenetic meta-analysis ----
## check species names ----

length(setdiff(dat_scholer$tip_label, tree$tip.label))
length(setdiff(tree$tip.label, dat_scholer$tip_label))
plot(tree)

### prune data/tree to match the tree/data  ----
data_species <- unique(dat_scholer$tip_label)
dat_scholer <- dat_scholer |> 
  filter(tip_label %in% tree$tip.label)

drop_species <- setdiff(tree$tip.label, unique(dat_scholer$tip_label))
tree <- drop.tip(tree, drop_species)

A <- vcv.phylo(tree, corr = TRUE)

head(dat_scholer)


## BM model ----
### single tree ----
BM <- rma.mv(logit_survival, se^2,
                random = list(~ 1 | ref, 
                              ~ 1 | effect_id, 
                              ~ 1 | non_phylo, 
                              ~ 1 | tip_label),
                R = list(tip_label = A), 
                data = dat_scholer,
                sparse = TRUE,
                verbose = TRUE)

summary(BM)
# Multivariate Meta-Analysis Model (k = 821; method: REML)
# Multivariate Meta-Analysis Model (k = 821; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -607.6516  1215.3032  1225.3032  1248.8498  1225.3769   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed     factor    R 
# sigma^2.1  0.1213  0.3483    189     no        ref   no 
# sigma^2.2  0.1177  0.3430    821     no  effect_id   no 
# sigma^2.3  0.0093  0.0966    530     no   non_phylo   no 
# sigma^2.4  0.3852  0.6207    530     no  tip_label  yes 
# 
# Test for Heterogeneity:
#   Q(df = 820) = 910083.3991, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# 0.7559  0.1742  4.3393  <.0001  0.4145  1.0973  *** 

orchaRd::i2_ml(BM)
# I2_Total       I2_ref I2_effect_id  I2_non_phyo I2_tip_label 
# 99.890015    19.127454    18.550809     1.471464    60.740288 

confint(BM)

dat_scholer$cMass <- scale(log(dat_scholer$mass), center = TRUE, scale = FALSE)

BM_1 <- rma.mv(logit_survival, se^2,
               mod = ~ + cMass,
             random = list(~ 1 | ref, 
                           ~ 1 | effect_id, 
                           ~ 1 | non_phylo, 
                           ~ 1 | tip_label),
             R = list(tip_label = A), 
             data = dat_scholer)

summary(BM_1)
orchaRd::i2_ml(BM_1)


### multiple tree ----
head(dat_scholer)
phy_model <-  function(cor_tree = vcv_tree){
  model <- rma.mv(logit_survival, se^2,
                  random = list(~ 1 | ref, 
                                ~ 1 | effect_id, 
                                ~ 1 | non_phylo, 
                                ~ 1 | tip_label),
                  R = list(tip_label = cor_tree), 
                  data = dat_scholer)
  }


tree_50 <- tr_scholer[1:50]
species_keep <- unique(dat_scholer$tip_label)

tree_50_pruned <- lapply(tree_50, function(tr) {
  drop_species <- setdiff(tr$tip.label, species_keep)
  ape::drop.tip(tr, drop_species)
})

class(tree_50_pruned) <- "multiPhylo"

vcv_tree50 <- purrr::map(tree_50_pruned, ~vcv(.x, corr = TRUE))
ma_50 <- parallel::mclapply(vcv_tree50, phy_model, mc.cores = 4)

saveRDS(ma_50, here("Rdata", "ma_50.RDS")) 

length(unique(tree_50_pruned[[1]]$tip.label))
length(unique(dat_scholer$tip_label))

## spatial exponential model to get alpha value ----

# add a constant column for specifying a shared correlation structure
dat_scholer$Const <- 1

#  make distance matrix
## by subtracting the correlation matrix from the identity matrix (i.e., distance = 1 - correlation)
I <- matrix(1, nrow = 530, ncol = 530)
## convert the phylogenetic correlation matrix into a distance matrix
D <- I - A 
D[1:5, 1:5]  # Check the first few rows and columns
dat_scholer$const <- 1

head(dat_scholer)


spatial <- rma.mv(logit_survival, se^2,
                  random = list(~ 1 | ref,  # random effect for study
                            ~ 1 | effect_id, # random effect for each effect size
                            ~ 1 | non_phylo, # random effect for species
                            ~ tip_label | const), # random effect for phylogeny
                  dist = list(D),  # phylogenetic distance matrix
                  struct = "SPEXP", # use spatial exponential correlation structure
                  control = list(rho.init = 1), 
                  data = dat_scholer)
summary(spatial)
# Multivariate Meta-Analysis Model (k = 821; method: REML)
# 
# logLik   Deviance        AIC        BIC       AICc   
# -607.6517  1215.3034  1227.3034  1255.5592  1227.4067   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed     factor 
# sigma^2.1  0.1213  0.3483    189     no        ref 
# sigma^2.2  0.1177  0.3430    821     no  effect_id 
# sigma^2.3  0.0093  0.0966    530     no  non_phylo 
# 
# outer factor: const      (nlvls = 1)
# inner term:   ~tip_label (nlvls = 530)
# 
# estim     sqrt  fixed 
# tau^2       4693.7019  68.5106     no 
# rho        12183.2344              no 
# 
# Test for Heterogeneity:
#   Q(df = 820) = 910083.3991, p-val < .0001
# 
# Model Results:
#   
#   estimate       se    zval    pval      ci.lb     ci.ub    
# 0.7559  68.5080  0.0110  0.9912  -133.5174  135.0291   



## OU model ----
# rho <- spatial$rho
alpha <- 0.9
A_OU <- exp(-alpha * D)
OU <- rma.mv(logit_survival, se^2,
                random = list(~ 1 | ref, 
                              ~ 1 | effect_id, 
                              ~ 1 | non_phylo, 
                              ~ 1 | tip_label),
                R = list(tip_label = A_OU), 
                data = dat_scholer
             )

summary(OU)
# ultivariate Meta-Analysis Model (k = 821; method: REML)
# 
# 
# logLik   Deviance        AIC        BIC       AICc   
# -607.6517  1215.3034  1225.3034  1248.8499  1225.3771   
# 
# Variance Components:
#   
#   estim     sqrt  nlvls  fixed     factor    R 
# sigma^2.1     0.1213   0.3483    189     no        ref   no 
# sigma^2.2     0.1177   0.3430    821     no  effect_id   no 
# sigma^2.3     0.0093   0.0966    530     no  non_phylo   no 
# sigma^2.4  4693.6375  68.5101    530     no  tip_label  yes 
# 
# Test for Heterogeneity:
#   Q(df = 820) = 910083.3991, p-val < .0001
# 
# Model Results:
#   
#   estimate       se    zval    pval      ci.lb     ci.ub    
# 0.7559  68.5075  0.0110  0.9912  -133.5164  135.0282    

I <- matrix(1, nrow = 530, ncol = 530)
D <- I - A 

alphas <- c(0.9, 0.7, 0.5, 0.3, 0.1, 0.05, 0.01)
out <- vector("list", length(alphas))

fmt_dur <- function(sec) sprintf("%02d:%02d:%02d",
                                floor(sec/3600),
                                floor((sec %% 3600)/60),
                                round(sec %% 60))

pb <- txtProgressBar(min = 0, max = length(alphas), style = 3)
t0 <- Sys.time()

for (i in seq_along(alphas)) {
  a <- alphas[i]
  t_iter0 <- Sys.time()

  A <- exp(-a * D)

  fit <- rma.mv(
    yi     = logit_survival,
    V      = se^2,
    random = list(~ 1 | ref,
                  ~ 1 | effect_id,
                  ~ 1 | non_phylo,
                  ~ 1 | tip_label),
    R      = list(tip_label = A),
    data   = dat_scholer,
    method = "REML"
  )

  out[[i]] <- data.frame(
    alpha        = a,
    k            = fit$k,
    logLik       = as.numeric(logLik(fit)),
    AIC          = AIC(fit),
    BIC          = BIC(fit),
    beta0        = unname(coef(fit))[1],
    se_beta0     = fit$se[1],
    sigma2_ref        = fit$sigma2[1],
    sigma2_effect_id  = fit$sigma2[2],
    sigma2_non_phylo  = fit$sigma2[3],
    sigma2_tip_label  = fit$sigma2[4],
    row.names = NULL, check.names = FALSE
  )

  setTxtProgressBar(pb, i)
  t_iter   <- as.numeric(difftime(Sys.time(), t_iter0, units = "secs"))
  t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eta_sec  <- (t_elapsed / i) * (length(alphas) - i)
  cat(sprintf("\r  [%d/%d] alpha=%.2f | iter=%s | elapsed=%s | ETA=%s",
              i, length(alphas), a, fmt_dur(t_iter),
              fmt_dur(t_elapsed), fmt_dur(eta_sec)))
  flush.console()

}
close(pb)

results_ou <- do.call(rbind, out)
library(knitr)
kable(results_ou)



### brms ----
norm_utf8 <- function(x) iconv(x, from = "", to = "UTF-8", sub = " ")

clean_whitespace <- function(x) {
  x <- norm_utf8(as.character(x))
    x <- gsub("[\\p{Z}\\p{C}\u00A0\u202F\u2007\u2009\u200A\u200B\u200C\u200D\u2060\uFEFF]+",
            " ", x, perl = TRUE)
  trimws(x)
}

dat_scholer <- transform(
  dat_scholer,
  ref       = factor(clean_whitespace(ref)),
  effect_id = factor(clean_whitespace(effect_id)),
  species   = factor(clean_whitespace(species)),
  tip_label = factor(clean_whitespace(tip_label))
)

fit_brms1 <- bf(logit_survival | se(se) ~ 1 + 
                  (1 | ref) + 
                  (1 | effect_id) + 
                  (1 | species) + 
                  (1 | gr(tip_label, cov = A))
)

prior_brms1 <- get_prior(
  formula = fit_brms1,
  data = dat_scholer,
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
  data = dat_scholer,
  data2 = list(A = A),
  prior = prior_brms1,
  iter = 10000, 
  warmup = 7000,  
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain), 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(m_brms1)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: logit_survival | se(se) ~ 1 + (1 | ref) + (1 | effect_id) + (1 | species) + (1 | gr(tip_label, cov = A)) 
#    Data: dat_scholer (Number of observations: 821) 
#   Draws: 2 chains, each with iter = 10000; warmup = 7000; thin = 1;
#          total post-warmup draws = 6000

# Multilevel Hyperparameters:
# ~effect_id (Number of levels: 821) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.34      0.01     0.32     0.37 1.00     1824     3868

# ~ref (Number of levels: 189) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.35      0.04     0.28     0.43 1.00     2753     3841

# ~species (Number of levels: 540) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.09      0.05     0.00     0.19 1.00      368      691

# ~tip_label (Number of levels: 530) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.62      0.05     0.53     0.72 1.00     1655     2404

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.76      0.18     0.42     1.10 1.00     5286     4903

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.00      0.00     0.00     0.00   NA       NA       NA
pairs(m_brms1)
# spatial meta-analysis ----
## organise data ----
## make Euclidian distance matrix ----
coords <- cbind(dat_scholer$lat, dat_scholer$long)

dat_sf <- st_as_sf(dat_scholer, coords = c("long", "lat"), crs = 4326)
dat_sf_proj <- st_transform(dat_sf, crs = 3857)

coords_m <- st_coordinates(dat_sf_proj)
dat_scholer$x_km <- coords_m[,1] / 1000
dat_scholer$y_km <- coords_m[,2] / 1000

coords_km <- cbind(dat_scholer$x_km, dat_scholer$y_km)
dist_matrix_euclid <- as.matrix(dist(coords_km))

rownames(dist_matrix_euclid) <- dat_scholer$datapt_id
colnames(dist_matrix_euclid) <- dat_scholer$datapt_id

# rename and modify the variables
dat_scholer <- dat_scholer |> 
  rename(study_id = ref,
         effect_id = X)

dat_scholer$const <- 1
dat_scholer$vi <- dat_scholer$se^2




# meta-analysis ----
## exponential model ----
### metafor ----

names(dat_scholer)

EXP_m1 <- rma.mv(logit_survival, vi, 
                 random = list(
                   ~ 1 | study_id,
                   ~ 1 | effect_id, 
                   ~ effect_id | const
                 ), 
                 struct = "SPEXP", 
                 data = dat_scholer,
                 dist = list(effect_id = dist_matrix_euclid)
                 )

summary(EXP_m1)

EXP_m2 <- rma.mv(logit_survival, vi, 
                 random = list(
                   ~ 1 | study_id,
                   ~ 1 | effect_id, 
                   ~ x_km + y_km | const
                 ), 
                 struct = "SPEXP", 
                 data = dat_scholer
                 )

### brms ----

fit_EXP <- bf(
  logit_survival | se(se) ~ 1 +
    (1 | study_id) + 
    (1 | effect_id) + 
    gp(x_km, y_km, 
       cov = "exponential",
       scale = FALSE) 
    
) 


