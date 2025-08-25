
pacman::p_load(brms, MCMCglmm, metafor, metadat, tidyverse, 
               dplyr, tidyr, purrr, stringr, readr, magrittr, 
               janitor, flextable, kableExtra, geosphere, here, rotl)

# Use the data from the metadat package ----
dat <- dat.moura2021$dat
head(dat)

dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)
tree <- dat.moura2021$tree
tree <- compute.brlen(tree, method="Grafen", power=1)

A_BM <- vcv(tree, corr=TRUE)
dat$phy <- dat$species.id
summary(dat)

# fit BM model ----
mod_BM <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                  R = list(phy = A_BM), 
                  data = dat)
summary(mod_BM)

# spatial model (exponential) ----
dat$const <- 1　# add a constant column for the spatial model

A_BM[1:5, 1:5]
I <- matrix(1, nrow = 341, ncol = 341)

D <- I - A_BM # make distance matrix
D[1:20, 1:20]
max(D) # maximum distance is 1
mod_spatial <- rma.mv(yi, vi,
                      random = list(~ 1 | study.id, 
                                    ~ 1 | effect.size.id, 
                                    ~ 1 | species.id, 
                                    ~ phy|const),
                      dist=list(D),
                      struct = "SPEXP",
                      control=list(rho.init=1), 
                      data = dat)
summary(mod_spatial)

# OU model ----
# Tmat <- cophenetic(tree) #computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
rho <- mod_spatial$rho
A_OU2 <- exp(-D / rho)
mod_OU <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                 R = list(phy = A_OU2), 
                 data = dat)
summary(mod_OU)
## -- ##
alpha <- 1/rho
A_OU <- exp(-alpha * D)
mod_OU2 <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                 R = list(phy = A_OU), 
                 data = dat)
summary(mod_OU2)

# spatial model (gaussian)
mod_gau <- rma.mv(yi, vi,
                  random = list(~ 1 | study.id, 
                                ~ 1 | effect.size.id, 
                                ~ 1 | species.id, 
                                ~ phy|const),
                  dist=list(D),
                  struct = "SPGAU",
                  control=list(rho.init=1), 
                  data = dat)
summary(mod_gau)
max(D)

# linear model
system.time(mod_linear <- rma.mv(yi, vi,
                      random = list(~ 1 | study.id, 
                                    ~ 1 | effect.size.id, 
                                    ~ 1 | species.id, 
                                    ~ phy|const),
                      dist=list(D),
                      struct = "SPLIN",
                      control=list(rho.init=1 ,
                                   optimizer = "optim"
                                   ), # Optimizer (nlminb) did not achieve convergence (convergence = 1).
                      data = dat)
)
summary(mod_linear)

system.time(mod_BM2 <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                 R = list(phy = A_BM), 
                 control = list(optimizer = "optim"), 
                 data = dat)
)
summary(mod_BM2)

# gaussian process model in brms ----
dat <- dat.moura2021$dat
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)
tree <- dat.moura2021$tree
tree <- compute.brlen(tree, method="Grafen", power=1)

A_BM <- vcv(tree, corr=TRUE)
dat$phy <- dat$species.id
summary(dat)

A_BM[1:5, 1:5]
I <- matrix(1, nrow = 341, ncol = 341)

D <- I - A_BM # make distance matrix
D[1:20, 1:20]

formula_1 <- bf(yi | se(vi) ~ 1 + gp(D, cov = "exp_quad", gr = TRUE))

prior_ma1 <- get_prior(
  formula = formula_1,
  data = dat,
  data2 = list(D = D),
  family = gaussian()
)

system.time(fit_brm1 <- brm(
  formula_1,
  data = dat,
  data2 = list(D = D),
  chains = 2, 
  iter = 3000, 
  warmup = 2000,
  prior = prior_ma1
))
