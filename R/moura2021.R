
pacman::p_load(brms, MCMCglmm, metafor, metadat, tidyverse, 
               dplyr, tidyr, purrr, stringr, readr, magrittr, 
               janitor, flextable, kableExtra, geosphere, here, rotl)

# Use the data from the metadat package 
dat <- dat.moura2021$dat
head(dat)

dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)
tree <- dat.moura2021$tree
tree <- compute.brlen(tree, method="Grafen", power=1)

A_BM <- vcv(tree, corr=TRUE)
dat$phy <- dat$species.id
summary(dat)

### fit BM model
mod_BM <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                  R = list(phy = A_BM), 
                  data = dat)
summary(mod_BM)

### spatial model (exponential)
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

## OU model
alpha <- 1/0.0182　# alpha = 1/rho - rho = 0.0182
Tmat <- cophenetic(tree) #computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
format_back <- function(x) {
  x <- tolower(x)
  x <- gsub(" ", "_", x)
  return(x)
}
rownames(Tmat) <- format_back(rownames(Tmat))
colnames(Tmat) <- format_back(colnames(Tmat))
Tmat <- Tmat[rownames(A_BM), colnames(A_BM)]

A_OU <- exp(-alpha * Tmat)

mod_OU <- rma.mv(yi, vi,
                 random = list(~ 1 | study.id, 
                               ~ 1 | effect.size.id, 
                               ~ 1 | species.id, 
                               ~ 1 | phy),
                 R = list(phy = A_OU), 
                 data = dat)
summary(mod_OU)

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
