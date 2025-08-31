## FID data - did not converge the spatial model… :/

dat_FID <- read.csv(here("data", "FID.csv"))
names(dat_FID)

dat_FID$EffectID <- seq_len(nrow(dat_FID))
dat_FID$Phylo <- dat_FID$Species

tr_FIDs <- read.nexus(here("data", "FID.nex"))
tr_FID <- tr_FIDs[[1]]

dat_FID$Phylo <- gsub(" ", "_", dat_FID$Phylo)
length(setdiff(dat_FID$Phylo, tr_FID$tip.label))
length(setdiff(tr_FID$tip.label, dat_FID$Phylo))

A_FID <- vcv.phylo(tr_FID, corr = TRUE)


# BM model ----
## metafor ----
BM_FID <- rma.mv(EffectSize, Var,
                random = list(~ 1 | EffectID,
                            ~ 1 | StudyID,
                            ~ 1 | Species,
                            ~ 1 | Phylo),
                R = list(Phylo = A_FID), 
                data = dat_FID,
                sparse = TRUE,
                verbose = TRUE
                )
summary(BM_FID)                
i2 <- orchaRd::i2_ml(BM_FID)
round(i2, 2)

## OU model ----
# spatial exp model
dat_FID$Const <- 1
# length(unique(dat_FID$Species)) # 99 species
I_FID <- matrix(1, nrow = 99, ncol = 99)
## convert the phylogenetic correlation matrix into a distance matrix
D_FID <- I_FID - A_FID

spatial_FID <- rma.mv(EffectSize, Var,
                  random = list(~ 1 | StudyID,  # random effect for study
                            ~ 1 | EffectID, # random effect for each effect size
                            ~ Phylo | Const), # random effect for phylogeny
                  dist = list(D_FID),  # phylogenetic distance matrix
                  struct = "SPEXP", # use spatial exponential correlation structure
                  control = list(rho.init = 1), 
                  data = dat_FID,
                  sparse = TRUE,
                  verbose = TRUE,
                  method = "REML"
                  )
