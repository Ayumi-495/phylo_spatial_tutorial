pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape, sf,
                purrr, stringr, readr, lubridate, magrittr, janitor, 
              　geosphere, rotl)
rm(list = ls())
dat_scholer <- read.csv(here("data", "Scholer_2020", "Scholer_2020.csv"))
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

dat_scholer$effect_id <- seq_len(nrow(dat_scholer))
head(dat_scholer)

## BM model ----
BM <- rma.mv(logit_survival, se^2,
                random = list(~ 1 | ref, 
                              ~ 1 | effect_id, 
                              ~ 1 | species, 
                              ~ 1 | tip_label),
                R = list(tip_label = A), 
                data = dat_scholer)

summary(BM)

## spatial exponential model to get alpha value ----

# add a constant column for specifying a shared correlation structure
dat_scholer$Const <- 1

#  make distance matrix
## by subtracting the correlation matrix from the identity matrix (i.e., distance = 1 - correlation)
I <- matrix(1, nrow = 530, ncol = 530)
length(D[1, ])
## convert the phylogenetic correlation matrix into a distance matrix
D <- I - A 
D_dist <- cophenetic(tree)
D[1:5, 1:5]  # Check the first few rows and columns
D_dist[1:5, 1:5]
dat_scholer$const <- 1

head(dat_scholer)


spatial <- rma.mv(logit_survival, se^2,
                  random = list(~ 1 | ref,  # random effect for study
                            ~ 1 | effect_id, # random effect for each effect size
                            ~ 1 | species, # random effect for species
                            ~ tip_label | const), # random effect for phylogeny
                  dist = list(D),  # phylogenetic distance matrix
                  struct = "SPEXP", # use spatial exponential correlation structure
                  control = list(rho.init = 1), 
                  data = dat_scholer)
summary(spatial)

## OU model ----
rho <- spatial$rho
alpha <- 1/rho
A_OU <- exp(-alpha * D)

OU <- rma.mv(logit_survival, se^2,
                random = list(~ 1 | ref, 
                              ~ 1 | effect_id, 
                              ~ 1 | species, 
                              ~ 1 | tip_label),
                R = list(tip_label = A_OU), 
                data = dat_scholer
             )

summary(OU)




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


