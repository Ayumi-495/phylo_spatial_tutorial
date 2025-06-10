pacman::p_load(brms, MCMCglmm, metafor, metadat, tidyverse, 
               dplyr, tidyr, purrr, stringr, readr, magrittr, 
               janitor, flextable, kableExtra, geosphere, here, ape)

Huberty_Denno_2004 <- read_csv(here("data", "Huberty_Denno_2004", "Huberty_Denno_2004_cleaned.csv"))
HubDen_tree <- read.tree(here("data", "Huberty_Denno_2004", "Huberty_Denno_2004_tree.txt"))

tip_labels <- HubDen_tree$tip.label
print(tip_labels)
unique(Huberty_Denno_2004$Herbivore)

setdiff(tip_labels, Huberty_Denno_2004$Herbivore)
# [1] "Semiothesia_sp."     "Cameraria_sp."       "Multireoides_sp."    "Phytocoris_nigrinus" "Parthenicus_sp."    
# [6] "Aphalara_sp."  

Huberty_Denno_2004$Herbivore <- gsub("(?<=_sp)$", ".", Huberty_Denno_2004$Herbivore, perl = TRUE)
tips_to_keep <- intersect(Huberty_Denno_2004$Herbivore, tip_labels)
HubDen_tree_pruned <- drop.tip(HubDen_tree, setdiff(tip_labels, tips_to_keep))

setdiff(HubDen_tree_pruned$tip.label, Huberty_Denno_2004$Herbivore)

# make variance covariance matrix
HubDen_vcv <- vcv.phylo(HubDen_tree_pruned, corr = TRUE)
Huberty_Denno_2004_pruned <- Huberty_Denno_2004 %>%
  filter(Herbivore %in% rownames(HubDen_vcv))

length(unique(Huberty_Denno_2004_pruned$Herbivore)) 
length(HubDen_tree_pruned$tip.label)

HubDen_2004 <- Huberty_Denno_2004 %>%
  filter(Herbivore %in% rownames(HubDen_vcv)) %>%
  mutate(Herbivore = factor(Herbivore)) 
HubDen_2004$Phylo <- HubDen_2004$Herbivore

HubDen_2004$effect_id <- seq_len(nrow(HubDen_2004))
View(HubDen_2004)
names(HubDen_2004)
# meta-analysis ----

## BM model ----
HD_BM <- rma.mv(d, var,
                 random = list(~ 1 | Reference, 
                               ~ 1 | effect_id, 
                               ~ 1 | Herbivore, 
                               ~ 1 | Phylo),
                 R = list(Phylo = HubDen_vcv), 
                 data = HubDen_2004)

summary(HD_BM)

## spatial model (exponential) ----
HubDen_2004$Const <- 1
I <- matrix(1, nrow = 54, ncol = 54)

D_HD <- I - HubDen_vcv # make distance matrix

HD_spatial <- rma.mv(d, var,
                     random = list(~ 1 | Reference, 
                                   ~ 1 | effect_id, 
                                   ~ 1 | Herbivore, 
                                   ~ Phylo | Const),
                     dist=list(D_HD), 
                     struct = "SPEXP",
                     control=list(rho.init=1), 
                     data = HubDen_2004)

summary(HD_spatial)

## OU model ----
rho <- HD_spatial$rho
alpha <- 1/rho
HD_OU <- exp(-alpha * D_HD)

HD_OU <- rma.mv(d, var,
                random = list(~ 1 | Reference, 
                              ~ 1 | effect_id, 
                              ~ 1 | Herbivore, 
                              ~ 1 | Phylo),
                  R = list(Phylo = HD_OU), 
                  data = HubDen_2004
                )

summary(HD_OU)

## GP model ----
HD_GP <- rma.mv(d, var,
                random = list(~ 1 | Reference, 
                              ~ 1 | effect_id, 
                              ~ 1 | Herbivore, 
                          ~ Phylo | Const),
                     dist=list(D_HD), 
                     struct = "SPGAU",
                     control=list(rho.init=1), 
                     data = HubDen_2004)

summary(HD_GP)

## linear model ----
HD_LIN <- rma.mv(d, var,
                random = list(~ 1 | Reference, 
                              ~ 1 | effect_id, 
                              ~ 1 | Herbivore, 
                              ~ Phylo | Const),
                dist=list(D_HD), 
                struct = "SPLIN",
                control=list(rho.init=1,
                             optimizer = "optim"), 
                data = HubDen_2004)

summary(HD_LIN)

