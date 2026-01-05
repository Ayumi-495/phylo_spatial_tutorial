pacman::p_load(brms, metafor, metadat, tidyverse, data.table, crayon, here, ape,
               purrr, stringr, readr, lubridate, magrittr, janitor, glmmTMB, rotl)

## this is necessary the first time you run meta-analysis by glmmTMB
remove.packages("glmmTMB")
remotes::install_github("coraliewilliams/glmmTMB",
                        ref = "equalto_covstruc",
                        subdir = "glmmTMB")
library(glmmTMB)
# # check the available covstruc (you should see equalto in number 13)
glmmTMB:::.valid_covstruct


rm(list = ls())


# scholer 2020 ----
dat_scholer <- read.csv(here("data", "Scholer_2020", "Scholer_2020.csv"))
dat_scholer$effect_id <- seq_len(nrow(dat_scholer))
dat_scholer$non_phylo <- dat_scholer$tip_label

tr_scholer <- read.nexus(here("data", "Scholer_2020", "Scholer.nex"))
tree <- tr_scholer[[1]]

## check species names ----

length(setdiff(dat_scholer$tip_label, tree$tip.label))
length(setdiff(tree$tip.label, dat_scholer$tip_label))
# plot(tree)

### prune data/tree to match the tree/data  ----
data_species <- unique(dat_scholer$tip_label)
dat_scholer <- dat_scholer |> 
  filter(tip_label %in% tree$tip.label)

drop_species <- setdiff(tree$tip.label, unique(dat_scholer$tip_label))
tree <- drop.tip(tree, drop_species)

A <- vcv.phylo(tree, corr = TRUE)
A <- A[sort(rownames(A)), sort(rownames(A))]

dat_scholer$g <- 1 
dat_scholer$var <- (dat_scholer$se)^2

VCV <- diag(dat_scholer$var, nrow = nrow(dat_scholer)) 
rownames(VCV)<- colnames(VCV)<- dat_scholer$effect_id

dat_scholer$effect_id <- as.factor(dat_scholer$effect_id)
head(dat_scholer)


BM_tmb <- glmmTMB(logit_survival ~ 1 + equalto(0 + effect_id|g, VCV) +
                         (1|ref) + 
                         (1|non_phylo) +
                         propto(0 + tip_label|g, A),
                       data = dat_scholer,
                       REML = T)

head(confint(BM_tmb), 10)
sigma(BM_tmb)^2
summary(BM_tmb)
BM_tmb_varcor <- VarCorr(BM_tmb)$cond
head(BM_tmb_varcor, 10)
# 
est = unlist(fixef(BM_tmb))[[1]] #overal mean
se = as.numeric(sqrt(vcov(BM_tmb)[[1]])) #overall mean SE
sigma2.u = BM_tmb_varcor$ref[1] ##among study variance estimate
sigma2.non.phylo = BM_tmb_varcor$non_phylo[1] ##non-phylo variance estimate
sigma2.phylo = BM_tmb_varcor$g.1[1] ##phylo variance component
sigma2.m = sigma(BM_tmb) ##within study variance estimate
sigma2.total = sum(sigma2.u+sigma2.non.phylo+sigma2.phylo+sigma2.m)


metafor_1 <- data.frame(model = "BM_metafor", 
                      logLik = logLik(BM),
                      est = BM$b[[1]], 
                      se = BM$se[[1]], 
                      sigma2.u = sqrt(BM$sigma2[1]), ## among study variance
                      sigma2.m = sqrt(BM$sigma2[2]), ## within study variance
                      sigma2.non.phylo = sqrt(BM$sigma2[3]),
                      sigma2.phylo = sqrt(BM$sigma2[4])
                      )


glmmTMB_1 <- data.frame(model = "BM_tmb",
                      logLik = logLik(BM_tmb)[1],
                      est = unlist(fixef(BM_tmb))[[1]], #overall mean
                      se = as.numeric(sqrt(vcov(BM_tmb)[[1]])), #overall mean SE
                      sigma2.u = sqrt(BM_tmb_varcor$ref[1]), ##among study variance estimate
                      sigma2.m = sigma(BM_tmb), ##within study variance estimate
                      sigma2.non.phylo = sqrt(BM_tmb_varcor$non_phylo[1]), ##non-phylo variance estimate
                      sigma2.phylo = sqrt(BM_tmb_varcor$g.1[1]) ##phylo variance component
                      )

output <- rbind(metafor_1, glmmTMB_1)
knitr::kable(output)            
  
# Huberty_Denno_2004 ----
Huberty_Denno_2004 <- read_csv(here("data", "Huberty_Denno_2004", "Huberty_Denno_2004_cleaned.csv"))
HubDen_tree <- read.tree(here("data", "Huberty_Denno_2004", "Huberty_Denno_2004_tree.txt"))

tip_labels <- HubDen_tree$tip.label
print(tip_labels)
unique(Huberty_Denno_2004$Herbivore)

setdiff(tip_labels, Huberty_Denno_2004$Herbivore)
# [1] "Semiothesia_sp."     "Cameraria_sp."       "Multireoides_sp."    "Phytocoris_nigrinus" "Parthenicus_sp."  "Aphalara_sp."  

Huberty_Denno_2004$Herbivore <- gsub("(?<=_sp)$", ".", Huberty_Denno_2004$Herbivore, perl = TRUE)
tips_to_keep <- intersect(Huberty_Denno_2004$Herbivore, tip_labels)
HubDen_tree_pruned <- drop.tip(HubDen_tree, setdiff(tip_labels, tips_to_keep))

setdiff(HubDen_tree_pruned$tip.label, Huberty_Denno_2004$Herbivore)


set.seed(42)
tr_bin <- multi2di(HubDen_tree_pruned, random = TRUE)

if (is.null(tr_bin$edge.length)) {
  tr_bin$edge.length <- rep(1, nrow(tr_bin$edge))
}

pos_min <- min(tr_bin$edge.length[tr_bin$edge.length > 0], na.rm = TRUE)
if (!is.finite(pos_min)) pos_min <- 1
tr_bin$edge.length[is.na(tr_bin$edge.length) | tr_bin$edge.length == 0] <- pos_min * 1e-6

is.ultrametric(tr_bin)

## ultrametric
tr_um <- compute.brlen(tr_bin, method = "Grafen")
any(tr_um$edge.length <= 0, na.rm = TRUE) 

## check
is.binary(tr_um)
is.ultrametric(tr_um)
plot(tr_um)

A <- vcv(tr_um, corr = TRUE)

# change the variable name (sub-guild) to sub_guild

Huberty_Denno_2004 <- Huberty_Denno_2004 |> 
  rename(
    Study_id = Reference, 
    Species = Herbivore,
    Sub_guild = `Sub-guild`
  ) |> 
  mutate(Phylo = Species,
         Effect_id = row_number(),
         Feeding_mode = ifelse(Sub_guild == "Free-living", "Free-living", "Endophytic")
  )


Huberty_Denno_2004_1 <- Huberty_Denno_2004 %>%
  filter(Phylo %in% rownames(A)) %>%
  mutate(Phylo = factor(Phylo))

sp_in_both <- intersect(rownames(A), Huberty_Denno_2004$Phylo) # 49 species

Huberty_Denno_2004_1 <- Huberty_Denno_2004 %>%
  filter(Phylo %in% sp_in_both) %>%
  droplevels()
A <- A[sp_in_both, sp_in_both]
Huberty_Denno_2004_1$Phylo <- factor(Huberty_Denno_2004_1$Phylo, levels = rownames(A))

## run model ----
Huberty_Denno_2004_1$g <- 1 
V <- vcalc(var, cluster = Study_id, obs = Effect_id, data = Huberty_Denno_2004_1, rho=0.2) #'[how to define rho? is this correct way??]

colnames(V) <- 1:nrow(Huberty_Denno_2004_1)
row.names(V) <- 1:nrow(Huberty_Denno_2004_1)

dat1_BM_tmb <- glmmTMB(d ~ 1 + equalto(0 + Effect_id|g, var) + 
                         (1|Study_id) + 
                         (1|Species) +
                         propto(0 + Phylo|g, A),
                       data = Huberty_Denno_2004_1,
                       REML = T)

