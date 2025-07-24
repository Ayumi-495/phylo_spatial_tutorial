
pacman::p_load(brms, MCMCglmm, metafor, metadat, tidyverse, 
               dplyr, tidyr, purrr, stringr, readr, magrittr, 
               janitor, flextable, kableExtra, geosphere, here, rotl)

dat_pietro <- read.csv(here("data", "pietro_2022", "raw_data.csv"))
source(here("data", "pietro_2022","calculate_smd.R")) # function that calculate effect sizes

mod_data <- calculate_smd(dat_pietro)
mod_data <-
  mod_data %>%
  filter(!(m_quality_status == "assumed" &
             m_trait_table1 == "column2"))
hyp_1 <- 
  mod_data %>% 
  filter(m_quality_status == "verified")

hyp_1$const <- 1　# add a constant column for the spatial model

# make correlation matrix for the phylogenetic tree ----
# 1. get species names and match them with Open Tree of Life
taxa <- tnrs_match_names(unique(hyp_1$species),
                         context_name = "Animals",
                         do_approximate_matching = TRUE)

# 2. fix taxon names that are not matched correctly
# dermestes_maculatus
taxa$ott_id[taxa$search_string == "dermestes_maculatus"] <- 1036466

# drosophila_melanogaster
taxa$ott_id[taxa$search_string == "drosophila_melanogaster"] <- 505714

# 3. connect taxon names with species names
taxon_map <- structure(taxa$search_string,
                       names = taxa$unique_name)

# 4. get the phylogenetic tree
tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
tr$tip.label <- taxon_map[otl_tips]
tr$node.label <- NULL
tr <- compute.brlen(tr)

# 5. make a correlation matrix for the phylogenetic tree
A <- vcv(tr, cor = TRUE)
rownames(A) <- tolower(gsub(" ", "_", rownames(A)))
colnames(A) <- tolower(gsub(" ", "_", colnames(A)))

# 6. label the chips
tr$tip.label <- str_replace(str_to_title(taxon_map[otl_tips]), "_", " ")

D <- 1 - A # make 
D[1:5, 1:5]  # Check the first few rows and columns

# BM model ----
mod1 <- rma.mv(yi = smd,
       V = var,
       random = list(~1|effect_id,
                     ~1|study_id,
                     ~1|experiment_id,
                     ~1|measure_id,
                     ~1|species_id,
                     ~1|species),
       R = list(species = A),
       method = "REML",
       data = hyp_1)
mod1

# spatial exponential model----
mod2 <- rma.mv(yi = smd,
               V = var,
               random = list(~1|effect_id,
                             ~1|study_id,
                             ~1|experiment_id,
                             ~1|measure_id,
                             ~1|species_id,
                             ~ species|const),
               dist=list(D),
               struct = "SPEXP",
               # control=list(rho.init=10), 
               data = hyp_1)
mod2

# OU model ----
alpha <- 0.0166
Tmat <- cophenetic(tr) #computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
format_back <- function(x) {
  x <- tolower(x)
  x <- gsub(" ", "_", x)
  return(x)
}
rownames(Tmat) <- format_back(rownames(Tmat))
colnames(Tmat) <- format_back(colnames(Tmat))
Tmat <- Tmat[rownames(A), colnames(A)]

OU_cor <- exp(-alpha * Tmat)

mod3 <- rma.mv(yi = smd,
               V = var,
               random = list(~1|effect_id,
                             ~1|study_id,
                             ~1|experiment_id,
                             ~1|measure_id,
                             ~1|species_id,
                             ~1|species),
               R = list(species = OU_cor),
               method = "REML",
               data = hyp_1)
mod3
