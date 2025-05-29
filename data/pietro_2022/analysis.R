# load packages and necessary functions ----
library(metafor)
library(orchaRd)
library(tidyverse)

source("phylogeny.R") # function that builds phylogenetic tree
source("calculate_smd.R") # function that calculate effect sizes
source("significance.R") # function that performs pairwise comparisons

# NOTES

# 1) numbers in the significance function represent the order of the factors 
# seen in the summary of models; e.g. for model hyp_1_categories, 2 - 1 
# represents the difference in effect size means between high- and low-quality
# males, as high-quality is the first and low-quality is the second row in the
# model's summary

# 2) all models with moderators contain a -1 in the argument mods to remove the
# intercept from the model. If you want to verify the model output with an 
# intercept, simply remove the -1.

# load data ----
raw_data <- read_csv("raw_data.csv")

# calculate effect sizes ----
mod_data <- calculate_smd(raw_data)

# remove effect sizes that belong to "None" category in table 1 ----
mod_data <-
  mod_data %>%
  filter(!(m_quality_status == "assumed" &
             m_trait_table1 == "column2"))

# separating into two datasets, one for each hypothesis ----

# hypothesis 1: male mate choice expression depends on male quality
hyp_1 <- 
  mod_data %>% 
  filter(m_quality_status == "verified")

# hypothesis 2: male mate choice expression depends on male phenotype
hyp_2 <-
  mod_data %>%
  filter(m_trait_table1 == "column1")

# meta-analytical model with all possible random factors ----
set.seed(42)

full.random <- function (df){
  phylo.tree(df)
  
  rma.mv(yi = smd,
         V = var,
         random = list(~1|effect_id,
                       ~1|study_id,
                       ~1|experiment_id,
                       ~1|measure_id,
                       ~1|species_id,
                       ~1|species),
         R = list(species = cor),
         method = "REML",
         data = df)
}

# heterogeneity hypothesis 1
hyp_1_full_random_model <- full.random(hyp_1)
i2_ml(hyp_1_full_random_model, method = "wv") * 100

# heterogeneity hypothesis 2
hyp_2_full_random_model <- full.random(hyp_2)
i2_ml(hyp_2_full_random_model, method = "wv") * 100

# models from now on do not take phylogeny into consideration ----

# hypothesis 1 test ----
# hypothesis 1 general estimate
hyp_1_no_moderators <- rma.mv(yi = smd,
                              V = var,
                              random = list(~1|effect_id,
                                            ~1|study_id,
                                            ~1|experiment_id,
                                            ~1|measure_id,
                                            ~1|species_id),
                              method = "REML",
                              data = hyp_1)
summary(hyp_1_no_moderators)

#heterogeneity
i2_ml(hyp_1_no_moderators, method = "wv") * 100

# hypothesis 1 each male quality category
hyp_1_categories <- rma.mv(yi = smd,
                            V = var,
                            random = list(~1|effect_id,
                                          ~1|study_id,
                                          ~1|experiment_id,
                                          ~1|measure_id,
                                          ~1|species_id),
                            mods = ~ m_quality_cat - 1,
                            method = "REML",
                            data = hyp_1)
summary(hyp_1_categories)

# r2 hypothesis 1
r2_ml(hyp_1_categories)

# differences between males of distinct male quality 
significance.test(hyp_1_categories)


# experimental design influence on hypothesis 1 results ----

# general effect of physical contact
hyp_1_contact_general <- rma.mv(yi = smd,
                                V = var,
                                random = list(~1|effect_id,
                                              ~1|study_id,
                                              ~1|experiment_id,
                                              ~1|measure_id,
                                              ~1|species_id),
                                mods = ~ physical_contact - 1,
                                method = "REML",
                                data = hyp_1)
summary(hyp_1_contact_general)

# interaction between physical contact and male quality category
hyp_1_contact_interaction <- rma.mv(yi = smd,
                                    V = var,
                                    random = list(~1|effect_id,
                                                  ~1|study_id,
                                                  ~1|experiment_id,
                                                  ~1|measure_id,
                                                  ~1|species_id),
                                    mods = ~ paste(physical_contact,
                                                   m_quality_cat) - 1,
                                    method = "REML",
                                    data = hyp_1)
summary(hyp_1_contact_interaction)

# differences between males of distinct quality when females were allowed to interact with males
significance.test(hyp_1_contact_interaction)

# differences between males of distinct quality when females were precluded from interacting with males
significance.test(hyp_1_contact_interaction, option = 2)


# general effect of number of females available
hyp_1_f_available_general <- rma.mv(yi = smd,
                                    V = var,
                                    random = list(~1|effect_id,
                                                  ~1|study_id,
                                                  ~1|experiment_id,
                                                  ~1|measure_id,
                                                  ~1|species_id),
                                    mods = ~ f_available - 1,
                                    method = "REML",
                                    data = hyp_1)
summary(hyp_1_f_available_general)

# interaction between number of females available and male quality category
hyp_1_f_available_interaction <- rma.mv(yi = smd,
                                        V = var,
                                        random = list(~1|effect_id,
                                                      ~1|study_id,
                                                      ~1|experiment_id,
                                                      ~1|measure_id,
                                                      ~1|species_id),
                                        mods = ~ paste(f_available,
                                                       m_quality_cat) - 1,
                                        method = "REML",
                                        data = hyp_1)
summary(hyp_1_f_available_interaction)

# differences between males of distinct quality in multiple-choice tests
significance.test(hyp_1_f_available_interaction)

# differences between males of distinct quality in no-choice tests
significance.test(hyp_1_f_available_interaction, option = 2)

# general effect of female quality status
hyp_1_f_status_general <- rma.mv(yi = smd,
                                 V = var,
                                 random = list(~1|effect_id,
                                               ~1|study_id,
                                               ~1|experiment_id,
                                               ~1|measure_id,
                                               ~1|species_id),
                                 mods = ~ f_quality_status - 1,
                                 method = "REML",
                                 data = hyp_1)
summary(hyp_1_f_status_general)

# female quality status
hyp_1_f_status_interaction <- rma.mv(yi = smd,
                                 V = var,
                                 random = list(~1|effect_id,
                                               ~1|study_id,
                                               ~1|experiment_id,
                                               ~1|measure_id,
                                               ~1|species_id),
                                 mods = ~ paste(f_quality_status,
                                                m_quality_cat) - 1,
                                 method = "REML",
                                 data = hyp_1)
summary(hyp_1_f_status_interaction)

# differences between males of distinct quality when female quality is assumed
significance.test(hyp_1_f_status_interaction)

# differences between males of distinct quality when female quality is verified
significance.test(hyp_1_f_status_interaction, option = 2)

# r2 all experimental design moderators and their interaction with male quality category
all_exp_design <- rma.mv(yi = smd,
                         V = var,
                         random = list(~1|effect_id,
                                       ~1|study_id,
                                       ~1|experiment_id,
                                       ~1|measure_id,
                                       ~1|species_id),
                         mods = ~ m_quality_cat * physical_contact +
                           m_quality_cat * f_available +
                           m_quality_cat * f_quality_status,
                         method = "REML",
                         data = hyp_1)

r2_ml(all_exp_design)

# hypothesis 2 test ----

# hypothesis 2 general estimate
hyp_2_no_moderators <- rma.mv(yi = smd,
                        V = var,
                        random = list(~1|effect_id,
                                      ~1|study_id,
                                      ~1|experiment_id,
                                      ~1|measure_id,
                                      ~1|species_id),
                        method = "REML",
                        data = hyp_2)
summary(hyp_2_no_moderators)

# heterogeneity hypothesis 2
i2_ml(hyp_2_no_moderators, method = "wv") * 100

# interaction between male trait and male trait value category
set.seed(42)
hyp_2_categories <- rma.mv(yi = smd,
                           V = var,
                           random = list(~1|effect_id,
                                         ~1|study_id,
                                         ~1|experiment_id,
                                         ~1|measure_id,
                                         ~1|species_id),
                           mods = ~ paste(m_trait, m_trait_cat) - 1,
                           method = "REML",
                           data = hyp_2)
summary(hyp_2_categories)

# r2 hypothesis 2
r2_ml(hyp_2_categories)

# differences between males of different age
significance.test(hyp_2_categories, hypothesis = 2)

# differences between males of different body condition
significance.test(hyp_2_categories, hypothesis = 2, option = 2)

# differences between males of different body size (linear measures)
significance.test(hyp_2_categories, hypothesis = 2, option = 3)

# differences between males of different body mass
significance.test(hyp_2_categories, hypothesis = 2, option = 4)

# publication bias ----

# calculate squared-root of inverse of effective sample size
df_bias <- 
  mod_data %>%
  mutate(sqrt_inv_e_n = sqrt((n1 + n2)/(n1 * n2))) 

bias_model <- rma.mv(yi = smd,
                     V = var,
                     random = list(~1|effect_id,
                                   ~1|study_id,
                                   ~1|experiment_id,
                                   ~1|measure_id,
                                   ~1|species_id),
                     mods = ~ 1 +
                       scale(sqrt_inv_e_n) +
                       scale(year),
                     method = "REML",
                     data = df_bias)

summary(bias_model)
