pacman::p_load("coda", "tidyverse", "here", "metafor",
               "MCMCglmm", "brms", "MASS", "patchwork", 
               "phytools", "patchwork", "bayesplot", "tidybayes", "orchaRd")


# Figure 1 ----
# grid of (normalized) distances
d <- seq(0, 1, length.out = 501)

# parameters
rho <- 0.2          # spatial range parameter
alpha <- 1 / rho    # OU–exponential equivalence: alpha = 1/rho
bm_floor <- 1e-3    # tiny floor so BM-like curve doesn't hit zero (illustrative)

# kernel functions (illustrative BM; exact BM depends on shared path lengths on a tree)
k_linear      <- pmax(0, 1 - d / rho)                   # spatial linear with hard cutoff at d = rho
k_exponential <- exp(-d / rho)                          # spatial exponential
k_squaredexp  <- exp(-(d^2) / (rho^2))                  # spatial squared exponential (Gaussian)
k_ou          <- exp(-alpha * d)                        # phylogenetic OU; identical to exponential with alpha = 1/rho
k_bm_like     <- pmax(bm_floor, 1 - d)                  # BM-like linear decay without hard cutoff (illustrative)

# assemble tidy data
df <- tibble(
  distance = rep(d, 5),
  correlation = c(k_bm_like, k_ou, k_linear, k_exponential, k_squaredexp),
  model = factor(rep(c("BM",
                       "OU (α = 1/ρ)",
                       "Linear (spatial)",
                       "Exponential (spatial)",
                       "Squared exponential (spatial)"),
                     each = length(d)),
                 levels = c("BM",
                            "OU (α = 1/ρ)",
                            "Linear (spatial)",
                            "Exponential (spatial)",
                            "Squared exponential (spatial)"))
)

# plot
ggplot(df, aes(distance, correlation, color = model)) +
  geom_line(size = 1) +
  labs(x = "Distance (normalised)", y = "Correlation",
       color = NULL,
       title = "Comparison of correlation–distance relationships") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")


# Worked example figure ----

## phylo eg. 1 ----
### intercept-only model ----

#### metafor ----
summary(moura_BM_metafor1)
## fixed effect
metafor_p1 <- orchaRd::orchard_plot(moura_BM_metafor1, 
                           group = "study.id",
                           xlab = "Effect size",
                           angle = 45) + 
  scale_x_discrete(labels = c("Overall effect")) +
  scale_color_manual(values = "#CDAD00") +
  scale_fill_manual(values = "#FFD700") + 
  # scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

## random effect
ci_var <- confint(moura_BM_metafor1, level = 0.95) 
# estimate  ci.lb  ci.ub 
# sigma^2.1   0.0192 0.0108 0.0325 
# sigma.1     0.1384 0.1038 0.1802 
# 
# estimate  ci.lb  ci.ub 
# sigma^2.2   0.0145 0.0121 0.0172 
# sigma.2     0.1202 0.1099 0.1311 
# 
# estimate  ci.lb  ci.ub 
# sigma^2.3   0.0557 0.0334 0.0788 
# sigma.3     0.2359 0.1827 0.2807 
# 
# estimate  ci.lb  ci.ub 
# sigma^2.4   0.0512 0.0179 0.1792 
# sigma.4     0.2263 0.1336 0.4233 

tbl <- tibble::as_tibble(ci_var, rownames = "term")

label_map <- c("Effect_id", "Study_id", "Species", "Phylo")
tbl_var <- tbl %>%
  filter(str_detect(term, "^sigma\\^2\\.")) %>%
  mutate(
    idx   = as.integer(str_match(term, "\\.(\\d+)$")[,2]),
    label = case_when(
      idx == 1 ~ "Study_id",
      idx == 2 ~ "Effect_id",
      idx == 3 ~ "Species",
      idx == 4 ~ "Phylo"
    ),
    label = factor(label, levels = label_map) 
  )

metafor_p2 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
  geom_point(size = 3.0, color = "#528B8B") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#528B8B") +
  coord_flip() +
  labs(x = NULL, y = "Variance 95% CI") +
  # scale_y_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  theme_classic(base_size = 12)

metafor_eg1_1 <- metafor_p1 / metafor_p2
metafor_eg1_1

#### brms ----
moura_BM_brms1 <- readRDS(here("Rdata", "moura2021_BM_brms.rds"))
get_variables(moura_BM_brms1)

## fixed effects
fixed_effects_samples_brms <- moura_BM_brms1 %>%
  spread_draws(b_Intercept)
fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill           = "#FFD700",
    color          = "#CDAD00"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-1.0, 2.0, 1)) + 
  labs(y = "Overall effect"
  ) +
  theme_classic()
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- moura_BM_brms1 %>%
  spread_draws(sd_species.id.phy__Intercept, sd_effect.size.id__Intercept, sd_species.id__Intercept, sd_study.id__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2,
         .variable = case_when(
           .variable == "sd_species.id.phy__Intercept" ~ "Phylo",
           .variable == "sd_study.id__Intercept" ~ "Study",
           .variable == "sd_species.id__Intercept" ~ "Species",
           .variable == "sd_effect.size.id__Intercept" ~ "Effect_id"
         ),
         .variable = factor(.variable, 
                            levels = c("Effect_id", "Study", "Species", "Phylo")
                            )
         )

head(random_effects_samples_brms)

brms_p2 <- ggplot(random_effects_samples_brms, aes(x = .value^2, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "#98F5FF",
    color = "#528B8B"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
　scale_x_continuous(breaks = seq(0, 0.08, 0.02), limits = c(0, 0.08)) + 
  labs(
    # title = "Posterior distributions of random effects - brms",
    y = "Random effects (variance)"
  ) +
  theme_classic() 
brms_eg1_1 <- brms_p1 /brms_p2
brms_eg1_1

### meta-regression ----
#### metafor ----
## fixed effect
res1 <- orchaRd::mod_results(HD_BM_2, mod = "Feeding_mode",  group = "Study_id", subset = TRUE)
res1

metafor_p3 <- orchard_plot(res1, 
                           mod = "Feeding_mode",
                           group = "Study_id",
                           xlab = "Effect size (Hedges' d)",
                           angle = 45) + 
  # scale_x_discrete(labels = c("Overall effect")) +
  # scale_color_manual(values = "#CDBE70") +
  # scale_fill_manual(values = "#EEDC82") + 
  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

## random effect
ci_var2 <- confint(HD_B_2, level = 0.95) 
tbl <- tibble::as_tibble(ci_var2, rownames = "term")

label_map <- c("Phylo", "Study_id", "Species", "Effect_id")
tbl_var <- tbl %>%
  filter(str_detect(term, "^sigma\\^2\\.")) %>%
  mutate(
    idx   = as.integer(str_match(term, "\\.(\\d+)$")[,2]),
    label = case_when(
      idx == 1 ~ "Study_id",
      idx == 2 ~ "Effect_id",
      idx == 3 ~ "Species",
      idx == 4 ~ "Phylo"
    ),
    label = factor(label, levels = label_map) 
  )

metafor_p4 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
  geom_point(size = 3.0, color = "#528B8B") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#528B8B") +
  coord_flip() +
  labs(x = NULL, y = "Variance 95% CI") +
  scale_y_continuous(breaks = seq(0, 11.0, 1), limits = c(0, 11.0)) + 
  theme_classic(base_size = 12)

metafor_eg1_2 <- metafor_p3 / metafor_p4

#### brms ----
get_variables(dat1_BM_brms2)

## fixed effects
fixed_effects_samples_brms <- dat1_BM_brms2 %>%
  spread_draws(b_Feeding_modeEndophytic, b_Feeding_modeFreeMliving)
fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p3 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill           = "#FFD700",
    color          = "#CDAD00"
    ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effect (variance)"
  ) +
  theme_classic()
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- dat1_BM_brms2 %>%
  spread_draws(sd_species.id.phy__Intercept, sd_effect.size.id__Intercept, sd_species.id.phy__Intercept, sd_study.id__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2,
         .variable = case_when(
           .variable == "sd_species.id.phy__Intercept" ~ "Phylo",
           .variable == "sd_study.id__Intercept" ~ "Study",
           .variable == "sd_species.id.phy__Intercept" ~ "Species",
           .variable == "sd_effect.size.id__Intercept" ~ "Effect_id"
         ),
         .variable = factor(.variable, 
                            levels = c("Phylo", "Study", "Species", "Effect_id")
         )
  )

head(random_effects_samples_brms)

brms_p4 <- ggplot(random_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "#79CDCD",
    color = "#528B8B"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects (variance)"
  ) +
  theme_classic() 
brms_eg1_2 <- brms_p3/brms_p4

metafor_eg1_2 |brms_eg1_2

## spatial e.g. 1 ----
# # use metafor_eg3_gau, metafor_eg3_exp, 
# # m_exp, m_gau
# 
# # metafor
# summary(metafor_eg3_gau)
# s_metafor_p1 <- orchaRd::orchard_plot(metafor_eg3_gau, 
#                            group = "site",
#                            xlab = "Effect size",
#                            angle = 45) + 
#   scale_x_discrete(labels = c("Overall effect")) +
#   scale_color_manual(values = "#CDBE70") +
#   scale_fill_manual(values = "#EEDC82") + 
# #  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
#   theme_classic()
# 
# ## random effect
# ci_var <- confint(HD_BM, level = 0.95) 
# tbl <- tibble::as_tibble(ci_var, rownames = "term")
# 
# label_map <- c("Phylo", "Study_id", "Species", "Effect_id")
# tbl_var <- tbl %>%
#   filter(str_detect(term, "^sigma\\^2\\.")) %>%
#   mutate(
#     idx   = as.integer(str_match(term, "\\.(\\d+)$")[,2]),
#     label = case_when(
#       idx == 1 ~ "Study_id",
#       idx == 2 ~ "Effect_id",
#       idx == 3 ~ "Species",
#       idx == 4 ~ "Phylo"
#     ),
#     label = factor(label, levels = label_map) 
#   )
# 
# metafor_p2 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
#   geom_point(size = 3.0, color = "#9AC0CD") +
#   geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#9AC0CD") +
#   coord_flip() +
#   labs(x = NULL, y = "Variance 95% CI") +
#   scale_y_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
#   theme_classic(base_size = 12)
# 
# metafor_eg1_1 <- metafor_p1 / metafor_p2

## spatial e.g. 2 ----
dat_Roger <- read.csv(here("data", "examples", "Roger_etal_2024", "Roger_etal_2024.csv"))
sp_eg1_exp_mf <- readRDS(here("Rdata", "tutorial_v2", "sp_eg1_exp_mf.rds"))
### metafor ----
summary(sp_eg1_exp_mf)
s_metafor_p1 <- orchaRd::orchard_plot(sp_eg1_exp_mf, 
                                      group = "study_id",
                                      xlab = "Effect size",
                                      angle = 45,
                                      twig.size = 0.3,
                                      trunk.size = 0.5
) + 
  scale_x_discrete(labels = c("Overall effect")) +
  scale_color_manual(values = "#FFD700") +
  scale_fill_manual(values = "#CDAD00") + 
  #  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

ci_var2 <- confint(sp_eg1_exp_mf, level = 0.95)
tbl2 <- tibble::as_tibble(ci_var2, rownames = "term")
# A tibble: 5 × 4
# term    estimate  ci.lb ci.ub
# <chr>      <dbl>  <dbl> <dbl>
#   1 sigma^2    0.793 0.726  0.868
# 2 sigma      0.891 0.852  0.932
# 3 tau^2      1.23  1.01   1.50 
# 4 tau        1.11  1.01   1.22 
# 5 rho        0.174 0.0174 0.857

tbl_var2 <- tbl2 %>% 
  filter(term %in% c("tau^2","sigma^2")) %>% 
  mutate(term = factor(term, levels = c("tau^2","sigma^2")))

s_metafor_p1var <- ggplot(tbl_var2, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "#528B8B") +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), height = 0.25, color = "#528B8B") +
  # scale_x_continuous(limits = c(0.0, 1.80)) + 
  labs(x = "Variance with with 95% CI", y = NULL) +
  theme_classic(base_size = 12)

tbl_rho2 <- tbl2 %>% 
  filter(term == "rho")

s_metafor_p1rho <- ggplot(tbl_rho2, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "#CD6090") +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), 
                height = 0.25, color = "#CD6090") +
  labs(x = "rho with 95% CI", y = NULL) +
  theme_classic(base_size = 12)

#### brms ----

summary(sp_eg1_exp_brms)

get_variables_dynamic <- function(model, pattern) {
  variables <- get_variables(model)
  variables[grep(pattern, variables)]
}

get_variables(sp_eg1_exp_brms)

rename_vars_exp4 <- function(variable) {
  variable <- as.character(variable)

  # fixed
  variable <- gsub("^b_Intercept$", "Overall effect", variable)

  # residual SD (sigma in brms)
  variable <- gsub("^sigma$", "Residual SD", variable)

  # GP SD
  variable <- gsub("^sdgp_gpx_kmy_km$", "SD of gp_x", variable)

  # GP length-scale
  variable <- gsub("^lscale_gpx_kmy_km$", "lscale", variable)

  variable
}

### fixed effect ----
visualize_fixed_effects_exp4 <- function(model) {
  fixed_effect_vars <- get_variables_dynamic(model, "^b_")
  if (length(fixed_effect_vars) == 0) {
    message("No fixed effects found")
    return(NULL)
  }
  
  tryCatch({
    fixed_effects_samples <- model %>%
      spread_draws(!!!syms(fixed_effect_vars)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(fixed_effect_vars),
        names_to  = ".variable",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.variable = rename_vars_exp4(.variable))
    
    ggplot(fixed_effects_samples, aes(x = .value, y = .variable)) +
      ggdist::stat_halfeye(
        normalize      = "xy",
        point_interval = "mean_qi",
        fill           = "#FFD700",
        color          = "#CDAD00"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Fixed effects", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects_exp4: ", e$message)
    return(NULL)
  })
}

### random effect ----
sd_to_var <- function(df) {
  df %>%
    dplyr::mutate(
      .value    = (.value)^2,
      .variable = gsub("^SD of ", "Var of ", .variable)
    )
}

#### random effect (Var) ----
visualize_random_and_gp_var_exp4 <- function(model) {
  
  re_vars <- get_variables_dynamic(model, "^(sigma|sdgp_)")
  
  if (length(re_vars) == 0) {
    message("No random or GP SD parameters found")
    return(NULL)
  }
  
  tryCatch({
    samples <- model %>%
      spread_draws(!!!syms(re_vars)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(re_vars),
        names_to  = ".variable",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.variable = rename_vars_exp4(.variable)) %>%
      sd_to_var()
    
    ggplot(samples, aes(x = .value, y = .variable)) +
      ggdist::stat_halfeye(
        normalize      = "xy",
        point_interval = "mean_qi",
        fill           = "#B0E2FF",
        color          = "#4F94CD"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(
        y = "Random and GP effects (variance)",
        x = "Posterior values (variance)"
      ) +
      # scale_x_continuous(limits = c(0.7, 2.0)) + 
      theme_classic()
    
  }, error = function(e) {
    message("Error in visualize_random_and_gp_var_exp4: ", e$message)
    return(NULL)
  })
}

visualize_gp_lscale_exp4 <- function(model) {
  gp_ls_vars <- get_variables_dynamic(model, "^lscale_")
  if (length(gp_ls_vars) == 0) {
    message("No GP length-scale parameters found")
    return(NULL)
  }
  
  tryCatch({
    gp_ls_samples <- model %>%
      spread_draws(!!!syms(gp_ls_vars)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(gp_ls_vars),
        names_to  = ".variable",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.variable = rename_vars_exp4(.variable))
    
    ggplot(gp_ls_samples, aes(x = .value, y = .variable)) +
      ggdist::stat_halfeye(
        normalize      = "xy",
        point_interval = "mean_qi",
        fill           = "#CD6090",
        color          = "#8B3A62"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "lscale", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_gp_lscale_exp4: ", e$message)
    return(NULL)
  })
}

s_brms_p1  <- visualize_fixed_effects_exp4(sp_eg1_exp_brms)
s_metafor_p1sd <- visualize_random_and_gp_var_exp4(sp_eg1_exp_brms)
s_metafor_p1lscale <- visualize_gp_lscale_exp4(sp_eg1_exp_brms)


metafor_sp2.1 <- s_metafor_p1 / (s_metafor_p1var + s_metafor_p1rho)
metafor_sp2.1
brms_sp2.1 <- s_brms_p1 / (s_metafor_p1sd + s_metafor_p1lscale)
brms_sp2.1
