pacman::p_load("coda", "tidyverse", "here", 
               "MCMCglmm", "brms", "MASS", "patchwork", 
               "phytools", "patchwork", "bayesplot", "tidybayes")


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
## fixed effect
metafor_p1 <- orchard_plot(HD_BM, 
                           group = "Study_id",
                           xlab = "Effect size (Hedges' d)",
                           angle = 45) + 
  scale_x_discrete(labels = c("Overall effect")) +
  scale_color_manual(values = "#CDBE70") +
  scale_fill_manual(values = "#EEDC82") + 
  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

## random effect
ci_var <- confint(HD_BM, level = 0.95) 
tbl <- tibble::as_tibble(ci_var, rownames = "term")

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

metafor_p2 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
  geom_point(size = 3.0, color = "#9AC0CD") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#9AC0CD") +
  coord_flip() +
  labs(x = NULL, y = "Variance 95% CI") +
  scale_y_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  theme_classic(base_size = 12)

metafor_eg1_1 <- metafor_p1 / metafor_p2

#### brms ----
get_variables(dat1_BM_brms)

## fixed effects
fixed_effects_samples_brms <- dat1_BM_brms %>%
  spread_draws(b_Intercept)
fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "#EEDC82", 
    color = "#CDBE70"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effect (variance)"
  ) +
  theme_classic()
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- dat1_BM_brms %>%
  spread_draws(sd_Phylo__Intercept, sd_Effect_id__Intercept, sd_Species__Intercept, sd_Study_id__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2,
         .variable = case_when(
           .variable == "sd_Phylo__Intercept"      ~ "Phylo",
           .variable == "sd_Study_id__Intercept"   ~ "Study",
           .variable == "sd_Species__Intercept"    ~ "Species",
           .variable == "sd_Effect_id__Intercept"  ~ "Effect_id"
         ),
         .variable = factor(.variable, 
                            levels = c("Phylo", "Study", "Species", "Effect_id")
                            )
         )

head(random_effects_samples_brms)

brms_p2 <- ggplot(random_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "#B2DFEE", 
    color = "#9AC0CD"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects (variance)"
  ) +
  theme_classic() 
brms_eg1_1 <- brms_p1/brms_p2

metafor_eg1_1 |brms_eg1_1

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
  geom_point(size = 3.0, color = "#9AC0CD") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#9AC0CD") +
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
    fill = "#EEDC82", 
    color = "#CDBE70"
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
  spread_draws(sd_Phylo__Intercept, sd_Effect_id__Intercept, sd_Species__Intercept, sd_Study_id__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2,
         .variable = case_when(
           .variable == "sd_Phylo__Intercept"      ~ "Phylo",
           .variable == "sd_Study_id__Intercept"   ~ "Study",
           .variable == "sd_Species__Intercept"    ~ "Species",
           .variable == "sd_Effect_id__Intercept"  ~ "Effect_id"
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
    fill = "#B2DFEE", 
    color = "#9AC0CD"
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
# use metafor_eg3_gau, metafor_eg3_exp, 
# m_exp, m_gau

# metafor
summary(metafor_eg3_gau)
s_metafor_p1 <- orchaRd::orchard_plot(metafor_eg3_gau, 
                           group = "site",
                           xlab = "Effect size",
                           angle = 45) + 
  scale_x_discrete(labels = c("Overall effect")) +
  scale_color_manual(values = "#CDBE70") +
  scale_fill_manual(values = "#EEDC82") + 
#  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

## random effect
ci_var <- confint(HD_BM, level = 0.95) 
tbl <- tibble::as_tibble(ci_var, rownames = "term")

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

metafor_p2 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
  geom_point(size = 3.0, color = "#9AC0CD") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#9AC0CD") +
  coord_flip() +
  labs(x = NULL, y = "Variance 95% CI") +
  scale_y_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  theme_classic(base_size = 12)

metafor_eg1_1 <- metafor_p1 / metafor_p2

## spatial e.g. 2 ----
# EXP_eg4_metafor, m_exp4_brms
dat_Roger <- read.csv(here("data", "examples", "Roger_etal_2024", "Roger_etal_2024.csv"))

summary(EXP_eg4_metafor)
s_metafor_p1 <- orchaRd::orchard_plot(EXP_eg4_metafor, 
                                      group = "effect_id",
                                      xlab = "Effect size",
                                      angle = 45) + 
  scale_x_discrete(labels = c("Overall effect")) +
  scale_color_manual(values = "#CDBE70") +
  scale_fill_manual(values = "#EEDC82") + 
  #  scale_y_continuous(breaks = seq(-4.0, 4.0, 1), limits = c(-4.0, 4.0)) + 
  theme_classic()

ci_var <- confint(EXP_eg4_metafor, level = 0.95) 
tbl <- tibble::as_tibble(ci_var, rownames = "term")

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

metafor_p2 <- ggplot(tbl_var, aes(x = label, y = estimate)) +
  geom_point(size = 3.0, color = "#9AC0CD") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, color = "#9AC0CD") +
  coord_flip() +
  labs(x = NULL, y = "Variance 95% CI") +
  scale_y_continuous(breaks = seq(0, 10.0, 1), limits = c(0, 10.0)) + 
  theme_classic(base_size = 12)

#### brms ----
summary(m_exp4_brms)

get_variables_dynamic <- function(model, pattern) {
  variables <- get_variables(model)
  variables[grep(pattern, variables)]
}

rename_vars_exp4 <- function(variable) {
  # fixed
  variable <- gsub("b_Intercept", "mu (overall mean)", variable)
  
  # random
  variable <- gsub("sd_effect_id__Intercept", "_effect (SD of effect_id)", variable)
  
  # GP hyperparameters
  variable <- gsub("sdgp_gpx_kmy_km", "SD_GP (spatial GP)", variable)
  variable <- gsub("lscale_gpx_kmy_km", "ℓ_GP (length-scale, km)", variable)
  
  # #
  # variable <- gsub("^sigma$", "σ_resid", variable)
  
  return(variable)
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
        fill           = "lightcyan3",
        color          = "lightcyan4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Fixed effects", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects_exp4: ", e$message)
    return(NULL)
  })
}

visualize_random_effects_exp4 <- function(model) {
  
  random_effect_vars <- get_variables_dynamic(model, "^sd_")
  # random_effect_vars <- setdiff(random_effect_vars, "sd…")
  
  if (length(random_effect_vars) == 0) {
    message("No random effects found")
    return(NULL)
  }
  
  tryCatch({
    random_effects_samples <- model %>%
      spread_draws(!!!syms(random_effect_vars)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(random_effect_vars),
        names_to  = ".variable",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.variable = rename_vars_exp4(.variable))
    
    ggplot(random_effects_samples, aes(x = .value, y = .variable)) +
      ggdist::stat_halfeye(
        normalize      = "xy",
        point_interval = "mean_qi",
        fill           = "olivedrab3",
        color          = "olivedrab4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Random effects (SD)", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_random_effects_exp4: ", e$message)
    return(NULL)
  })
}

visualize_gp_sd_exp4 <- function(model) {
  gp_sd_vars <- get_variables_dynamic(model, "^sdgp_")
  if (length(gp_sd_vars) == 0) {
    message("No GP SD parameters found")
    return(NULL)
  }
  
  tryCatch({
    gp_sd_samples <- model %>%
      spread_draws(!!!syms(gp_sd_vars)) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(gp_sd_vars),
        names_to  = ".variable",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.variable = rename_vars_exp4(.variable))
    
    ggplot(gp_sd_samples, aes(x = .value, y = .variable)) +
      ggdist::stat_halfeye(
        normalize      = "xy",
        point_interval = "mean_qi",
        fill           = "#FF6347",
        color          = "#8B3626"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "GP hyperparameters", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_gp_sd_exp4: ", e$message)
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
        fill           = "#FF6347",
        color          = "#8B3626"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "GP hyperparameters", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_gp_lscale_exp4: ", e$message)
    return(NULL)
  })
}

p_fix  <- visualize_fixed_effects_exp4(m_exp4_brms)
p_re   <- visualize_random_effects_exp4(m_exp4_brms)
p_gp_sd  <- visualize_gp_sd_exp4(m_exp4_brms)
p_gp_lsc <- visualize_gp_lscale_exp4(m_exp4_brms)

p_fix
p_re
p_gp_sd
p_gp_lsc
