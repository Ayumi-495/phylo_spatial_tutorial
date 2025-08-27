pacman::p_load("coda", "tidyverse", "here", 
               "MCMCglmm", "brms", "MASS", "patchwork", 
               "phytools", "patchwork", "bayesplot", "tidybayes")



# worked example figure ----

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
