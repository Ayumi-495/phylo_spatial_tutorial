#!/usr/bin/env Rscript

# Recreate static tutorial figures without running any model fit. Moura and Lim
# values are the validated frequentist results printed in tutorial_v2.qmd.
# Spain values are read from saved regional-audit output. This script does not
# fit brms, metafor, or glmmTMB models and does not alter saved RDS objects.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript revision_checks/restore_legacy_tutorial_figures.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
figure_dir <- file.path(root, "figs", "tutorial")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

result_figure <- function(fixed, components, title, filename, component_note) {
  p_fixed <- ggplot(fixed, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.16,
                   colour = "#CDAD00", na.rm = TRUE) +
    geom_point(size = 3, colour = "#CDAD00") +
    labs(title = title, x = "Fixed effect (95% confidence interval)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))

  p_components <- ggplot(components, aes(x = estimate, y = component)) +
    geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.16,
                   colour = "#528B8B", na.rm = TRUE) +
    geom_point(size = 3, colour = "#528B8B") +
    facet_grid(component ~ ., scales = "free_x", space = "free_y", switch = "y") +
    labs(x = component_note, y = NULL) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 0),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())

  ggsave(file.path(figure_dir, filename), p_fixed / p_components,
         width = 8.5, height = 7.2, dpi = 220, bg = "white")
}

with_profile_ci <- function(component, estimate, ci_lb, ci_ub) {
  tibble(component = component, estimate = estimate, ci_lb = ci_lb, ci_ub = ci_ub)
}

point_only <- function(component, estimate) {
  tibble(component = component, estimate = estimate, ci_lb = NA_real_, ci_ub = NA_real_)
}

# Moura et al. BM meta-analysis (metafor): source values and profile CIs shown
# in the metafor section of tutorial_v2.qmd.
result_figure(
  fixed = tibble(term = "Pooled mean", estimate = 0.3682, ci_lb = 0.1131, ci_ub = 0.6232),
  components = bind_rows(
    with_profile_ci("Study", 0.0192, 0.0108, 0.0325),
    with_profile_ci("Effect size", 0.0145, 0.0121, 0.0172),
    with_profile_ci("Species (non-phylogenetic)", 0.0557, 0.0334, 0.0788),
    with_profile_ci("Species (phylogenetic)", 0.0512, 0.0179, 0.1792)
  ),
  title = "Moura et al.: BM meta-analysis (metafor)",
  filename = "metafor_moura2021.png",
  component_note = "Variance (component-specific 95% profile confidence interval)"
)

# Moura et al. BM meta-analysis (glmmTMB): source values displayed in the
# glmmTMB section. The figure uses its reported Wald CI for the pooled mean;
# variance components are fitted point estimates.
result_figure(
  fixed = tibble(term = "Pooled mean", estimate = 0.3681658, ci_lb = 0.1132610, ci_ub = 0.6230707),
  components = bind_rows(
    point_only("Study", 0.1384142^2),
    point_only("Effect size", 0.01445014),
    point_only("Species (non-phylogenetic)", 0.2359271^2),
    point_only("Species (phylogenetic)", 0.05122353)
  ),
  title = "Moura et al.: BM meta-analysis (glmmTMB)",
  filename = "tmb_eg1_1.png",
  component_note = "Variance (fitted point estimate; component-specific horizontal scale)"
)

# Lim et al. BM meta-regression (metafor): source values printed in the model
# section. Component CIs were not retained in the tutorial output, so the
# component panel deliberately shows point estimates only.
result_figure(
  fixed = tibble(
    term = c("Intercept", "Wild-environment contrast"),
    estimate = c(-0.1395, 0.0166), ci_lb = c(-0.3921, -0.1460), ci_ub = c(0.1131, 0.1792)
  ),
  components = bind_rows(
    point_only("Effect size", 0.0630),
    point_only("Species (phylogenetic)", 0.0521),
    point_only("Species (non-phylogenetic)", 0.0792)
  ),
  title = "Lim et al.: BM meta-regression (metafor)",
  filename = "metafor_lim.png",
  component_note = "Variance (fitted point estimate; component-specific horizontal scale)"
)

# Lim et al. BM meta-regression (glmmTMB): fixed-effect point estimates and
# standard errors are printed in tutorial_v2.qmd. The interval is the usual
# estimate +/- 1.96 standard errors; components are fitted point estimates.
lim_tmb_fixed <- tibble(
  term = c("Intercept", "Wild-environment contrast"),
  estimate = c(-0.13951, 0.01661), se = c(0.12893, 0.08469)
) |>
  mutate(ci_lb = estimate - 1.96 * se, ci_ub = estimate + 1.96 * se)
result_figure(
  fixed = select(lim_tmb_fixed, term, estimate, ci_lb, ci_ub),
  components = bind_rows(
    point_only("Effect size", 0.0630428),
    point_only("Species (phylogenetic)", 0.0521124),
    point_only("Species (non-phylogenetic)", 0.0792370)
  ),
  title = "Lim et al.: BM meta-regression (glmmTMB)",
  filename = "tmb_eg2_mr.png",
  component_note = "Variance (fitted point estimate; component-specific horizontal scale)"
)

# Spain spatial-only metafor figure: use the saved regional-audit result and
# the saved profile-likelihood CI files. No model fitting or profiling occurs.
spain_dir <- file.path(root, "revision_checks", "regional_cross_package_audit_outputs")
spain <- read.csv(file.path(spain_dir, "metafor_spatial_only_result.csv"), check.names = FALSE)
tau_ci <- read.csv(file.path(spain_dir, "metafor_profile_ci_tau2.csv"), row.names = 1, check.names = FALSE)[1, ]
rho_ci <- read.csv(file.path(spain_dir, "metafor_profile_ci_rho.csv"), row.names = 1, check.names = FALSE)[1, ]

p_spain_mean <- ggplot(tibble(term = "Pooled mean", estimate = spain$mean,
                              ci_lb = spain$ci_lb, ci_ub = spain$ci_ub),
                       aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.16, colour = "#CDAD00") +
  geom_point(size = 3, colour = "#CDAD00") +
  labs(title = "Spain regional subset: spatial-only model (metafor)",
       x = "Pooled effect (95% t confidence interval)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12))
p_spain_variance <- ggplot(bind_rows(
  point_only("IID effect-size variance", spain$iid_effect_variance),
  with_profile_ci("Spatial variance", tau_ci$estimate, tau_ci$ci.lb, tau_ci$ci.ub)
), aes(x = estimate, y = component)) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.16, colour = "#528B8B", na.rm = TRUE) +
  geom_point(size = 3, colour = "#528B8B") +
  labs(x = "Variance (spatial: 95% profile confidence interval)", y = NULL) +
  theme_classic(base_size = 12)
p_spain_range <- ggplot(tibble(term = "Exponential range", estimate = rho_ci$estimate,
                               ci_lb = rho_ci$ci.lb, ci_ub = rho_ci$ci.ub), aes(x = estimate, y = term)) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.16, colour = "#528B8B") +
  geom_point(size = 3, colour = "#528B8B") +
  labs(x = "Range in km (95% profile confidence interval)", y = NULL) +
  theme_classic(base_size = 12)
ggsave(file.path(figure_dir, "spain_metafor_spatial_only.png"),
       p_spain_mean / p_spain_variance / p_spain_range,
       width = 8.5, height = 8.2, dpi = 220, bg = "white")

writeLines(c(
  "Static tutorial figures regenerated without fitting models.",
  paste("R:", R.version.string),
  "Moura/Lim values: existing validated frequentist output printed in tutorial_v2.qmd.",
  "Spain values: saved regional metafor result and saved tau2/rho profile CIs.",
  "No brms, metafor, or glmmTMB model was fitted or modified."
), file.path(root, "revision_checks", "legacy_tutorial_figure_restore_manifest.txt"))

message("Wrote tutorial figures to: ", figure_dir)
