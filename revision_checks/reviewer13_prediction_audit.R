#!/usr/bin/env Rscript

# Reviewer 13 audit. This script reads validated saved fits only: it does not
# fit, profile, or overwrite a model.

suppressPackageStartupMessages({
  library(metafor)
  library(posterior)
  library(readr)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript revision_checks/reviewer13_prediction_audit.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
out_dir <- file.path(root, "revision_checks", "reviewer13_prediction_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
q_summary <- function(x) {
  q <- stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  tibble(lower_95 = q[[1]], lower_50 = q[[2]], median = q[[3]],
         upper_50 = q[[4]], upper_95 = q[[5]])
}

baseline_path <- file.path(root, "revision_checks", "ou_correctness_outputs", "baseline_fit_objects.rds")
brms_path <- file.path(root, "Rdata", "tutorial_v2", "moura2021_BM_brms.rds")
assert(file.exists(baseline_path), paste("Missing validated BM object:", baseline_path))
assert(file.exists(brms_path), paste("Missing validated Moura brms RDS:", brms_path))

# This exact saved fit is the intercept-only BM example. Do not substitute the
# separately saved Moura meta-regression object, which has a moderator.
fit <- readRDS(baseline_path)$fit_bm
assert(inherits(fit, "rma.mv"), "baseline$fit_bm is not an rma.mv object.")
assert(fit$k == 1828L && fit$p == 1L && identical(fit$test, "t"),
       "Saved object is not the expected intercept-only Moura BM fit.")
assert(identical(unname(fit$s.names), c("study.id", "effect.size.id", "species.id", "species.id.phy")),
       "Unexpected random-effect ordering.")
assert(identical(as.character(fit$dfs), "residual") && length(fit$ddf) == 1L,
       "Unexpected degrees-of-freedom rule.")
assert("species.id.phy" %in% names(fit$R) &&
         isTRUE(all.equal(unname(diag(fit$R[["species.id.phy"]])), rep(1, 341L), tolerance = 1e-12)),
       "Phylogenetic covariance is not a unit-diagonal correlation matrix.")

# The target is a new study, new effect-size, and marginal new species
# realization. Each random-effect covariance has unit diagonal, so each fitted
# component enters once. No new vi is added: the target is latent, not an
# observed effect-size estimate.
component_names <- c("study", "effect_size", "species_nonphylogenetic", "species_phylogenetic")
component_variance <- stats::setNames(unname(fit$sigma2), component_names)
assert(all(is.finite(component_variance)) && all(component_variance >= 0), "Invalid fitted variances.")
mean_variance <- unname(fit$vb[1L, 1L])
heterogeneity_variance <- sum(component_variance)
prediction_variance <- mean_variance + heterogeneity_variance
df <- unname(fit$ddf[[1L]])
t_critical <- stats::qt(0.975, df)
pi_manual <- unname(fit$b[[1L]]) + c(-1, 1) * t_critical * sqrt(prediction_variance)

# Independent API check: metafor defaults newvi to zero and must agree.
pi_metafor <- metafor::predict.rma(fit)
assert(isTRUE(all.equal(as.numeric(pi_metafor$pi.lb), pi_manual[[1L]], tolerance = 1e-10)) &&
         isTRUE(all.equal(as.numeric(pi_metafor$pi.ub), pi_manual[[2L]], tolerance = 1e-10)),
       "Manual PI does not match metafor::predict.rma().")

component_table <- tibble(
  component = c("Pooled mean estimation", names(component_variance), "Total latent prediction variance"),
  variance = c(mean_variance, unname(component_variance), prediction_variance),
  included_in_target = TRUE,
  role = c("Variance of fitted intercept", rep("Marginal new random-effect realization", 4L),
           "Fitted-mean variance plus four marginal random-effect variances")
)
interval_table <- tibble(
  target = "new latent true effect; new study and marginal new species realization",
  scale = "Fisher's Z",
  estimate = unname(fit$b[[1L]]),
  confidence_interval_lower = unname(fit$ci.lb),
  confidence_interval_upper = unname(fit$ci.ub),
  prediction_interval_lower = pi_manual[[1L]],
  prediction_interval_upper = pi_manual[[2L]],
  pooled_mean_variance = mean_variance,
  study_variance = component_variance[["study"]],
  effect_size_variance = component_variance[["effect_size"]],
  species_nonphylogenetic_variance = component_variance[["species_nonphylogenetic"]],
  species_phylogenetic_variance = component_variance[["species_phylogenetic"]],
  heterogeneity_variance = heterogeneity_variance,
  latent_prediction_variance = prediction_variance,
  degrees_of_freedom = df,
  t_critical = t_critical,
  future_sampling_variance = 0,
  variance_component_uncertainty_propagated = FALSE
)
write_csv(component_table, file.path(out_dir, "frequentist_latent_prediction_components.csv"))
write_csv(interval_table, file.path(out_dir, "frequentist_latent_prediction_interval.csv"))

# Preserve posterior pairing by using the full brms RDS. The precomputed long
# CSV has marginal values only and cannot support this draw-wise construction.
fit_brms <- readRDS(brms_path)
assert(inherits(fit_brms, "brmsfit"), "Saved Moura object is not a brmsfit.")
draws <- posterior::as_draws_df(fit_brms)
required <- c("b_Intercept", "sd_study.id__Intercept", "sd_effect.size.id__Intercept",
              "sd_species.id__Intercept", "sd_species.id.phy__Intercept")
assert(all(required %in% names(draws)), "Validated brms RDS lacks a required draw variable.")
set.seed(20260912L)
latent_draws <- tibble(
  draw = seq_len(nrow(draws)),
  pooled_mean = draws$b_Intercept,
  new_study = stats::rnorm(nrow(draws), 0, draws$sd_study.id__Intercept),
  new_effect_size = stats::rnorm(nrow(draws), 0, draws$sd_effect.size.id__Intercept),
  new_species_nonphylogenetic = stats::rnorm(nrow(draws), 0, draws$sd_species.id__Intercept),
  new_species_phylogenetic = stats::rnorm(nrow(draws), 0, draws$sd_species.id.phy__Intercept)
) |>
  mutate(latent_true_effect = pooled_mean + new_study + new_effect_size +
           new_species_nonphylogenetic + new_species_phylogenetic)
posterior_summary <- q_summary(latent_draws$latent_true_effect) |>
  mutate(target = "new latent true effect; new study and marginal new species realization",
         scale = "Fisher's Z", draws = nrow(latent_draws), seed = 20260912L,
         future_sampling_error_included = FALSE)
write_csv(latent_draws, file.path(out_dir, "brms_latent_prediction_draws.csv"))
write_csv(posterior_summary, file.path(out_dir, "brms_latent_prediction_summary.csv"))

metadata <- tibble(
  item = c("frequentist_fit", "bayesian_fit", "legacy_thin_line_function", "legacy_thin_line_package",
           "legacy_thin_line_target", "future_sampling_error", "posterior_simulation_seed"),
  value = c("baseline_fit_objects.rds$fit_bm", "moura2021_BM_brms.rds",
            "orchaRd::pred_interval_esmeans() via orchard_plot()", "orchaRd 2.2.1",
            "same marginal latent target for this sigma2-only BM model",
            "excluded because no future sampling variance is specified", "20260912")
)
write_csv(metadata, file.path(out_dir, "audit_metadata.csv"))
writeLines(c(
  "# Reviewer 13 prediction audit",
  "",
  "Target: the latent true effect for a new effect-size observation from a new study and a marginal new species-level realization. It includes fitted study, effect-size, non-phylogenetic species, and phylogenetic species heterogeneity, plus uncertainty in the pooled mean. It excludes future sampling error.",
  "",
  "The frequentist interval is a plug-in t prediction interval conditional on estimated variance components. It does not separately propagate uncertainty in those component estimates.",
  "",
  "The legacy Figure 3 thin line came from orchaRd::orchard_plot() through pred_interval_esmeans(), not a direct metafor call. For this simple BM model it agrees with metafor::predict.rma() with default newvi = 0.",
  "",
  "The Bayesian distribution was constructed explicitly from joint posterior draws and four new marginal random-effect draws. brms::posterior_predict() was not used because it simulates observed outcomes with sampling error."
), file.path(out_dir, "README.md"))

cat("REVIEWER13_PREDICTION_AUDIT_CHECKS_PASSED\n")
