#!/usr/bin/env Rscript

# Reproducible audit of the saved profile-likelihood confidence intervals for
# the final Moura intercept-only BM model. The authoritative tutorial source
# retains the direct metafor::confint(phylo_eg1_meta_ma_BM) output (lines
# 943--958 and 1643--1658 of tutorial_v2.qmd). Re-running its default profile
# search is unnecessary and computationally expensive for k = 1,828; this
# script verifies the saved estimates against the final fitted model and saves
# the retained, rounded printed profile results as a tabular audit artifact.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript revision_checks/reviewer13_moura_variance_profiles.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
out_dir <- file.path(root, "revision_checks", "reviewer13_prediction_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

fit <- readRDS(file.path(root, "revision_checks", "ou_correctness_outputs", "baseline_fit_objects.rds"))$fit_bm
assert(inherits(fit, "rma.mv") && fit$k == 1828L && fit$p == 1L && length(fit$sigma2) == 4L,
       "Expected saved final intercept-only Moura BM fit with four sigma2 components.")

# Rounded values copied from the retained final-model metafor::confint() output
# cited above, not full-precision bounds, a Wald calculation, or a new model fit.
profile_ci <- tribble(
  ~component_index, ~component, ~estimate, ~ci_lb, ~ci_ub,
  1L, "Study variance", 0.0192, 0.0108, 0.0325,
  2L, "Effect-size variance", 0.0145, 0.0121, 0.0172,
  3L, "Species variance, non-phylogenetic", 0.0557, 0.0334, 0.0788,
  4L, "Species variance, phylogenetic", 0.0512, 0.0179, 0.1792
) |>
  mutate(
    saved_fit_estimate = unname(fit$sigma2[component_index]),
    estimate_agrees_with_saved_fit = abs(estimate - saved_fit_estimate) < 5e-4,
    interval_method = "retained final-model metafor profile-likelihood 95% CI (rounded printed confint bounds)",
    source = "tutorial_v2.qmd lines 943-958 and 1643-1658",
    fit_changed = FALSE
  )

assert(all(profile_ci$estimate_agrees_with_saved_fit),
       "A retained profile estimate no longer agrees with the saved final BM fit.")
assert(all(is.finite(profile_ci$ci_lb)) && all(is.finite(profile_ci$ci_ub)) &&
         all(profile_ci$ci_lb >= 0) && all(profile_ci$ci_ub > profile_ci$ci_lb),
       "Retained profile intervals are not finite, non-negative, and ordered.")

write_csv(profile_ci, file.path(out_dir, "moura_bm_variance_profile_ci.csv"))
saveRDS(profile_ci, file.path(out_dir, "moura_bm_variance_profile_ci.rds"))
cat("REVIEWER13_MOURA_SAVED_PROFILE_AUDIT_PASSED\n")
