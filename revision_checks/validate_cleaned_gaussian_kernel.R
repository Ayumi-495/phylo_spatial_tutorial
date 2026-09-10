# Focused validator for the cleaned Gaussian audit. It does not refit models.
# Usage: Rscript revision_checks/validate_cleaned_gaussian_kernel.R <inputs|fits|identification|i2|handoff>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("inputs", "fits", "identification", "i2", "handoff")) {
  stop("Usage: Rscript revision_checks/validate_cleaned_gaussian_kernel.R <inputs|fits|identification|i2|handoff>")
}
stage <- args[[1L]]
root <- normalizePath(".")
input_dir <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")
out_dir <- file.path(root, "revision_checks", "cleaned_gaussian_kernel_outputs")

read_one <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

if (stage == "inputs") {
  p <- readRDS(file.path(input_dir, "published_cleaned_prepared.rds"))
  d <- p$dat; D <- p$distance_km
  record <- read_one(file.path(out_dir, "cleaned_gaussian_input_validation.csv"))
  stopifnot(
    nrow(d) == 2355L, nlevels(d$effect_id) == 2355L, nlevels(d$study_id) == 390L,
    nlevels(d$site_id) == 380L, nlevels(d$const) == 1L,
    identical(rownames(D), levels(d$site_id)), identical(colnames(D), levels(d$site_id)),
    nrow(record) == 1L, isTRUE(record$same_yi_V_X_as_saved_cleaned_SPEXP[[1L]]),
    isTRUE(record$distance_rows_match_site_levels[[1L]]), isTRUE(record$distance_cols_match_site_levels[[1L]])
  )
  cat("CLEANED_GAUSSIAN_INPUTS_VALIDATED\n")
}

if (stage == "fits") {
  sp <- read_one(file.path(out_dir, "spatial_only_spgau_result.csv"))
  multi <- read_one(file.path(out_dir, "combined_spgau_targeted_multistart.csv"))
  stopifnot(
    nrow(sp) == 1L, nrow(multi) == 3L,
    identical(sort(multi$label), sort(c("short_312_km", "intermediate_800_km", "long_3000_km"))),
    all(c(sp$kernel, multi$kernel) == "SPGAU"),
    all(c(sp$n_effects, multi$n_effects) == 2355L),
    all(c(sp$n_studies, multi$n_studies) == 390L),
    all(c(sp$n_sites, multi$n_sites) == 380L),
    all(is.finite(c(sp$pooled_mean, sp$ci_lb, sp$ci_ub, sp$iid_effect_variance, sp$spatial_variance, sp$rho_km, sp$logLik_REML, sp$AIC_REML))),
    all(is.finite(as.matrix(multi[, c("pooled_mean", "ci_lb", "ci_ub", "iid_effect_variance", "study_variance", "spatial_variance", "rho_km", "logLik_REML", "AIC_REML")]))),
    all(!grepl("not_converged|error", c(sp$convergence_status, multi$convergence_status), ignore.case = TRUE))
  )
  cat("CLEANED_GAUSSIAN_FITS_VALIDATED\n")
}

if (stage == "identification") {
  multi <- read_one(file.path(out_dir, "combined_spgau_targeted_multistart.csv"))
  zero <- read_one(file.path(out_dir, "combined_spgau_tau2_zero_result.csv"))
  summary <- read_one(file.path(out_dir, "cleaned_gaussian_identification_summary.csv"))
  old <- readRDS(file.path(input_dir, "unstructured_only.rds"))
  old_ll <- as.numeric(old$fit.stats["ll", "REML"])
  stopifnot(
    nrow(zero) == 1L, zero$fixed_spatial_variance[[1L]] == 0,
    is.finite(zero$logLik_REML[[1L]]), is.finite(summary$combined_logLik_spread[[1L]]),
    summary$combined_logLik_spread[[1L]] >= 0,
    abs(zero$logLik_REML[[1L]] - old_ll) < 1e-6,
    abs(summary$tau_zero_vs_saved_unstructured_logLik_difference[[1L]]) < 1e-6,
    summary$best_vs_tau2_zero_logLik_loss[[1L]] >= 0
  )
  write.csv(data.frame(
    best_combined_label = summary$best_combined_label[[1L]],
    combined_logLik_spread = summary$combined_logLik_spread[[1L]],
    rho_min_km = min(multi$rho_km), rho_max_km = max(multi$rho_km),
    best_vs_tau2_zero_logLik_loss = summary$best_vs_tau2_zero_logLik_loss[[1L]],
    tau_zero_equals_unstructured_within = abs(zero$logLik_REML[[1L]] - old_ll),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "cleaned_gaussian_validation_summary.csv"), row.names = FALSE)
  cat("CLEANED_GAUSSIAN_IDENTIFICATION_VALIDATED\n")
}

if (stage == "i2") {
  x <- read_one(file.path(out_dir, "cleaned_gaussian_generalized_i2.csv"))
  expected_v <- 0.110937235977794
  stopifnot(all(abs(x$v_tilde - expected_v) < 1e-12))
  split_rows <- split(x, interaction(x$model, x$solution, drop = TRUE))
  for (z in split_rows) {
    total <- z[z$component == "total", , drop = FALSE]
    parts <- z[z$component != "total", , drop = FALSE]
    stopifnot(nrow(total) == 1L, abs(total$variance - sum(parts$variance)) < 1e-10,
              abs(total$I2_percent - sum(parts$I2_percent)) < 1e-8)
  }
  cat("CLEANED_GAUSSIAN_I2_VALIDATED\n")
}

if (stage == "handoff") {
  handoff <- file.path(root, "revision_checks", "SPATIAL_AUDIT_HANDOFF_2026-09-09.md")
  stopifnot(file.exists(handoff), file.exists(file.path(root, "revision_checks", "cleaned_gaussian_kernel_audit_2026-09-10.md")))
  text <- paste(readLines(handoff, warn = FALSE), collapse = "\n")
  new_paths <- c(
    "revision_checks/cleaned_gaussian_kernel_audit_2026-09-10.md",
    "revision_checks/cleaned_gaussian_kernel_audit.R",
    "revision_checks/cleaned_gaussian_kernel_outputs/cleaned_gaussian_input_validation.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/spatial_only_spgau_result.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/combined_spgau_targeted_multistart.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/combined_spgau_tau2_zero_result.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/cleaned_gaussian_generalized_i2.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/cleaned_gaussian_identification_summary.csv",
    "revision_checks/cleaned_gaussian_kernel_outputs/cleaned_gaussian_validation_summary.csv"
  )
  stopifnot(all(file.exists(file.path(root, new_paths))), all(vapply(new_paths, grepl, logical(1), x = text, fixed = TRUE)))
  forbidden <- c("tutorial_v2.qmd", "manuscript", "response letter")
  changed <- system2("git", c("diff", "--name-only", "9eea010..HEAD"), stdout = TRUE)
  # Before the commit, compare both index and worktree to the explicit checkpoint.
  changed <- unique(c(changed, system2("git", c("diff", "--name-only", "9eea010"), stdout = TRUE),
                      system2("git", c("diff", "--cached", "--name-only", "9eea010"), stdout = TRUE)))
  stopifnot(!any(basename(changed) %in% forbidden), all(!nzchar(changed) | startsWith(changed, "revision_checks/")))
  cat("CLEANED_GAUSSIAN_HANDOFF_VALIDATED\n")
}
