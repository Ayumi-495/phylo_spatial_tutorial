# Focused SPGAU sensitivity on the validated published_cleaned Grau-Andres data.
#
# This deliberately reuses the saved cleaned prepared object; it does not rebuild
# the R18 exclusion or coordinate pipeline. It fits exactly one spatial-only
# model and three targeted combined starts (short/intermediate/~3000 km). The
# tau2 = 0 restricted likelihood is recovered exactly from the saved cleaned
# unstructured fit, which is mathematically the same model. No profile grid or
# additional spatial model is run.
#
# Usage from repository root:
#   Rscript revision_checks/cleaned_gaussian_kernel_audit.R fit

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || args[[1L]] != "fit") {
  stop("Usage: Rscript revision_checks/cleaned_gaussian_kernel_audit.R fit")
}

root <- normalizePath(".")
input_dir <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")
out_dir <- file.path(root, "revision_checks", "cleaned_gaussian_kernel_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

capture_conditions <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(force(expr), error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

fit_status <- function(fit, warnings) {
  if (!is.null(fit$converged)) return(if (isTRUE(fit$converged)) "converged" else "not_converged")
  if (!is.null(fit$optres$convergence)) {
    return(if (identical(as.integer(fit$optres$convergence), 0L)) "converged" else "not_converged")
  }
  if (any(grepl("converg|Hessian|optim", warnings, ignore.case = TRUE))) return("optimizer_warning")
  "completed_no_explicit_optimizer_status"
}

prepared_path <- file.path(input_dir, "published_cleaned_prepared.rds")
if (!file.exists(prepared_path)) stop("Missing validated cleaned object: ", prepared_path)
prepared <- readRDS(prepared_path)
dat <- prepared$dat
distance_km <- prepared$distance_km

stopifnot(
  nrow(dat) == 2355L,
  is.factor(dat$effect_id), nlevels(dat$effect_id) == 2355L,
  is.factor(dat$study_id), nlevels(dat$study_id) == 390L,
  is.factor(dat$site_id), nlevels(dat$site_id) == 380L,
  is.factor(dat$const), nlevels(dat$const) == 1L,
  all(is.finite(dat$d_Hedges)), all(is.finite(dat$var_Hedges)), all(dat$var_Hedges > 0),
  identical(rownames(distance_km), levels(dat$site_id)),
  identical(colnames(distance_km), levels(dat$site_id)),
  identical(rownames(distance_km), colnames(distance_km)),
  isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
  all(is.finite(distance_km)), all(distance_km >= 0), all(abs(diag(distance_km)) < 1e-10)
)

# Assert the exact same cleaned yi, V, and X used by the saved primary SPEXP fits.
saved_spexp <- lapply(c("unstructured_only", "spatial_only", "combined"), function(x) {
  path <- file.path(input_dir, paste0(x, ".rds"))
  if (!file.exists(path)) stop("Missing saved cleaned SPEXP fit: ", path)
  readRDS(path)
})
names(saved_spexp) <- c("unstructured", "spatial_only", "combined")
for (fit in saved_spexp) {
  stopifnot(
    length(fit$yi) == nrow(dat),
    isTRUE(all.equal(as.numeric(fit$yi), as.numeric(dat$d_Hedges), tolerance = 0)),
    isTRUE(all.equal(as.numeric(fit$vi), as.numeric(dat$var_Hedges), tolerance = 0)),
    ncol(fit$X) == 1L, isTRUE(all.equal(as.numeric(fit$X[, 1L]), rep(1, nrow(dat)), tolerance = 0))
  )
}

input_record <- data.frame(
  n_effects = nrow(dat), n_studies = nlevels(dat$study_id), n_sites = nlevels(dat$site_id),
  n_const_levels = nlevels(dat$const), distance_units = "WGS84 ellipsoidal geodesic kilometres",
  distance_rows_match_site_levels = identical(rownames(distance_km), levels(dat$site_id)),
  distance_cols_match_site_levels = identical(colnames(distance_km), levels(dat$site_id)),
  same_yi_V_X_as_saved_cleaned_SPEXP = TRUE,
  stringsAsFactors = FALSE
)
write.csv(input_record, file.path(out_dir, "cleaned_gaussian_input_validation.csv"), row.names = FALSE)

fit_spgau <- function(random, sigma2.init = NULL, tau2.init = NULL, rho.init = NULL,
                      tau2_fixed = NULL, rho_fixed = NULL) {
  # The same assertion is deliberately repeated immediately before every fit.
  stopifnot(
    identical(rownames(distance_km), levels(dat$site_id)),
    identical(colnames(distance_km), levels(dat$site_id))
  )
  call <- list(
    yi = dat$d_Hedges, V = dat$var_Hedges, mods = ~ 1,
    random = random, struct = "SPGAU", dist = list(site_id = distance_km),
    data = dat, method = "REML", test = "t", sparse = TRUE
  )
  if (!is.null(sigma2.init) || !is.null(tau2.init) || !is.null(rho.init)) {
    call$control <- list(
      sigma2.init = sigma2.init,
      tau2.init = tau2.init,
      rho.init = rho.init
    )
  }
  if (!is.null(tau2_fixed)) call$tau2 <- tau2_fixed
  if (!is.null(rho_fixed)) call$rho <- rho_fixed
  do.call(metafor::rma.mv, call)
}

result_row <- function(fit, label, model, elapsed, warnings,
                       start_effect_variance = NA_real_, start_study_variance = NA_real_,
                       start_spatial_variance = NA_real_, start_rho_km = NA_real_,
                       fixed_spatial_variance = NA_real_, fixed_rho_km = NA_real_) {
  data.frame(
    model = model, label = label, kernel = "SPGAU",
    n_effects = nrow(dat), n_studies = nlevels(dat$study_id), n_sites = nlevels(dat$site_id),
    distance_units = "km", correlation = "exp(-d^2/rho^2)",
    start_effect_variance = start_effect_variance, start_study_variance = start_study_variance,
    start_spatial_variance = start_spatial_variance, start_rho_km = start_rho_km,
    fixed_spatial_variance = fixed_spatial_variance, fixed_rho_km = fixed_rho_km,
    pooled_mean = as.numeric(fit$b[1L]), ci_lb = as.numeric(fit$ci.lb[1L]), ci_ub = as.numeric(fit$ci.ub[1L]),
    iid_effect_variance = as.numeric(fit$sigma2[1L]),
    study_variance = if (length(fit$sigma2) >= 2L) as.numeric(fit$sigma2[2L]) else NA_real_,
    spatial_variance = as.numeric(fit$tau2[1L]), rho_km = as.numeric(fit$rho[1L]),
    logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]), AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
    convergence_status = fit_status(fit, warnings), elapsed_seconds = elapsed,
    warnings = paste(warnings, collapse = " | "), stringsAsFactors = FALSE
  )
}

fit_and_save <- function(label, model, random, ..., path_stem) {
  started <- proc.time()[["elapsed"]]
  captured <- capture_conditions(fit_spgau(random, ...))
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(captured$value, "error")) stop(label, " failed: ", conditionMessage(captured$value))
  fit <- captured$value
  saveRDS(fit, file.path(out_dir, paste0(path_stem, ".rds")))
  list(fit = fit, elapsed = elapsed, warnings = captured$warnings)
}

spatial_random <- list(~ 1 | effect_id, ~ site_id | const)
combined_random <- list(~ 1 | effect_id, ~ 1 | study_id, ~ site_id | const)

# One primary spatial-only Gaussian fit.
spatial_only <- fit_and_save("spatial-only", "spatial_only", spatial_random,
                             path_stem = "spatial_only_spgau")
spatial_row <- result_row(spatial_only$fit, "primary", "spatial_only", spatial_only$elapsed,
                          spatial_only$warnings)
write.csv(spatial_row, file.path(out_dir, "spatial_only_spgau_result.csv"), row.names = FALSE)

# Exactly three informative combined starts: short, intermediate, and the prior ~3000-km branch.
starts <- list(
  short_312_km = list(sigma2 = c(0.751, 1.140), tau2 = 0.080, rho = 312),
  intermediate_800_km = list(sigma2 = c(0.751, 1.160), tau2 = 0.060, rho = 800),
  long_3000_km = list(sigma2 = c(0.751, 1.210), tau2 = 0.040, rho = 3000)
)
combined_rows <- list()
for (label in names(starts)) {
  s <- starts[[label]]
  item <- fit_and_save(
    label, "combined", combined_random,
    sigma2.init = s$sigma2, tau2.init = s$tau2, rho.init = s$rho,
    path_stem = paste0("combined_spgau_", label)
  )
  combined_rows[[label]] <- result_row(
    item$fit, label, "combined", item$elapsed, item$warnings,
    start_effect_variance = s$sigma2[1L], start_study_variance = s$sigma2[2L],
    start_spatial_variance = s$tau2, start_rho_km = s$rho
  )
}
combined_rows <- do.call(rbind, combined_rows)
combined_rows$delta_logLik_from_best <- max(combined_rows$logLik_REML) - combined_rows$logLik_REML
combined_rows$delta_AIC_from_best <- combined_rows$AIC_REML - min(combined_rows$AIC_REML)
write.csv(combined_rows, file.path(out_dir, "combined_spgau_targeted_multistart.csv"), row.names = FALSE)

best_index <- which.max(combined_rows$logLik_REML)[1L]
best_label <- combined_rows$label[[best_index]]
writeLines(best_label, file.path(out_dir, "combined_spgau_best_label.txt"))
best_file <- file.path(out_dir, paste0("combined_spgau_", best_label, ".rds"))
file.copy(best_file, file.path(out_dir, "combined_spgau_best.rds"), overwrite = TRUE)

# A tau2=0 Gaussian combined model is exactly the saved cleaned unstructured
# model: rho is nonidentified when its variance is zero. Reusing its likelihood
# avoids an unnecessary numerical optimization of a degenerate spatial term.
unstructured <- saved_spexp$unstructured
tau_zero_row <- data.frame(
  model = "combined_tau2_zero", label = "tau2_fixed_zero", kernel = "SPGAU",
  n_effects = nrow(dat), n_studies = nlevels(dat$study_id), n_sites = nlevels(dat$site_id),
  distance_units = "km", correlation = "exp(-d^2/rho^2)",
  start_effect_variance = NA_real_, start_study_variance = NA_real_,
  start_spatial_variance = NA_real_, start_rho_km = NA_real_,
  fixed_spatial_variance = 0, fixed_rho_km = NA_real_,
  pooled_mean = as.numeric(unstructured$b[1L]), ci_lb = as.numeric(unstructured$ci.lb[1L]),
  ci_ub = as.numeric(unstructured$ci.ub[1L]),
  iid_effect_variance = as.numeric(unstructured$sigma2[1L]),
  study_variance = as.numeric(unstructured$sigma2[2L]), spatial_variance = 0,
  rho_km = NA_real_,
  logLik_REML = as.numeric(unstructured$fit.stats["ll", "REML"]),
  AIC_REML = as.numeric(unstructured$fit.stats["AIC", "REML"]),
  convergence_status = "reused_exact_unstructured_restriction", elapsed_seconds = 0,
  warnings = "No refit: tau2=0 makes rho unidentified and yields the saved cleaned unstructured model.",
  stringsAsFactors = FALSE
)
write.csv(tau_zero_row, file.path(out_dir, "combined_spgau_tau2_zero_result.csv"), row.names = FALSE)

# The clean generalized-I2 v_tilde is prevalidated by reviewer18_cleaned_i2.R.
v_tilde <- 0.110937235977794
i2_rows <- rbind(
  data.frame(model = "spatial_only", solution = "primary", component = c("iid_effect", "spatial", "total"),
             variance = c(spatial_row$iid_effect_variance, spatial_row$spatial_variance,
                          spatial_row$iid_effect_variance + spatial_row$spatial_variance)),
  data.frame(model = "combined", solution = combined_rows$label[best_index], component = c("iid_effect", "study", "spatial", "total"),
             variance = c(combined_rows$iid_effect_variance[best_index], combined_rows$study_variance[best_index],
                          combined_rows$spatial_variance[best_index], sum(combined_rows[best_index, c("iid_effect_variance", "study_variance", "spatial_variance")])))
)
i2_rows$v_tilde <- v_tilde
i2_rows$I2_percent <- 100 * i2_rows$variance / (i2_rows$variance + v_tilde)
total_indices <- which(i2_rows$component == "total")
for (idx in total_indices) {
  model_rows <- i2_rows$model == i2_rows$model[idx] & i2_rows$solution == i2_rows$solution[idx]
  components <- i2_rows$component[model_rows] != "total"
  denom <- i2_rows$variance[idx] + v_tilde
  i2_rows$I2_percent[which(model_rows)[components]] <- 100 * i2_rows$variance[which(model_rows)[components]] / denom
  i2_rows$I2_percent[idx] <- 100 * i2_rows$variance[idx] / denom
}
write.csv(i2_rows, file.path(out_dir, "cleaned_gaussian_generalized_i2.csv"), row.names = FALSE)

summary_record <- data.frame(
  best_combined_label = best_label,
  combined_logLik_spread = max(combined_rows$logLik_REML) - min(combined_rows$logLik_REML),
  best_vs_tau2_zero_logLik_loss = combined_rows$logLik_REML[best_index] - tau_zero_row$logLik_REML,
  tau_zero_vs_saved_unstructured_logLik_difference = tau_zero_row$logLik_REML - as.numeric(saved_spexp$unstructured$fit.stats["ll", "REML"]),
  stringsAsFactors = FALSE
)
write.csv(summary_record, file.path(out_dir, "cleaned_gaussian_identification_summary.csv"), row.names = FALSE)

message("CLEANED_GAUSSIAN_FITS_SAVED")
