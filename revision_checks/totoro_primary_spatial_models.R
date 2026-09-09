# Primary full-dataset spatial meta-analysis audit for Totoro.
# Each completed model is saved immediately. This script does not modify the tutorial qmd.

suppressPackageStartupMessages({
  library(metafor)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript totoro_primary_spatial_models.R <data_csv> <output_dir>")
}
data_csv <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(out_dir))

capture_result <- function(expr) {
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
    return(if (identical(fit$optres$convergence, 0L)) "converged" else "not_converged")
  }
  if (any(grepl("converg", warnings, ignore.case = TRUE))) return("optimizer_warning")
  "completed_no_explicit_optimizer_status"
}

fit_stat <- function(fit, statistic) {
  as.numeric(fit$fit.stats[statistic, "REML"])
}

append_result <- function(row) {
  path <- file.path(out_dir, "primary_model_results.csv")
  old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  new <- rbind(old, row)
  write.csv(new, path, row.names = FALSE)
}

save_fit <- function(name, result, elapsed_seconds) {
  if (inherits(result$value, "error")) {
    row <- data.frame(
      model = name, mean = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
      sigma2_effect = NA_real_, sigma2_study = NA_real_, tau2_spatial = NA_real_, rho_km = NA_real_,
      logLik_REML = NA_real_, AIC_REML = NA_real_, convergence_status = "error",
      elapsed_seconds = elapsed_seconds,
      warnings = paste(c(conditionMessage(result$value), result$warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )
    append_result(row)
    return(invisible(NULL))
  }

  fit <- result$value
  saveRDS(fit, file.path(out_dir, paste0(name, ".rds")))
  sigma2 <- fit$sigma2
  row <- data.frame(
    model = name,
    mean = as.numeric(fit$b[1]),
    ci_lb = as.numeric(fit$ci.lb[1]),
    ci_ub = as.numeric(fit$ci.ub[1]),
    sigma2_effect = if (length(sigma2) >= 1L) sigma2[1] else NA_real_,
    sigma2_study = if (length(sigma2) >= 2L) sigma2[2] else NA_real_,
    tau2_spatial = if (length(fit$tau2)) fit$tau2[1] else NA_real_,
    rho_km = if (length(fit$rho)) fit$rho[1] else NA_real_,
    logLik_REML = fit_stat(fit, "ll"),
    AIC_REML = fit_stat(fit, "AIC"),
    convergence_status = fit_status(fit, result$warnings),
    elapsed_seconds = elapsed_seconds,
    warnings = paste(result$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
  append_result(row)
  invisible(fit)
}

dat <- read.csv(data_csv, stringsAsFactors = FALSE)
dat <- dat[complete.cases(dat$latitude, dat$longitude), ]
stopifnot(nrow(dat) == 2361L, all(dat$var_Hedges > 0))

# A site is a unique recorded coordinate pair. It is not an assertion about field-site identity.
dat$effect_id <- factor(seq_len(nrow(dat)))
dat$study_id <- factor(dat$study_id)
dat$site_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
site_levels <- sort(unique(dat$site_key))
dat$site_id <- factor(dat$site_key, levels = site_levels)
dat$const <- factor("all_sites")
stopifnot(nlevels(dat$study_id) == 393L, nlevels(dat$site_id) == 383L)

site_lookup <- unique(dat[c("site_id", "site_key", "latitude", "longitude")])
site_lookup <- site_lookup[match(site_levels, site_lookup$site_key), ]
stopifnot(identical(as.character(site_lookup$site_id), site_levels))
stopifnot(identical(levels(dat$site_id), site_levels))

# WGS84 ellipsoidal great-circle distances: geosphere::distGeo() returns metres.
coordinates_lonlat <- as.matrix(site_lookup[c("longitude", "latitude")])
distance_km <- geosphere::distm(coordinates_lonlat, fun = geosphere::distGeo) / 1000
rownames(distance_km) <- colnames(distance_km) <- site_levels

# Required ID/matrix assertions before model fitting.
stopifnot(identical(rownames(distance_km), levels(dat$site_id)))
stopifnot(identical(colnames(distance_km), levels(dat$site_id)))
stopifnot(isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)))
stopifnot(all(abs(diag(distance_km)) < 1e-10), all(distance_km >= 0), all(is.finite(distance_km)))

write.csv(site_lookup, file.path(out_dir, "site_lookup.csv"), row.names = FALSE)
write.csv(distance_km, file.path(out_dir, "great_circle_distance_km.csv"), row.names = TRUE)
writeLines(c(
  paste("R:", R.version.string),
  paste("metafor:", as.character(packageVersion("metafor"))),
  paste("geosphere:", as.character(packageVersion("geosphere"))),
  "Distance method: geosphere::distm(..., fun = geosphere::distGeo) / 1000; WGS84 ellipsoidal great-circle kilometres.",
  "Sampling-error assumption: V = var_Hedges, hence diagonal and independent.",
  "Spatial outer group: const has one level, so spatial covariance is allowed across studies."
), file.path(out_dir, "primary_audit_metadata.txt"))

run_and_save <- function(name, expr) {
  started <- proc.time()[["elapsed"]]
  result <- capture_result(expr)
  elapsed <- proc.time()[["elapsed"]] - started
  save_fit(name, result, elapsed)
  message(name, " completed in ", round(elapsed, 1), " seconds")
}

# Same yi, diagonal V, intercept-only fixed effect, and effect-level heterogeneity in all models.
run_and_save("unstructured_only", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~ 1 | effect_id, ~ 1 | study_id),
  data = dat, method = "REML", test = "t", sparse = TRUE
))

run_and_save("spatial_only", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~ 1 | effect_id, ~ site_id | const),
  struct = "SPEXP", dist = list(site_id = distance_km),
  data = dat, method = "REML", test = "t", sparse = TRUE
))

run_and_save("combined", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~ 1 | effect_id, ~ 1 | study_id, ~ site_id | const),
  struct = "SPEXP", dist = list(site_id = distance_km),
  data = dat, method = "REML", test = "t", sparse = TRUE
))
