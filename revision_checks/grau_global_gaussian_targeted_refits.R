# Limited targeted free-refit audit for the full Grau-Andres combined SPGAU model.
#
# This script deliberately tests only three informative starts:
#   1. the optimized fixed-rho=200 km solution;
#   2. an intermediate rho=800 km start;
#   3. the existing primary rho~3091 km solution.
# It is not a broad multi-start search.
#
# Usage:
#   Rscript grau_global_gaussian_targeted_refits.R fixed200 <data.csv> <output_dir>
#   Rscript grau_global_gaussian_targeted_refits.R free <data.csv> <output_dir> <fixed200|intermediate800|primary3091>
#   Rscript grau_global_gaussian_targeted_refits.R compile <data.csv> <output_dir>

suppressPackageStartupMessages({
  library(metafor)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) stop("See usage in file header.")
stage <- args[[1L]]
data_csv <- args[[2L]]
out_dir <- args[[3L]]
start_label <- if (length(args) == 4L) args[[4L]] else NA_character_
stopifnot(stage %in% c("fixed200", "free", "compile"))
if (stage == "free") stopifnot(start_label %in% c("fixed200", "intermediate800", "primary3091"))

target_dir <- file.path(out_dir, "targeted_free_refits")
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

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

prepare_data <- function() {
  dat <- read.csv(data_csv, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat$latitude, dat$longitude), , drop = FALSE]
  stopifnot(nrow(dat) == 2361L, all(is.finite(dat$d_Hedges)),
            all(is.finite(dat$var_Hedges)), all(dat$var_Hedges > 0))
  dat$effect_id <- factor(seq_len(nrow(dat)))
  dat$study_id <- factor(dat$study_id)
  dat$site_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
  site_levels <- sort(unique(dat$site_key))
  dat$site_id <- factor(dat$site_key, levels = site_levels)
  dat$const <- factor("all_sites")
  stopifnot(nlevels(dat$study_id) == 393L, nlevels(dat$site_id) == 383L,
            nlevels(dat$const) == 1L)

  sites <- unique(dat[c("site_key", "latitude", "longitude")])
  sites <- sites[match(site_levels, sites$site_key), , drop = FALSE]
  distance_km <- geosphere::distm(as.matrix(sites[c("longitude", "latitude")]),
                                  fun = geosphere::distGeo) / 1000
  rownames(distance_km) <- colnames(distance_km) <- site_levels
  stopifnot(identical(rownames(distance_km), levels(dat$site_id)),
            identical(colnames(distance_km), levels(dat$site_id)),
            isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
            all(abs(diag(distance_km)) < 1e-10), all(distance_km >= 0),
            all(is.finite(distance_km)))
  list(dat = dat, distance_km = distance_km)
}

fit_combined <- function(prepared, sigma2.init, tau2.init, rho.init, rho_fixed = NULL) {
  dat <- prepared$dat
  distance_km <- prepared$distance_km
  stopifnot(identical(rownames(distance_km), levels(dat$site_id)),
            identical(colnames(distance_km), levels(dat$site_id)))
  call <- list(
    yi = dat$d_Hedges, V = dat$var_Hedges, mods = ~ 1,
    random = list(~ 1 | effect_id, ~ 1 | study_id, ~ site_id | const),
    struct = "SPGAU", dist = list(site_id = distance_km),
    data = dat, method = "REML", test = "t", sparse = TRUE,
    control = list(sigma2.init = sigma2.init,
                   tau2.init = tau2.init,
                   rho.init = rho.init)
  )
  if (!is.null(rho_fixed)) call$rho <- rho_fixed
  do.call(metafor::rma.mv, call)
}

result_row <- function(fit, label, start_sigma2, start_tau2, start_rho,
                       elapsed, warnings, fixed_rho = NA_real_) {
  data.frame(
    start_label = label,
    fixed_rho_km = fixed_rho,
    start_effect_variance = start_sigma2[1],
    start_study_variance = start_sigma2[2],
    start_spatial_variance = start_tau2,
    start_rho_km = start_rho,
    pooled_mean = as.numeric(fit$b[1]),
    ci_lb = as.numeric(fit$ci.lb[1]), ci_ub = as.numeric(fit$ci.ub[1]),
    iid_effect_variance = as.numeric(fit$sigma2[1]),
    study_variance = as.numeric(fit$sigma2[2]),
    spatial_variance = as.numeric(fit$tau2[1]),
    rho_km = as.numeric(fit$rho[1]),
    logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]),
    AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
    convergence_status = fit_status(fit, warnings),
    elapsed_seconds = elapsed,
    warnings = paste(warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

primary_path <- file.path(out_dir, "combined_spgau.rds")
if (!file.exists(primary_path)) stop("Missing primary fit: ", primary_path)
primary <- readRDS(primary_path)
stopifnot(inherits(primary, "rma.mv"), length(primary$sigma2) == 2L,
          length(primary$tau2) == 1L, length(primary$rho) == 1L)

if (stage == "fixed200") {
  p <- prepare_data()
  started <- proc.time()[["elapsed"]]
  captured <- capture_conditions(fit_combined(
    p, sigma2.init = primary$sigma2, tau2.init = primary$tau2,
    rho.init = 200, rho_fixed = 200
  ))
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(captured$value, "error")) stop(captured$value)
  fit <- captured$value
  saveRDS(fit, file.path(target_dir, "fixed_rho_200.rds"))
  row <- result_row(fit, "fixed_rho_200", primary$sigma2, primary$tau2,
                    200, elapsed, captured$warnings, fixed_rho = 200)
  write.csv(row, file.path(target_dir, "fixed_rho_200.csv"), row.names = FALSE)
  message("FIXED_RHO_200_SAVED")
}

if (stage == "free") {
  p <- prepare_data()
  fixed_path <- file.path(target_dir, "fixed_rho_200.rds")
  if (!file.exists(fixed_path)) stop("Run fixed200 first: missing ", fixed_path)
  fixed200 <- readRDS(fixed_path)
  starts <- switch(start_label,
    fixed200 = list(sigma2 = fixed200$sigma2, tau2 = fixed200$tau2, rho = 200),
    intermediate800 = list(
      sigma2 = (fixed200$sigma2 + primary$sigma2) / 2,
      tau2 = (fixed200$tau2 + primary$tau2) / 2,
      rho = 800
    ),
    primary3091 = list(sigma2 = primary$sigma2, tau2 = primary$tau2, rho = primary$rho)
  )
  started <- proc.time()[["elapsed"]]
  captured <- capture_conditions(fit_combined(
    p, sigma2.init = starts$sigma2, tau2.init = starts$tau2,
    rho.init = starts$rho
  ))
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(captured$value, "error")) stop(captured$value)
  fit <- captured$value
  saveRDS(fit, file.path(target_dir, paste0("free_", start_label, ".rds")))
  row <- result_row(fit, start_label, starts$sigma2, starts$tau2,
                    starts$rho, elapsed, captured$warnings)
  write.csv(row, file.path(target_dir, paste0("free_", start_label, ".csv")), row.names = FALSE)
  message("FREE_REFIT_SAVED: ", start_label)
}

if (stage == "compile") {
  paths <- file.path(target_dir, paste0("free_",
    c("fixed200", "intermediate800", "primary3091"), ".csv"))
  if (!all(file.exists(paths))) stop("Missing free-refit CSVs: ",
                                     paste(paths[!file.exists(paths)], collapse = ", "))
  rows <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
  rows$delta_logLik_from_best <- max(rows$logLik_REML) - rows$logLik_REML
  rows$delta_AIC_from_best <- rows$AIC_REML - min(rows$AIC_REML)
  write.csv(rows, file.path(target_dir, "targeted_free_refits_compiled.csv"), row.names = FALSE)
  message("TARGETED_FREE_REFITS_COMPILED")
}
