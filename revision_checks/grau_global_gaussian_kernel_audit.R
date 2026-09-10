# Full global Grau-Andres Gaussian-kernel sensitivity audit.
# This script does not modify tutorial_v2.qmd. It uses WGS84 ellipsoidal
# great-circle distances in kilometres and metafor's exact SPGAU correlation:
# Cor(d) = exp(-d^2 / rho^2). Thus rho is the e-folding distance.

suppressPackageStartupMessages({
  library(metafor)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L || !args[[1L]] %in% c("fit", "profile")) {
  stop("Usage: Rscript grau_global_gaussian_kernel_audit.R <fit|profile> <data_csv> <output_dir> [spatial_only|combined|both]")
}
stage <- args[[1L]]
data_csv <- args[[2L]]
out_dir <- args[[3L]]
profile_models <- if (length(args) == 4L) args[[4L]] else "both"
stopifnot(profile_models %in% c("spatial_only", "combined", "both"))
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
  if (!is.null(fit$optres$convergence)) return(if (identical(fit$optres$convergence, 0L)) "converged" else "not_converged")
  if (any(grepl("converg", warnings, ignore.case = TRUE))) return("optimizer_warning")
  "completed_no_explicit_optimizer_status"
}

prepare_data <- function() {
  dat <- read.csv(data_csv, stringsAsFactors = FALSE)
  dat <- dat[complete.cases(dat$latitude, dat$longitude), , drop = FALSE]
  stopifnot(nrow(dat) == 2361L, all(is.finite(dat$d_Hedges)), all(dat$var_Hedges > 0))
  dat$effect_id <- factor(seq_len(nrow(dat)))
  dat$study_id <- factor(dat$study_id)
  dat$site_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
  site_levels <- sort(unique(dat$site_key))
  dat$site_id <- factor(dat$site_key, levels = site_levels)
  dat$const <- factor("all_sites")
  stopifnot(nlevels(dat$study_id) == 393L, nlevels(dat$site_id) == 383L)

  sites <- unique(dat[c("site_key", "latitude", "longitude")])
  sites <- sites[match(site_levels, sites$site_key), , drop = FALSE]
  stopifnot(identical(sites$site_key, site_levels))
  distance_km <- geosphere::distm(as.matrix(sites[c("longitude", "latitude")]),
                                  fun = geosphere::distGeo) / 1000
  rownames(distance_km) <- colnames(distance_km) <- site_levels
  stopifnot(identical(rownames(distance_km), levels(dat$site_id)),
            identical(colnames(distance_km), levels(dat$site_id)),
            isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
            all(abs(diag(distance_km)) < 1e-10), all(distance_km >= 0))
  list(dat = dat, distance_km = distance_km, sites = sites)
}

save_result <- function(name, captured, elapsed) {
  if (inherits(captured$value, "error")) {
    stop("Fit failed for ", name, ": ", conditionMessage(captured$value))
  }
  fit <- captured$value
  saveRDS(fit, file.path(out_dir, paste0(name, ".rds")))
  row <- data.frame(
    model = name, kernel = "Gaussian: exp(-d^2/rho^2)",
    mean = as.numeric(fit$b[1]), ci_lb = as.numeric(fit$ci.lb[1]), ci_ub = as.numeric(fit$ci.ub[1]),
    iid_effect_variance = as.numeric(fit$sigma2[1]),
    study_variance = if (length(fit$sigma2) > 1L) as.numeric(fit$sigma2[2]) else NA_real_,
    spatial_variance = as.numeric(fit$tau2[1]),
    e_folding_range_km = as.numeric(fit$rho[1]),
    logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]),
    AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
    convergence_status = fit_status(fit, captured$warnings), elapsed_seconds = elapsed,
    warnings = paste(captured$warnings, collapse = " | "), stringsAsFactors = FALSE
  )
  write.csv(row, file.path(out_dir, paste0(name, "_result.csv")), row.names = FALSE)
  fit
}

if (stage == "fit") {
  p <- prepare_data()
  write.csv(p$sites, file.path(out_dir, "site_lookup.csv"), row.names = FALSE)
  write.csv(p$distance_km, file.path(out_dir, "great_circle_distance_km.csv"), row.names = TRUE)
  writeLines(c(
    paste("R:", R.version.string), paste("metafor:", packageVersion("metafor")),
    paste("geosphere:", packageVersion("geosphere")),
    "Distance: WGS84 ellipsoidal great-circle kilometres via geosphere::distGeo.",
    "Spatial Gaussian correlation: exp(-d^2/rho^2); rho is the e-folding distance.",
    "Same yi, diagonal V, intercept-only fixed effect, iid effect-level term, and constant spatial outer group as the audited SPEXP fits."
  ), file.path(out_dir, "metadata.txt"))

  started <- proc.time()[["elapsed"]]
  spatial <- capture_conditions(rma.mv(
    yi = d_Hedges, V = var_Hedges,
    random = list(~ 1 | effect_id, ~ site_id | const),
    struct = "SPGAU", dist = list(site_id = p$distance_km),
    data = p$dat, method = "REML", test = "t", sparse = TRUE
  ))
  save_result("spatial_only_spgau", spatial, proc.time()[["elapsed"]] - started)

  started <- proc.time()[["elapsed"]]
  combined <- capture_conditions(rma.mv(
    yi = d_Hedges, V = var_Hedges,
    random = list(~ 1 | effect_id, ~ 1 | study_id, ~ site_id | const),
    struct = "SPGAU", dist = list(site_id = p$distance_km),
    data = p$dat, method = "REML", test = "t", sparse = TRUE
  ))
  save_result("combined_spgau", combined, proc.time()[["elapsed"]] - started)
  message("Primary global Gaussian fits saved.")
}

if (stage == "profile") {
  profile_and_save <- function(fit, name, parameter) {
    captured <- if (parameter == "tau2") {
      capture_conditions(profile(fit, tau2 = 1, steps = 17, progbar = FALSE, plot = FALSE))
    } else {
      capture_conditions(profile(fit, rho = 1, steps = 17, progbar = FALSE, plot = FALSE))
    }
    if (inherits(captured$value, "error")) {
      write.csv(data.frame(parameter = parameter, error = conditionMessage(captured$value)),
                file.path(out_dir, paste0(name, "_profile_", parameter, ".csv")), row.names = FALSE)
      return(invisible(NULL))
    }
    values <- captured$value[[parameter]]
    result <- data.frame(value = as.numeric(values), logLik_REML = as.numeric(captured$value$ll))
    result$delta_logLik <- max(result$logLik_REML, na.rm = TRUE) - result$logLik_REML
    write.csv(result, file.path(out_dir, paste0(name, "_profile_", parameter, ".csv")), row.names = FALSE)
    saveRDS(captured$value, file.path(out_dir, paste0(name, "_profile_", parameter, ".rds")))
    invisible(result)
  }
  names_to_profile <- switch(profile_models,
    spatial_only = "spatial_only_spgau",
    combined = "combined_spgau",
    both = c("spatial_only_spgau", "combined_spgau")
  )
  for (name in names_to_profile) {
    fit_path <- file.path(out_dir, paste0(name, ".rds"))
    if (!file.exists(fit_path)) stop("Missing primary fit: ", fit_path)
    fit <- readRDS(fit_path)
    profile_and_save(fit, name, "tau2")
    profile_and_save(fit, name, "rho")
  }
  message("Gaussian profiles saved.")
}
