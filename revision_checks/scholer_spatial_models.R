# Corrected Scholer spatial metafor audit.
# Keeps all fitting and profile work outside tutorial_v2.qmd.
# Usage:
#   Rscript scholer_spatial_models.R primary <prepared.rds> <output_dir>
#   Rscript scholer_spatial_models.R profile_point <prepared.rds> <output_dir> <spatial_only|combined> <tau2|rho> <fixed_value>

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) == 3L || length(args) == 6L)) {
  stop("Invalid arguments; see usage in file header.")
}
stage <- args[[1L]]
prepared_path <- args[[2L]]
out_dir <- args[[3L]]
if (!(stage %in% c("primary", "profile_point"))) stop("Unknown stage: ", stage)
if (stage == "primary" && length(args) != 3L) stop("primary takes exactly three arguments.")
if (stage == "profile_point" && length(args) != 6L) stop("profile_point takes six arguments.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "profiles"), recursive = TRUE, showWarnings = FALSE)

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
  if (any(grepl("converg|Hessian|optim", warnings, ignore.case = TRUE))) return("optimizer_warning")
  "completed_no_explicit_optimizer_status"
}

fit_stat <- function(fit, statistic) as.numeric(fit$fit.stats[statistic, "REML"])

upsert_csv <- function(row, path, keys) {
  previous <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  if (!is.null(previous) && nrow(previous)) {
    keep <- rep(TRUE, nrow(previous))
    for (key in keys) keep <- keep & previous[[key]] != row[[key]][1L]
    previous <- previous[keep, , drop = FALSE]
  }
  write.csv(rbind(previous, row), path, row.names = FALSE)
}

prepared <- readRDS(prepared_path)
dat <- prepared$dat
distance_km <- prepared$distance_km
stopifnot(
  nrow(dat) == 949L, nlevels(dat$effect_id) == 949L,
  nlevels(dat$study_id) == 205L, nlevels(dat$site_id) == 454L,
  all(is.finite(dat$logit_survival)), all(is.finite(dat$vi)), all(dat$vi > 0),
  nlevels(dat$const) == 1L
)

assert_distance_order <- function() {
  stopifnot(
    identical(rownames(distance_km), levels(dat$site_id)),
    identical(colnames(distance_km), levels(dat$site_id)),
    isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
    all(abs(diag(distance_km)) < 1e-10), all(distance_km >= 0), all(is.finite(distance_km))
  )
}

write_metadata <- function() {
  writeLines(c(
    paste("R:", R.version.string),
    paste("metafor:", as.character(packageVersion("metafor"))),
    "Dataset: Scholer et al. full data; 949 effect sizes, 205 references, 454 recorded coordinate locations.",
    "Distance: WGS84 ellipsoidal great-circle kilometres from geosphere::distGeo(), inherited from the audited prepared data.",
    "Sampling variance: vi = se^2; V is diagonal under the retained independent-sampling-error assumption.",
    "Fixed effects: intercept only. Common estimation/test settings: method=REML, test=t, sparse=TRUE.",
    "Spatial outer group: const has one level, allowing covariance across studies.",
    "No iid location intercept is fitted."
  ), file.path(out_dir, "scholer_model_metadata.txt"))
}

fit_model <- function(model, tau2 = NULL, rho = NULL) {
  # Required immediately before every fit, including the unstructured model.
  assert_distance_order()
  args <- list(
    yi = dat$logit_survival, V = dat$vi, mods = ~ 1,
    data = dat, method = "REML", test = "t", sparse = TRUE
  )
  if (model == "unstructured_only") {
    args$random <- list(~ 1 | effect_id, ~ 1 | study_id)
  } else if (model == "spatial_only") {
    args$random <- list(~ 1 | effect_id, ~ site_id | const)
    args$struct <- "SPEXP"
    args$dist <- list(site_id = distance_km)
  } else if (model == "combined") {
    args$random <- list(~ 1 | effect_id, ~ 1 | study_id, ~ site_id | const)
    args$struct <- "SPEXP"
    args$dist <- list(site_id = distance_km)
  } else {
    stop("Unknown model: ", model)
  }
  if (!is.null(tau2)) args$tau2 <- tau2
  if (!is.null(rho)) args$rho <- rho
  do.call(metafor::rma.mv, args)
}

primary_row <- function(model, result, elapsed_seconds) {
  if (inherits(result$value, "error")) {
    return(data.frame(
      model = model, n = 949L, references = 205L, sites = 454L,
      pooled_mean = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
      iid_effect_variance = NA_real_, study_variance = NA_real_,
      spatial_variance = NA_real_, rho_km = NA_real_,
      logLik_REML = NA_real_, AIC_REML = NA_real_,
      convergence_status = "error", elapsed_seconds = elapsed_seconds,
      warnings = paste(c(conditionMessage(result$value), result$warnings), collapse = " | "),
      stringsAsFactors = FALSE
    ))
  }
  fit <- result$value
  sigma2 <- fit$sigma2
  data.frame(
    model = model, n = fit$k, references = nlevels(dat$study_id), sites = nlevels(dat$site_id),
    pooled_mean = as.numeric(fit$b[1]), ci_lb = as.numeric(fit$ci.lb[1]), ci_ub = as.numeric(fit$ci.ub[1]),
    iid_effect_variance = if (length(sigma2) >= 1L) as.numeric(sigma2[1]) else NA_real_,
    study_variance = if (length(sigma2) >= 2L) as.numeric(sigma2[2]) else NA_real_,
    spatial_variance = if (length(fit$tau2)) as.numeric(fit$tau2[1]) else NA_real_,
    rho_km = if (length(fit$rho)) as.numeric(fit$rho[1]) else NA_real_,
    logLik_REML = fit_stat(fit, "ll"), AIC_REML = fit_stat(fit, "AIC"),
    convergence_status = fit_status(fit, result$warnings), elapsed_seconds = elapsed_seconds,
    warnings = paste(result$warnings, collapse = " | "), stringsAsFactors = FALSE
  )
}

if (stage == "primary") {
  write_metadata()
  for (model in c("unstructured_only", "spatial_only", "combined")) {
    started <- proc.time()[["elapsed"]]
    result <- capture_conditions(fit_model(model))
    elapsed <- proc.time()[["elapsed"]] - started
    row <- primary_row(model, result, elapsed)
    if (!inherits(result$value, "error")) {
      saveRDS(result$value, file.path(out_dir, paste0(model, ".rds")))
    }
    upsert_csv(row, file.path(out_dir, "scholer_primary_model_results.csv"), "model")
    message(model, ": ", row$convergence_status, "; elapsed ", round(elapsed, 2), " sec")
  }
}

if (stage == "profile_point") {
  model <- args[[4L]]
  parameter <- args[[5L]]
  fixed_value <- as.numeric(args[[6L]])
  if (!(model %in% c("spatial_only", "combined"))) stop("Profile model must be spatial_only or combined.")
  if (!(parameter %in% c("tau2", "rho"))) stop("Profile parameter must be tau2 or rho.")
  if (!is.finite(fixed_value) || fixed_value < 0) stop("Profile value must be finite and non-negative.")
  started <- proc.time()[["elapsed"]]
  result <- if (parameter == "tau2") capture_conditions(fit_model(model, tau2 = fixed_value)) else capture_conditions(fit_model(model, rho = fixed_value))
  elapsed <- proc.time()[["elapsed"]] - started
  tag <- gsub("[^0-9A-Za-z._-]", "_", format(fixed_value, scientific = FALSE, trim = TRUE))
  output_path <- file.path(out_dir, "profiles", paste0(model, "_", parameter, "_", tag, ".rds"))
  if (inherits(result$value, "error")) {
    row <- data.frame(model = model, parameter = parameter, fixed_value = fixed_value,
                      estimated_tau2 = NA_real_, estimated_rho_km = NA_real_,
                      logLik_REML = NA_real_, AIC_REML = NA_real_, elapsed_seconds = elapsed,
                      convergence_status = "error",
                      warnings = paste(c(conditionMessage(result$value), result$warnings), collapse = " | "),
                      stringsAsFactors = FALSE)
  } else {
    fit <- result$value
    saveRDS(fit, output_path)
    row <- data.frame(model = model, parameter = parameter, fixed_value = fixed_value,
                      estimated_tau2 = if (length(fit$tau2)) as.numeric(fit$tau2[1]) else NA_real_,
                      estimated_rho_km = if (length(fit$rho)) as.numeric(fit$rho[1]) else NA_real_,
                      logLik_REML = fit_stat(fit, "ll"), AIC_REML = fit_stat(fit, "AIC"),
                      elapsed_seconds = elapsed, convergence_status = fit_status(fit, result$warnings),
                      warnings = paste(result$warnings, collapse = " | "), stringsAsFactors = FALSE)
  }
  # Parallel profile workers always write their own CSV/RDS immediately. Set
  # PROFILE_APPEND=false to avoid concurrent writes to the compiled CSV; a
  # later serial compilation step can combine the per-point CSV files.
  if (Sys.getenv("PROFILE_APPEND", "true") != "false") {
    upsert_csv(row, file.path(out_dir, "scholer_profile_points.csv"), c("model", "parameter", "fixed_value"))
  }
  write.csv(row, file.path(out_dir, "profiles", paste0(model, "_", parameter, "_", tag, ".csv")), row.names = FALSE)
  message(model, " ", parameter, "=", fixed_value, ": ", row$convergence_status,
          "; elapsed ", round(elapsed, 2), " sec")
}
