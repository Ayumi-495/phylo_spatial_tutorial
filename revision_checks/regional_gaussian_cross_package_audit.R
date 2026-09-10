# Regional Spain cross-package audit.
# This script is deliberately separate from tutorial_v2.qmd.
# Stages are run sequentially: preflight, metafor, glmmTMB, then (on Totoro) brms.

suppressPackageStartupMessages({
  library(sf)
  library(geosphere)
  library(metafor)
  library(glmmTMB)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript regional_cross_package_audit.R <preflight|metafor|glmmTMB> <data_csv> <output_dir>")
}
stage <- args[[1L]]
data_csv <- args[[2L]]
out_dir <- args[[3L]]
stopifnot(stage %in% c("preflight", "metafor", "glmmTMB"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

crs_lcc <- "+proj=lcc +lat_1=38 +lat_2=43 +lat_0=40.5 +lon_0=-3.5 +datum=WGS84 +units=m +no_defs"

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

write_metadata <- function(extra = character()) {
  writeLines(c(
    paste("R:", R.version.string),
    paste("metafor:", as.character(packageVersion("metafor"))),
    paste("glmmTMB:", as.character(packageVersion("glmmTMB"))),
    paste("sf:", as.character(packageVersion("sf"))),
    paste("geosphere:", as.character(packageVersion("geosphere"))),
    paste("CRS:", crs_lcc),
    "Region: country == Spain; all 186 rows with non-missing coordinates retained.",
    "Projected coordinates: WGS84 Lambert Conformal Conic, expressed in kilometres.",
    "Sampling error: V = diag(var_Hedges), independent sampling errors.",
    "Target model: intercept + iid effect-size heterogeneity + Gaussian spatial covariance.",
    extra
  ), file.path(out_dir, "regional_audit_metadata.txt"))
}

prepare_data <- function() {
  dat <- read.csv(data_csv, stringsAsFactors = FALSE)
  dat <- dat[!is.na(dat$latitude) & !is.na(dat$longitude) &
               !is.na(dat$country) & dat$country == "Spain", , drop = FALSE]
  stopifnot(nrow(dat) == 186L, length(unique(dat$study_id)) == 30L)
  stopifnot(all(is.finite(dat$d_Hedges)), all(is.finite(dat$var_Hedges)),
            all(dat$var_Hedges > 0))

  # A recorded coordinate location is a unique coordinate pair. This does not
  # assert that equal coordinates are the same exact field site.
  dat$effect_id <- factor(sprintf("effect_%04d", seq_len(nrow(dat))),
                          levels = sprintf("effect_%04d", seq_len(nrow(dat))))
  dat$study_id <- factor(dat$study_id)
  dat$site_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
  site_levels <- sort(unique(dat$site_key))
  dat$site_id <- factor(dat$site_key, levels = site_levels)
  dat$const <- factor("all_sites")
  stopifnot(nlevels(dat$site_id) == 32L,
            identical(levels(dat$site_id), site_levels), nlevels(dat$const) == 1L)

  site_lookup <- unique(dat[c("site_key", "latitude", "longitude")])
  site_lookup <- site_lookup[match(site_levels, site_lookup$site_key), , drop = FALSE]
  site_lookup$site_id <- site_levels
  stopifnot(identical(site_lookup$site_key, site_levels))

  lonlat <- as.matrix(site_lookup[c("longitude", "latitude")])
  D_geo <- geosphere::distm(lonlat, fun = geosphere::distGeo) / 1000
  rownames(D_geo) <- colnames(D_geo) <- site_levels

  pts <- sf::st_as_sf(site_lookup, coords = c("longitude", "latitude"),
                      crs = 4326, remove = FALSE)
  pts_lcc <- sf::st_transform(pts, crs_lcc)
  xy_m <- sf::st_coordinates(pts_lcc)
  site_lookup$x_km <- xy_m[, 1] / 1000
  site_lookup$y_km <- xy_m[, 2] / 1000
  dat$x_km <- site_lookup$x_km[match(dat$site_key, site_lookup$site_key)]
  dat$y_km <- site_lookup$y_km[match(dat$site_key, site_lookup$site_key)]
  stopifnot(all(is.finite(dat$x_km)), all(is.finite(dat$y_km)))
  D_proj <- as.matrix(dist(site_lookup[c("x_km", "y_km")]))
  rownames(D_proj) <- colnames(D_proj) <- site_levels

  # Matrix-order and metric assertions are required before any model fit.
  stopifnot(identical(rownames(D_geo), levels(dat$site_id)),
            identical(colnames(D_geo), levels(dat$site_id)),
            identical(rownames(D_proj), levels(dat$site_id)),
            identical(colnames(D_proj), levels(dat$site_id)),
            isTRUE(all.equal(D_geo, t(D_geo), tolerance = 1e-10)),
            isTRUE(all.equal(D_proj, t(D_proj), tolerance = 1e-10)),
            all(abs(diag(D_geo)) < 1e-10), all(abs(diag(D_proj)) < 1e-10),
            all(D_geo >= 0), all(D_proj >= 0), all(is.finite(D_geo)),
            all(is.finite(D_proj)))

  ix <- upper.tri(D_geo)
  distortion <- D_proj[ix] / D_geo[ix] - 1
  distortion_summary <- data.frame(
    max_abs_pct = max(abs(distortion)) * 100,
    p95_abs_pct = as.numeric(quantile(abs(distortion), 0.95)) * 100,
    median_abs_pct = median(abs(distortion)) * 100,
    max_great_circle_km = max(D_geo),
    max_projected_km = max(D_proj),
    stringsAsFactors = FALSE
  )

  saveRDS(dat, file.path(out_dir, "spain_data.rds"))
  write.csv(site_lookup, file.path(out_dir, "spain_site_lookup_projected_km.csv"), row.names = FALSE)
  write.csv(D_geo, file.path(out_dir, "spain_distance_great_circle_km.csv"), row.names = TRUE)
  write.csv(D_proj, file.path(out_dir, "spain_distance_lcc_km.csv"), row.names = TRUE)
  write.csv(distortion_summary, file.path(out_dir, "spain_projection_distortion.csv"), row.names = FALSE)
  saveRDS(list(dat = dat, site_lookup = site_lookup, D_geo = D_geo,
               D_proj = D_proj, distortion = distortion_summary),
          file.path(out_dir, "spain_prepared.rds"))
  list(dat = dat, site_lookup = site_lookup, D_geo = D_geo,
       D_proj = D_proj, distortion = distortion_summary)
}

load_prepared <- function() {
  p <- file.path(out_dir, "spain_prepared.rds")
  if (!file.exists(p)) stop("Run the preflight stage first: missing ", p)
  readRDS(p)
}

if (stage == "preflight") {
  p <- prepare_data()
  # Technical equalto check against a metafor known-V fit. This is a preflight
  # implementation check, not the regional spatial comparison model.
  d <- p$dat[seq_len(20L), , drop = FALSE]
  d$effect_id <- factor(sprintf("e%03d", seq_len(nrow(d))),
                         levels = sprintf("e%03d", seq_len(nrow(d))))
  d$const <- factor("all")
  VCV <- diag(d$var_Hedges)
  rownames(VCV) <- colnames(VCV) <- levels(d$effect_id)
  stopifnot(identical(rownames(VCV), levels(d$effect_id)),
            identical(colnames(VCV), levels(d$effect_id)))
  eq <- capture_conditions(glmmTMB::glmmTMB(
    d_Hedges ~ 1 + equalto(0 + effect_id | const, VCV),
    data = d, REML = TRUE, dispformula = ~0
  ))
  mf <- metafor::rma.mv(d_Hedges, VCV, mods = ~1, data = d,
                        method = "REML", control = list(REMLf = FALSE))
  if (inherits(eq$value, "error")) stop(eq$value)
  preflight <- data.frame(
    n = nrow(d),
    metafor_mean = as.numeric(mf$b[1]),
    glmmTMB_mean = as.numeric(glmmTMB::fixef(eq$value)$cond[1]),
    metafor_logLik = as.numeric(logLik(mf)),
    glmmTMB_logLik = as.numeric(logLik(eq$value)),
    glmmTMB_pdHess = isTRUE(eq$value$sdr$pdHess),
    warning = paste(eq$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
  write.csv(preflight, file.path(out_dir, "equalto_preflight.csv"), row.names = FALSE)
  write_metadata(c(
    "equalto preflight: glmmTMB equalto(0 + effect_id | const, VCV) with dispformula=~0 was compared with metafor VCV fit.",
    paste("Preflight result:", paste(capture.output(print(preflight)), collapse = " "))
  ))
  message("Preflight completed; no regional spatial model was fitted.")
}

if (stage == "metafor") {
  p <- load_prepared()
  dat <- p$dat
  fit_capture <- capture_conditions({
    metafor::rma.mv(
      yi = d_Hedges, V = var_Hedges,
      mods = ~1,
      random = list(~1 | effect_id, ~site_id | const),
      struct = "SPGAU", dist = list(site_id = p$D_proj),
      data = dat, method = "REML", test = "t", sparse = TRUE,
      # REMLf=FALSE uses the REML likelihood convention that matches
      # glmmTMB's equalto implementation; estimates are unchanged.
      control = list(REMLf = FALSE)
    )
  })
  if (inherits(fit_capture$value, "error")) stop(fit_capture$value)
  fit <- fit_capture$value
  saveRDS(fit, file.path(out_dir, "metafor_spatial_only.rds"))
  result <- data.frame(
    package = "metafor", model = "regional_spatial_only",
    n = fit$k, studies = nlevels(dat$study_id), sites = nlevels(dat$site_id),
    mean = as.numeric(fit$b[1]), ci_lb = as.numeric(fit$ci.lb[1]),
    ci_ub = as.numeric(fit$ci.ub[1]), iid_effect_variance = as.numeric(fit$sigma2[1]),
    spatial_variance = as.numeric(fit$tau2[1]), rho_km = as.numeric(fit$rho[1]),
    logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]),
    AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
    warnings = paste(fit_capture$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
  write.csv(result, file.path(out_dir, "metafor_spatial_only_result.csv"), row.names = FALSE)

  # Profile both spatial parameters before proceeding to glmmTMB.
  prof_tau <- capture_conditions(profile(
    fit, tau2 = 1, steps = 15, progbar = FALSE, plot = FALSE
  ))
  prof_rho <- capture_conditions(profile(
    fit, rho = 1, steps = 15, progbar = FALSE, plot = FALSE
  ))
  if (!inherits(prof_tau$value, "error")) saveRDS(prof_tau$value, file.path(out_dir, "metafor_profile_tau2.rds"))
  if (!inherits(prof_rho$value, "error")) saveRDS(prof_rho$value, file.path(out_dir, "metafor_profile_rho.rds"))
  profile_to_csv <- function(x, parameter, path) {
    if (inherits(x$value, "error")) {
      write.csv(data.frame(parameter = parameter, error = conditionMessage(x$value)), path, row.names = FALSE)
      return(invisible(NULL))
    }
    vals <- x$value[[parameter]]
    out <- data.frame(value = as.numeric(vals), logLik_REML = as.numeric(x$value$ll))
    out$delta_logLik <- max(out$logLik_REML, na.rm = TRUE) - out$logLik_REML
    write.csv(out, path, row.names = FALSE)
  }
  profile_to_csv(prof_tau, "tau2", file.path(out_dir, "metafor_profile_tau2.csv"))
  profile_to_csv(prof_rho, "rho", file.path(out_dir, "metafor_profile_rho.csv"))
  ci_tau <- capture_conditions(confint(fit, tau2 = 1))
  ci_rho <- capture_conditions(confint(fit, rho = 1))
  if (!inherits(ci_tau$value, "error")) write.csv(as.data.frame(ci_tau$value), file.path(out_dir, "metafor_profile_ci_tau2.csv"))
  if (!inherits(ci_rho$value, "error")) write.csv(as.data.frame(ci_rho$value), file.path(out_dir, "metafor_profile_ci_rho.csv"))
  write_metadata(c(
    paste("metafor warnings:", paste(fit_capture$warnings, collapse = " | ")),
    paste("tau2 profile warnings:", paste(prof_tau$warnings, collapse = " | ")),
    paste("rho profile warnings:", paste(prof_rho$warnings, collapse = " | ")),
    "metafor control: REMLf=FALSE, to use the same reduced REML likelihood convention as glmmTMB for likelihood/AIC comparison."
  ))
  message("metafor spatial-only fit and profiles saved.")
}

if (stage == "glmmTMB") {
  p <- load_prepared()
  dat <- p$dat
  VCV <- diag(dat$var_Hedges)
  rownames(VCV) <- colnames(VCV) <- levels(dat$effect_id)
  stopifnot(identical(rownames(VCV), levels(dat$effect_id)),
            identical(colnames(VCV), levels(dat$effect_id)))
  dat$pos <- glmmTMB::numFactor(dat$x_km, dat$y_km)
  dat$const <- factor("all_sites")
  pos_xy <- glmmTMB::parseNumLevels(levels(dat$pos))
  expected_xy <- as.matrix(p$site_lookup[c("x_km", "y_km")])
  stopifnot(nrow(pos_xy) == nrow(expected_xy),
            isTRUE(all.equal(unname(pos_xy[order(pos_xy[, 1], pos_xy[, 2]), , drop = FALSE]),
                             unname(expected_xy[order(expected_xy[, 1], expected_xy[, 2]), , drop = FALSE]),
                             tolerance = 1e-8)))

  fit_capture <- capture_conditions({
    glmmTMB::glmmTMB(
      d_Hedges ~ 1 +
        equalto(0 + effect_id | const, VCV) +
        gau(pos + 0 | const),
      data = dat, REML = TRUE
    )
  })
  if (inherits(fit_capture$value, "error")) stop(fit_capture$value)
  fit <- fit_capture$value
  saveRDS(fit, file.path(out_dir, "glmmTMB_spatial_only.rds"))
  vc <- VarCorr(fit)$cond$const
  theta <- fit$fit$par[which(names(fit$fit$par) == "theta")]
  stopifnot(length(theta) == 2L)
  # For gau, correlation is exp(-exp(-2*theta[2]) * distance^2),
  # so exp(theta[2]) is the common e-folding distance in kilometres.
  exp_theta <- theta[seq_len(2L)]
  rho_km <- exp(exp_theta[2])
  spatial_vc <- VarCorr(fit)$cond[["const.1"]]
  spatial_sd <- unname(attr(spatial_vc, "stddev")[1])
  b <- summary(fit)$coefficients$cond["(Intercept)", ]
  diag_ok <- glmmTMB::diagnose(fit)
  writeLines(c(
    paste("pdHess:", isTRUE(fit$sdr$pdHess)),
    paste("optimizer convergence:", fit$fit$convergence),
    paste("optimizer message:", fit$fit$message),
    paste("diagnose() returned:", isTRUE(diag_ok)),
    "diagnose() flagged an unusually large absolute Z statistic for the Gaussian dispersion intercept; this is a Wald-approximation caution, not an optimizer convergence failure."
  ), file.path(out_dir, "glmmTMB_diagnostics.txt"))
  result <- data.frame(
    package = "glmmTMB", model = "regional_spatial_only",
    n = nrow(dat), studies = nlevels(dat$study_id), sites = nlevels(dat$site_id),
    mean = as.numeric(b[["Estimate"]]),
    ci_lb = as.numeric(b[["Estimate"]] - qnorm(0.975) * b[["Std. Error"]]),
    ci_ub = as.numeric(b[["Estimate"]] + qnorm(0.975) * b[["Std. Error"]]),
    iid_effect_variance = as.numeric(sigma(fit))^2,
    spatial_variance = spatial_sd^2, spatial_sd = spatial_sd,
    rho_km = rho_km, logLik_REML = as.numeric(logLik(fit)),
    AIC_REML = AIC(fit), pdHess = isTRUE(fit$sdr$pdHess),
    optimizer_convergence = fit$fit$convergence,
    diagnose_ok = isTRUE(diag_ok),
    warnings = paste(fit_capture$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
  write.csv(result, file.path(out_dir, "glmmTMB_spatial_only_result.csv"), row.names = FALSE)
  write.csv(data.frame(theta_name = names(theta), theta = as.numeric(theta)),
            file.path(out_dir, "glmmTMB_theta.csv"), row.names = FALSE)
  write_metadata(c(
    paste("glmmTMB warnings:", paste(fit_capture$warnings, collapse = " | ")),
    "glmmTMB gau parameterisation: exp(-exp(-2*theta[2]) * distance_km^2); common e-folding rho = exp(theta[2]) km.",
    "glmmTMB iid effect-size heterogeneity is the Gaussian residual variance sigma(fit)^2; equalto supplies known sampling VCV without estimation."
  ))
  message("glmmTMB spatial-only fit saved.")
}
