#!/usr/bin/env Rscript

# Re-estimate only the frequentist tutorial fits needed to obtain explicit
# uncertainty intervals for the package-specific result figures. This script
# is a precomputation step: it must not be called by a Quarto render and never
# overwrites the tutorial's distributed precomputed RDS objects.

suppressPackageStartupMessages({
  library(ape)
  library(glmmTMB)
  library(metadat)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[[1L]] else "all"
workers <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1L
if (!stage %in% c("all", "moura", "lim", "spain", "saved-glmmtmb") ||
    is.na(workers) || workers < 1L) {
  stop("Usage: Rscript precompute_visualisation_intervals.R {all|moura|lim|spain|saved-glmmtmb} [workers]", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
out_dir <- file.path(root, "revision_checks", "visualisation_interval_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# glmmTMB's TMB backend can use OpenMP. The number is capped by the caller;
# the metafor likelihood profiles remain sequential so that each interval is
# written immediately and can be resumed independently after an interruption.
Sys.setenv(OMP_NUM_THREADS = workers)
openmp_result <- try(glmmTMB::openmp(threads = workers), silent = TRUE)

assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(value, "error")) stop(conditionMessage(value), call. = FALSE)
  list(value = value, warnings = unique(warnings))
}

write_metadata <- function(name, lines) {
  writeLines(lines, file.path(out_dir, paste0(name, "_metadata.txt")))
}

fixed_rows <- function(fit, analysis, package) {
  if (inherits(fit, "rma.mv")) {
    ci_lb <- if (!is.null(fit$ci.lb)) as.numeric(fit$ci.lb) else as.numeric(fit$b) - qnorm(0.975) * fit$se
    ci_ub <- if (!is.null(fit$ci.ub)) as.numeric(fit$ci.ub) else as.numeric(fit$b) + qnorm(0.975) * fit$se
    data.frame(analysis = analysis, package = package, parameter = rownames(fit$b),
               estimate = as.numeric(fit$b), ci_lb = ci_lb, ci_ub = ci_ub,
               interval_method = "95% confidence interval reported by metafor", stringsAsFactors = FALSE)
  } else {
    beta <- glmmTMB::fixef(fit)$cond
    vc <- stats::vcov(fit)$cond
    se <- sqrt(diag(vc))[names(beta)]
    data.frame(analysis = analysis, package = package, parameter = names(beta),
               estimate = as.numeric(beta), ci_lb = as.numeric(beta - qnorm(0.975) * se),
               ci_ub = as.numeric(beta + qnorm(0.975) * se),
               interval_method = "95% Wald confidence interval", stringsAsFactors = FALSE)
  }
}

# glmmTMB stores Gaussian dispersion and the free random-effect scales on the
# log-SD scale. For a propto() term, its internally scaled covariance matrix
# means that exp(2 * theta) is not necessarily the reader-facing marginal
# variance. Use VarCorr() for that variance and retain the same fixed scale
# factor when transforming the Wald limits. This uses only the saved optimizer
# covariance; it does not optimise or refit the model again.
tmb_variance_rows <- function(fit, analysis, labels, groups) {
  par <- fit$fit$par
  cov_fixed <- fit$sdr$cov.fixed
  assert(identical(names(par), colnames(cov_fixed)),
         "glmmTMB parameter and covariance order disagree.")
  betadisp <- which(names(par) == "betadisp")
  theta <- which(names(par) == "theta")
  assert(length(betadisp) == 1L && length(theta) == length(groups) &&
           length(labels) == length(groups) + 1L,
         "Unexpected glmmTMB dispersion/theta component structure.")
  random <- glmmTMB::VarCorr(fit)$cond
  assert(all(groups %in% names(random)), "Expected glmmTMB random-effect group is absent.")
  estimate <- c(stats::sigma(fit)^2,
                vapply(groups, function(group) as.numeric(random[[group]][1L, 1L]), numeric(1)))
  indices <- c(betadisp, theta)
  se_log_sd <- sqrt(diag(cov_fixed)[indices])
  raw_variance <- exp(2 * par[indices])
  scale <- estimate / raw_variance
  lower_variance <- scale * exp(2 * (par[indices] - qnorm(0.975) * se_log_sd))
  upper_variance <- scale * exp(2 * (par[indices] + qnorm(0.975) * se_log_sd))
  out <- data.frame(
    analysis = analysis, package = "glmmTMB", parameter = labels,
    estimate = estimate, ci_lb = lower_variance, ci_ub = upper_variance,
    interval_method = "95% Wald confidence interval on log SD scale, squared to variance",
    stringsAsFactors = FALSE
  )
  assert(all(is.finite(as.matrix(out[c("estimate", "ci_lb", "ci_ub")]))) &&
           all(out$ci_lb >= 0) && all(out$estimate >= out$ci_lb) &&
           all(out$estimate <= out$ci_ub),
         "Invalid glmmTMB variance interval.")
  out
}

profile_rows <- function(fit, analysis, labels, type) {
  rows <- vector("list", length(labels))
  for (i in seq_along(labels)) {
    started <- proc.time()[["elapsed"]]
    result <- capture(do.call(
      metafor::confint.rma.mv,
      c(list(object = fit, level = 0.95, time = TRUE), stats::setNames(list(i), type))
    ))
    elapsed <- proc.time()[["elapsed"]] - started
    random <- as.data.frame(result$value$random[1L, , drop = FALSE])
    rows[[i]] <- data.frame(
      analysis = analysis, package = "metafor", parameter = labels[[i]],
      estimate = random$estimate, ci_lb = random$ci.lb, ci_ub = random$ci.ub,
      interval_method = "95% profile-likelihood confidence interval",
      elapsed_seconds = elapsed,
      warnings = paste(result$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
    write.csv(rows[[i]], file.path(out_dir, sprintf("%s_%s_%d_profile_ci.csv", analysis, type, i)), row.names = FALSE)
    saveRDS(result$value, file.path(out_dir, sprintf("%s_%s_%d_profile_ci.rds", analysis, type, i)))
    message(sprintf("%s %s profile %d/%d completed in %.1f seconds", analysis, type, i, length(labels), elapsed))
  }
  do.call(rbind, rows)
}

prepare_moura <- function() {
  dat <- metadat::dat.moura2021$dat
  dat$species.id.phy <- dat$species.id
  dat$effect.size.id <- factor(seq_len(nrow(dat)))
  dat <- metafor::escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  tree <- ape::compute.brlen(metadat::dat.moura2021$tree)
  A <- ape::vcv(tree, corr = TRUE)
  tips <- rownames(A)
  dat$study.id <- factor(dat$study.id)
  dat$species.id <- factor(dat$species.id)
  dat$species.id.phy <- factor(as.character(dat$species.id.phy), levels = tips)
  dat$g <- factor("all")
  V <- diag(dat$vi)
  rownames(V) <- colnames(V) <- levels(dat$effect.size.id)
  assert(!anyNA(dat$species.id.phy) && identical(levels(dat$species.id.phy), tips),
         "Moura phylogenetic factor is not aligned to A.")
  assert(identical(rownames(V), levels(dat$effect.size.id)), "Moura sampling VCV is misaligned.")
  list(dat = dat, A = A, V = V)
}

fit_moura <- function() {
  x <- prepare_moura()
  mf_capture <- capture(metafor::rma.mv(
    yi, vi,
    random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
    R = list(species.id.phy = x$A), data = x$dat, sparse = TRUE, method = "REML", test = "t"
  ))
  gt_capture <- capture(glmmTMB::glmmTMB(
    yi ~ 1 + equalto(0 + effect.size.id | g, x$V) + (1 | study.id) +
      (1 | species.id) + propto(0 + species.id.phy | g, x$A),
    data = x$dat, REML = TRUE
  ))
  mf <- mf_capture$value
  gt <- gt_capture$value
  saveRDS(mf, file.path(out_dir, "moura_metafor_bm.rds"))
  saveRDS(gt, file.path(out_dir, "moura_glmmTMB_bm.rds"))
  metafor_component_labels <- c("Study variance", "Effect-size variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic")
  glmmTMB_component_labels <- c("Effect-size variance", "Study variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic")
  prof <- profile_rows(mf, "moura", metafor_component_labels, "sigma2")
  out <- rbind(
    fixed_rows(mf, "moura", "metafor"),
    fixed_rows(gt, "moura", "glmmTMB"),
    prof[, c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")],
    tmb_variance_rows(gt, "moura", glmmTMB_component_labels,
                      c("study.id", "species.id", "g.1"))
  )
  write.csv(out, file.path(out_dir, "moura_frequentist_intervals.csv"), row.names = FALSE)
  write_metadata("moura", c(
    "Re-estimated on Totoro for package-specific tutorial figure intervals.",
    paste("metafor warnings:", paste(mf_capture$warnings, collapse = " | ")),
    paste("glmmTMB warnings:", paste(gt_capture$warnings, collapse = " | ")),
    paste("glmmTMB pdHess:", isTRUE(gt$sdr$pdHess)),
    paste("glmmTMB optimizer convergence:", gt$fit$convergence),
    paste("OpenMP request:", workers)
  ))
  assert(abs(mf$b[1] - glmmTMB::fixef(gt)$cond[[1L]]) < 0.005,
         "Moura pooled estimates do not agree across packages.")
  message("MOURA_FREQUENTIST_INTERVAL_PRECOMPUTE_PASSED")
}

prepare_lim <- function() {
  dat <- metadat::dat.lim2014$o_o_unadj
  tree <- ape::compute.brlen(metadat::dat.lim2014$o_o_unadj_tree)
  dat <- metafor::escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  A <- ape::vcv(tree, corr = TRUE)
  tips <- rownames(A)
  dat$species <- factor(as.character(dat$species))
  dat$phy <- factor(as.character(dat$species), levels = tips)
  dat$id <- factor(seq_len(nrow(dat)))
  dat$environment <- factor(dat$environment)
  dat$g <- factor("all")
  V <- diag(dat$vi)
  rownames(V) <- colnames(V) <- levels(dat$id)
  assert(!anyNA(dat$phy) && identical(levels(dat$phy), tips), "Lim phylogenetic factor is not aligned to A.")
  list(dat = dat, A = A, V = V)
}

fit_lim <- function() {
  x <- prepare_lim()
  mf_capture <- capture(metafor::rma.mv(
    yi, vi, mods = ~ environment,
    random = list(~ 1 | id, ~ 1 | phy, ~ 1 | species), R = list(phy = x$A),
    data = x$dat, sparse = TRUE, method = "REML", test = "z"
  ))
  gt_capture <- capture(glmmTMB::glmmTMB(
    yi ~ 1 + environment + equalto(0 + id | g, x$V) + (1 | species) +
      propto(0 + phy | g, x$A), data = x$dat, REML = TRUE
  ))
  mf <- mf_capture$value
  gt <- gt_capture$value
  saveRDS(mf, file.path(out_dir, "lim_metafor_bm_meta_regression.rds"))
  saveRDS(gt, file.path(out_dir, "lim_glmmTMB_bm_meta_regression.rds"))
  metafor_component_labels <- c("Effect-size variance", "Species variance, phylogenetic", "Species variance, non-phylogenetic")
  glmmTMB_component_labels <- c("Effect-size variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic")
  prof <- profile_rows(mf, "lim", metafor_component_labels, "sigma2")
  out <- rbind(
    fixed_rows(mf, "lim", "metafor"),
    fixed_rows(gt, "lim", "glmmTMB"),
    prof[, c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")],
    tmb_variance_rows(gt, "lim", glmmTMB_component_labels, c("species", "g.1"))
  )
  write.csv(out, file.path(out_dir, "lim_frequentist_intervals.csv"), row.names = FALSE)
  write_metadata("lim", c(
    "Re-estimated on Totoro for package-specific tutorial figure intervals.",
    paste("metafor warnings:", paste(mf_capture$warnings, collapse = " | ")),
    paste("glmmTMB warnings:", paste(gt_capture$warnings, collapse = " | ")),
    paste("glmmTMB pdHess:", isTRUE(gt$sdr$pdHess)),
    paste("glmmTMB optimizer convergence:", gt$fit$convergence),
    paste("OpenMP request:", workers)
  ))
  assert(max(abs(as.numeric(mf$b) - as.numeric(glmmTMB::fixef(gt)$cond))) < 0.005,
         "Lim fixed-effect estimates do not agree across packages.")
  message("LIM_FREQUENTIST_INTERVAL_PRECOMPUTE_PASSED")
}

fit_spain <- function() {
  input <- file.path(root, "revision_checks", "regional_cross_package_audit_outputs", "spain_prepared.rds")
  x <- readRDS(input)
  dat <- x$dat
  distance_km <- x$D_proj
  assert(nrow(dat) == 186L && nlevels(dat$study_id) == 30L && nlevels(dat$site_id) == 32L,
         "Unexpected prepared Spain data.")
  assert(is.matrix(distance_km) &&
           identical(rownames(distance_km), levels(dat$site_id)) &&
           identical(colnames(distance_km), levels(dat$site_id)),
         "Spain projected distance matrix is missing or misaligned.")
  fit_capture <- capture(metafor::rma.mv(
    yi = d_Hedges, V = var_Hedges,
    random = list(~ 1 | effect_id, ~ site_id | const), struct = "SPEXP",
    dist = list(site_id = distance_km), data = dat, method = "REML", test = "t", sparse = TRUE,
    control = list(REMLf = FALSE)
  ))
  fit <- fit_capture$value
  saveRDS(fit, file.path(out_dir, "spain_metafor_spatial_only.rds"))
  sigma <- profile_rows(fit, "spain", "IID effect-size variance", "sigma2")
  tau <- profile_rows(fit, "spain", "Spatial variance", "tau2")
  rho <- profile_rows(fit, "spain", "Exponential range (km)", "rho")
  out <- rbind(
    fixed_rows(fit, "spain", "metafor"),
    sigma[, c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")],
    tau[, c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")],
    rho[, c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")]
  )
  write.csv(out, file.path(out_dir, "spain_metafor_intervals.csv"), row.names = FALSE)
  write_metadata("spain", c(
    "Re-estimated on Totoro to add the missing iid variance profile interval.",
    paste("metafor warnings:", paste(fit_capture$warnings, collapse = " | ")),
    paste("fitted logLik REML:", as.numeric(fit$fit.stats["ll", "REML"]))
  ))
  message("SPAIN_METAFOR_INTERVAL_PRECOMPUTE_PASSED")
}

# Correctly derive reader-facing glmmTMB component variances from saved fits
# when profile/fitting artifacts were produced by an earlier interrupted run.
# This path is intentionally non-executing with respect to model fitting.
refresh_saved_glmmtmb_intervals <- function() {
  specs <- list(
    moura = list(
      fit = "moura_glmmTMB_bm.rds",
      csv = "moura_frequentist_intervals.csv",
      labels = c("Effect-size variance", "Study variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic"),
      groups = c("study.id", "species.id", "g.1")
    ),
    lim = list(
      fit = "lim_glmmTMB_bm_meta_regression.rds",
      csv = "lim_frequentist_intervals.csv",
      labels = c("Effect-size variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic"),
      groups = c("species", "g.1")
    )
  )
  for (name in names(specs)) {
    spec <- specs[[name]]
    fit_path <- file.path(out_dir, spec$fit)
    csv_path <- file.path(out_dir, spec$csv)
    assert(file.exists(fit_path) && file.exists(csv_path),
           paste("Missing saved glmmTMB interval artifact for", name))
    fit <- readRDS(fit_path)
    existing <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    replacement <- tmb_variance_rows(fit, name, spec$labels, spec$groups)
    existing <- existing[, names(replacement), drop = FALSE]
    retained <- existing[!(existing$package == "glmmTMB" &
                           existing$parameter %in% spec$labels), , drop = FALSE]
    updated <- rbind(retained, replacement)
    assert(sum(updated$package == "glmmTMB" &
                 updated$parameter %in% spec$labels) == length(spec$labels),
           paste("Could not replace all glmmTMB component intervals for", name))
    write.csv(updated, csv_path, row.names = FALSE)
  }
  message("SAVED_GLMMTMB_INTERVAL_REFRESH_PASSED")
}

if (stage %in% c("all", "moura")) fit_moura()
if (stage %in% c("all", "lim")) fit_lim()
if (stage %in% c("all", "spain")) fit_spain()
if (stage == "saved-glmmtmb") refresh_saved_glmmtmb_intervals()

writeLines(c(
  paste("stage:", stage),
  paste("requested OpenMP workers:", workers),
  paste("R:", R.version.string),
  paste("metafor:", as.character(packageVersion("metafor"))),
  paste("glmmTMB:", as.character(packageVersion("glmmTMB"))),
  paste("ape:", as.character(packageVersion("ape"))),
  paste("metadat:", as.character(packageVersion("metadat")))
), file.path(out_dir, "run_environment.txt"))
message("VISUALISATION_INTERVAL_PRECOMPUTE_PASSED")
