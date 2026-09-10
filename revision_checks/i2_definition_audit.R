# Reproducible I2 audit for models whose variance structures are already settled.
#
# This script does not refit any model. It reads saved fit objects and computes
# a generalized representative sampling variance from the known sampling V and
# fixed-effect design matrix X:
#
#   W = V^{-1}
#   P = W - W X (X' W X)^{-1} X' W
#   v_tilde = (k - p) / trace(P)
#
# I2 values are reported as percentages. Correlated random-effect variances are
# included only when their correlation matrices have unit diagonals, so the
# fitted scalar variance is also the component's per-observation marginal
# variance. Correlation between distinct observations is not summarized by I2.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(brms)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
checks_dir <- file.path(project_root, "revision_checks")
out_dir <- file.path(checks_dir, "i2_definition_audit_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

generalized_sampling_variance <- function(V, X, tolerance = 1e-12) {
  X <- as.matrix(X)
  k <- nrow(X)
  p <- ncol(X)
  stopifnot(k > p, qr(X)$rank == p)

  if (is.null(dim(V))) {
    stopifnot(length(V) == k, all(is.finite(V)), all(V > 0))
    v <- as.numeric(V)
    diagonal_V <- TRUE
  } else {
    V <- as.matrix(V)
    stopifnot(identical(dim(V), c(k, k)),
              isTRUE(all.equal(V, t(V), tolerance = tolerance)),
              all(is.finite(V)))
    off_diagonal <- V
    diag(off_diagonal) <- 0
    diagonal_V <- max(abs(off_diagonal)) <= tolerance
    if (diagonal_V) {
      v <- diag(V)
      stopifnot(all(v > 0))
    }
  }

  # Compute trace(P) without allocating P. For symmetric W,
  # trace[W X (X' W X)^-1 X' W]
  # = trace[(X' W X)^-1 X' W^2 X].
  if (diagonal_V) {
    w <- 1 / v
    WX <- X * w
    XtWX <- crossprod(X, WX)
    XtW2X <- crossprod(WX)
    trace_P <- sum(w) - sum(diag(solve(XtWX, XtW2X)))
  } else {
    chol_V <- chol(V)
    W <- chol2inv(chol_V)
    WX <- W %*% X
    XtWX <- crossprod(X, WX)
    XtW2X <- crossprod(WX)
    trace_P <- sum(diag(W)) - sum(diag(solve(XtWX, XtW2X)))
  }

  stopifnot(is.finite(trace_P), trace_P > 0)
  list(
    k = k,
    p = p,
    trace_P = trace_P,
    v_tilde = (k - p) / trace_P,
    diagonal_V = diagonal_V
  )
}

sampling_variance_summaries <- function(V, X) {
  v <- if (is.null(dim(V))) as.numeric(V) else diag(V)
  g <- generalized_sampling_variance(V, X)
  data.frame(
    k = g$k,
    p = g$p,
    diagonal_V = g$diagonal_V,
    arithmetic_mean_vi = mean(v),
    simple_harmonic_mean_vi = 1 / mean(1 / v),
    generalized_v_tilde = g$v_tilde,
    trace_P = g$trace_P,
    stringsAsFactors = FALSE
  )
}

# Independent algebra checks for the trace shortcut and the intercept-only
# Higgins-Thompson expression used by orchaRd's ratio method.
validate_generalized_formula <- function() {
  V_test <- matrix(c(0.20, 0.03, 0.01,
                     0.03, 0.35, 0.02,
                     0.01, 0.02, 0.50), 3, 3, byrow = TRUE)
  X_test <- cbind(1, c(-1, 0, 1))
  W_test <- solve(V_test)
  P_test <- W_test - W_test %*% X_test %*%
    solve(crossprod(X_test, W_test %*% X_test)) %*%
    crossprod(X_test, W_test)
  direct <- (nrow(X_test) - ncol(X_test)) / sum(diag(P_test))
  shortcut <- generalized_sampling_variance(V_test, X_test)$v_tilde
  stopifnot(isTRUE(all.equal(direct, shortcut, tolerance = 1e-12)))
  data.frame(non_diagonal_direct = direct,
             non_diagonal_trace_shortcut = shortcut,
             absolute_difference = abs(direct - shortcut))
}

i2_from_components <- function(components, v_tilde, model, package) {
  stopifnot(length(components) > 0L,
            all(is.finite(components)), all(components >= 0),
            is.finite(v_tilde), v_tilde > 0)
  denominator <- sum(components) + v_tilde
  out <- data.frame(
    package = package,
    model = model,
    component = c("total", names(components)),
    variance = c(sum(components), unname(components)),
    generalized_v_tilde = v_tilde,
    I2_percent = 100 * c(sum(components), unname(components)) / denominator,
    stringsAsFactors = FALSE
  )
  stopifnot(abs(out$I2_percent[out$component == "total"] -
                  sum(out$I2_percent[out$component != "total"])) < 1e-10)
  out
}

read_fit <- function(...) readRDS(file.path(checks_dir, ...))

# ---- Full global Grau-Andres models (metafor) ----
grau_u <- read_fit("totoro_spatial_audit_outputs", "unstructured_only.rds")
grau_s <- read_fit("totoro_spatial_audit_outputs", "spatial_only.rds")
grau_c <- read_fit("totoro_spatial_audit_outputs", "combined.rds")
grau_fits <- list(unstructured_only = grau_u, spatial_only = grau_s, combined = grau_c)

stopifnot(all(vapply(grau_fits, function(x) inherits(x, "rma.mv"), logical(1))),
          length(unique(vapply(grau_fits, function(x) x$k, integer(1)))) == 1L,
          all(vapply(grau_fits, function(x) x$k, integer(1)) == 2361L),
          all(vapply(grau_fits, function(x) isTRUE(all.equal(x$X, grau_u$X)), logical(1))),
          all(vapply(grau_fits, function(x) isTRUE(all.equal(x$vi, grau_u$vi)), logical(1))))

grau_sampling <- sampling_variance_summaries(grau_u$vi, grau_u$X)
grau_w <- 1 / grau_u$vi
grau_ratio_typical <- (grau_u$k - 1) * sum(grau_w) /
  (sum(grau_w)^2 - sum(grau_w^2))
stopifnot(isTRUE(all.equal(grau_ratio_typical,
                          grau_sampling$generalized_v_tilde, tolerance = 1e-12)))
grau_i2 <- rbind(
  i2_from_components(c(effect_size = grau_u$sigma2[1], study = grau_u$sigma2[2]),
                     grau_sampling$generalized_v_tilde, "unstructured_only", "metafor"),
  i2_from_components(c(effect_size = grau_s$sigma2[1], spatial = grau_s$tau2[1]),
                     grau_sampling$generalized_v_tilde, "spatial_only", "metafor"),
  i2_from_components(c(effect_size = grau_c$sigma2[1], study = grau_c$sigma2[2],
                       spatial = grau_c$tau2[1]),
                     grau_sampling$generalized_v_tilde, "combined", "metafor")
)

# ---- Spain matched exponential spatial-only models ----
spain_mf <- read_fit("regional_cross_package_audit_outputs", "metafor_spatial_only.rds")
spain_gt <- read_fit("regional_cross_package_audit_outputs", "glmmTMB_spatial_only.rds")
spain_br <- read_fit("regional_cross_package_audit_outputs", "brms_output", "brms_spatial_only.rds")
spain_prepared <- read_fit("regional_cross_package_audit_outputs", "spain_prepared.rds")

stopifnot(inherits(spain_mf, "rma.mv"), inherits(spain_gt, "glmmTMB"),
          inherits(spain_br, "brmsfit"), spain_mf$k == 186L,
          nrow(spain_prepared$dat) == 186L,
          isTRUE(all.equal(spain_mf$vi, spain_prepared$dat$var_Hedges)))

spain_sampling <- sampling_variance_summaries(spain_mf$V, spain_mf$X)
spain_w <- 1 / spain_mf$vi
spain_ratio_typical <- (spain_mf$k - 1) * sum(spain_w) /
  (sum(spain_w)^2 - sum(spain_w^2))
stopifnot(isTRUE(all.equal(spain_ratio_typical,
                          spain_sampling$generalized_v_tilde, tolerance = 1e-12)))
spain_mf_components <- c(effect_size = spain_mf$sigma2[1], spatial = spain_mf$tau2[1])

spain_gt_spatial <- VarCorr(spain_gt)$cond[["const.1"]]
spain_gt_components <- c(
  effect_size = sigma(spain_gt)^2,
  spatial = unname(attr(spain_gt_spatial, "stddev")[1])^2
)

# These two fits target the same likelihood and should agree numerically.
stopifnot(isTRUE(all.equal(unname(spain_mf_components),
                          unname(spain_gt_components), tolerance = 2e-6)))

spain_i2_frequentist <- rbind(
  i2_from_components(spain_mf_components, spain_sampling$generalized_v_tilde,
                     "Spain_exponential_spatial_only", "metafor"),
  i2_from_components(spain_gt_components, spain_sampling$generalized_v_tilde,
                     "Spain_exponential_spatial_only", "glmmTMB")
)

brms_draws <- posterior::as_draws_df(spain_br)
required_draws <- c("sigma", "sdgp_gpx_kmy_km")
stopifnot(all(required_draws %in% names(brms_draws)))
brms_iid_variance <- brms_draws$sigma^2
brms_spatial_variance <- brms_draws$sdgp_gpx_kmy_km^2
brms_total_variance <- brms_iid_variance + brms_spatial_variance
brms_denominator <- brms_total_variance + spain_sampling$generalized_v_tilde

summarise_draws <- function(x, component, variance_draws) {
  qs_i2 <- quantile(100 * x, c(0.025, 0.5, 0.975), names = FALSE)
  qs_var <- quantile(variance_draws, c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(
    package = "brms",
    model = "Spain_exponential_spatial_only",
    component = component,
    variance_q2.5 = qs_var[1],
    variance_median = qs_var[2],
    variance_q97.5 = qs_var[3],
    generalized_v_tilde = spain_sampling$generalized_v_tilde,
    I2_q2.5_percent = qs_i2[1],
    I2_median_percent = qs_i2[2],
    I2_q97.5_percent = qs_i2[3],
    stringsAsFactors = FALSE
  )
}

spain_i2_brms <- rbind(
  summarise_draws(brms_total_variance / brms_denominator, "total", brms_total_variance),
  summarise_draws(brms_iid_variance / brms_denominator, "effect_size", brms_iid_variance),
  summarise_draws(brms_spatial_variance / brms_denominator, "spatial", brms_spatial_variance)
)

# Draw-wise additivity holds exactly; marginal quantiles need not add because
# quantiles of separate posterior quantities are not algebraically additive.
stopifnot(max(abs((brms_iid_variance + brms_spatial_variance) / brms_denominator -
                    brms_total_variance / brms_denominator)) < 1e-15)

write.csv(grau_sampling, file.path(out_dir, "grau_sampling_variance_summary.csv"), row.names = FALSE)
write.csv(grau_i2, file.path(out_dir, "grau_i2.csv"), row.names = FALSE)
write.csv(spain_sampling, file.path(out_dir, "spain_sampling_variance_summary.csv"), row.names = FALSE)
write.csv(spain_i2_frequentist, file.path(out_dir, "spain_i2_metafor_glmmTMB.csv"), row.names = FALSE)
write.csv(spain_i2_brms, file.path(out_dir, "spain_i2_brms_posterior.csv"), row.names = FALSE)
validation_checks <- validate_generalized_formula()
validation_checks$grau_ratio_minus_generalized <-
  grau_ratio_typical - grau_sampling$generalized_v_tilde
validation_checks$spain_ratio_minus_generalized <-
  spain_ratio_typical - spain_sampling$generalized_v_tilde
write.csv(validation_checks, file.path(out_dir, "formula_validation_checks.csv"), row.names = FALSE)

writeLines(c(
  paste("R:", R.version.string),
  paste("metafor object package version:", grau_s$version),
  paste("glmmTMB:", as.character(packageVersion("glmmTMB"))),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("posterior:", as.character(packageVersion("posterior"))),
  "No model was refitted.",
  "Grau calculations use the saved models' actual vi and X.",
  "Spain calculations use the saved metafor model's V and X for one common v_tilde.",
  "brms I2 was transformed draw-by-draw, then summarized by posterior quantiles."
), file.path(out_dir, "session_info.txt"))

cat("Full Grau sampling variance summary:\n")
print(grau_sampling, digits = 12)
cat("\nFull Grau I2 (%):\n")
print(grau_i2, digits = 12, row.names = FALSE)
cat("\nSpain sampling variance summary:\n")
print(spain_sampling, digits = 12)
cat("\nSpain metafor/glmmTMB I2 (%):\n")
print(spain_i2_frequentist, digits = 12, row.names = FALSE)
cat("\nSpain brms posterior I2 (%):\n")
print(spain_i2_brms, digits = 12, row.names = FALSE)
