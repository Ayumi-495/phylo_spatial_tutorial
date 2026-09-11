#!/usr/bin/env Rscript

# Generalised I2 audit for the current Moura BM and OU Model A fits.
#
# This script reads the saved baseline fit objects from the OU correctness
# audit. It does not refit either model. Equations (11)-(13) in the revised
# tutorial define v_tilde from the sampling V and fixed-effect X, and then use
# the common denominator sum(h_j) + v_tilde for total and component I2.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "--run"
root <- normalizePath(".")
out_dir <- file.path(root, "revision_checks", "phylo_generalized_i2_outputs")
baseline_path <- file.path(root, "revision_checks", "ou_correctness_outputs",
                           "baseline_fit_objects.rds")

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
required_outputs <- c("sampling_summary.csv", "generalized_i2.csv",
                      "bm_phylogenetic_heritability.csv", "validation.csv",
                      "readout.txt")

if (identical(mode, "--check")) {
  absent <- required_outputs[!file.exists(file.path(out_dir, required_outputs))]
  assert(!length(absent), paste("Missing I2 audit outputs:", paste(absent, collapse = ", ")))
  sampling <- read.csv(file.path(out_dir, "sampling_summary.csv"), check.names = FALSE)
  i2 <- read.csv(file.path(out_dir, "generalized_i2.csv"), check.names = FALSE)
  h2 <- read.csv(file.path(out_dir, "bm_phylogenetic_heritability.csv"), check.names = FALSE)
  validation <- read.csv(file.path(out_dir, "validation.csv"), check.names = FALSE)
  assert(nrow(sampling) == 1L && isTRUE(sampling$same_sampling_V_and_X[[1L]]),
         "BM and OU Model A do not have the documented common V and X.")
  assert(isTRUE(sampling$same_generalized_v_tilde[[1L]]),
         "BM and OU Model A generalized v_tilde values differ.")
  for (model in unique(i2$model)) {
    x <- i2[i2$model == model, , drop = FALSE]
    total <- x[x$component == "total", , drop = FALSE]
    components <- x[x$component != "total", , drop = FALSE]
    assert(nrow(total) == 1L,
           paste("Expected one total I2 row for", model))
    assert(abs(total$marginal_variance - sum(components$marginal_variance)) < 1e-12,
           paste("Marginal variances do not sum for", model))
    assert(abs(total$I2_percent - sum(components$I2_percent)) < 1e-10,
           paste("Component I2 values do not sum for", model))
  }
  assert(nrow(h2) == 1L && h2$Hphylo2[[1L]] > 0 && h2$Hphylo2[[1L]] < 1,
         "Invalid BM phylogenetic heritability.")
  assert(all(validation$passed), "At least one saved I2 validation check failed.")
  cat("PHYLO_GENERALIZED_I2_AUDIT_CHECKS_PASSED\n")
  quit(status = 0L)
}

assert(identical(mode, "--run"), paste("Unknown mode:", mode))
assert(file.exists(baseline_path), paste("Missing baseline fit object:", baseline_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(baseline_path)
fits <- list(BM = baseline$fit_bm, OU_Model_A = baseline$fit_raw)
assert(all(vapply(fits, inherits, logical(1), what = "rma.mv")),
       "Saved BM and/or OU Model A object is not an rma.mv fit.")

# Equation (11), evaluated directly from V and X. V is diagonal in this
# example, but P is still formed from the full equation as an independent
# check of the trace shortcut used in the spatial audits.
generalized_sampling_variance <- function(vi, X) {
  X <- as.matrix(X)
  k <- length(vi)
  p <- ncol(X)
  assert(k == nrow(X) && k > p && all(is.finite(vi)) && all(vi > 0),
         "Invalid sampling variances or fixed-effect design matrix.")
  assert(qr(X)$rank == p, "Fixed-effect design matrix is rank deficient.")
  V <- diag(vi, nrow = k, ncol = k)
  W <- solve(V)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  trace_P_direct <- sum(diag(P))

  # Algebraically equivalent trace calculation, retained as a numerical check.
  w <- 1 / vi
  XtWX <- crossprod(X, X * w)
  XtW2X <- crossprod(X, X * (w^2))
  trace_P_shortcut <- sum(w) - sum(diag(solve(XtWX, XtW2X)))
  list(k = k, p = p, trace_P_direct = trace_P_direct,
       trace_P_shortcut = trace_P_shortcut,
       v_tilde = (k - p) / trace_P_direct,
       trace_difference = trace_P_direct - trace_P_shortcut,
       diagonal_V = isTRUE(all.equal(V, diag(diag(V)), tolerance = 0)))
}

# The two fits use identical observations, diagonal sampling variances, and
# intercept-only fixed-effect design. Verify this from the fitted objects.
same_vi <- isTRUE(all.equal(fits$BM$vi, fits$OU_Model_A$vi, tolerance = 0))
same_X <- isTRUE(all.equal(fits$BM$X, fits$OU_Model_A$X, tolerance = 0))
same_yi <- isTRUE(all.equal(fits$BM$yi, fits$OU_Model_A$yi, tolerance = 0))
assert(same_vi && same_X && same_yi,
       "BM and OU Model A do not use the same yi, sampling V, and X.")

sampling_bm <- generalized_sampling_variance(fits$BM$vi, fits$BM$X)
sampling_ou <- generalized_sampling_variance(fits$OU_Model_A$vi, fits$OU_Model_A$X)
same_v_tilde <- isTRUE(all.equal(sampling_bm$v_tilde, sampling_ou$v_tilde,
                                 tolerance = 1e-14))
assert(same_v_tilde, "BM and OU Model A generalized v_tilde values differ.")
assert(abs(sampling_bm$trace_difference) < 1e-8,
       "Direct and shortcut trace(P) calculations differ unexpectedly.")

# In both fits, each random-effect covariance has unit diagonal, so h_j in
# equations (12)-(13) is the fitted marginal variance on the diagonal.
components <- list(
  BM = c(study = fits$BM$sigma2[1L],
         effect_size = fits$BM$sigma2[2L],
         species_nonphylogenetic = fits$BM$sigma2[3L],
         species_phylogenetic = fits$BM$sigma2[4L]),
  OU_Model_A = c(study = fits$OU_Model_A$sigma2[1L],
                 effect_size = fits$OU_Model_A$sigma2[2L],
                 species_nonphylogenetic = fits$OU_Model_A$sigma2[3L],
                 species_phylogenetic = fits$OU_Model_A$tau2[1L])
)
assert(all(vapply(components, function(x) all(is.finite(x)) && all(x >= 0), logical(1))),
       "Invalid fitted heterogeneity variance component.")

i2_rows <- function(model, h, v_tilde) {
  total_h <- sum(h)
  denominator <- total_h + v_tilde
  out <- data.frame(
    model = model,
    component = c(names(h), "total"),
    marginal_variance = c(unname(h), total_h),
    generalized_v_tilde = v_tilde,
    denominator = denominator,
    I2_percent = 100 * c(unname(h), total_h) / denominator,
    stringsAsFactors = FALSE
  )
  assert(abs(out$I2_percent[out$component == "total"] -
               sum(out$I2_percent[out$component != "total"])) < 1e-10,
         paste("Component I2 values do not partition total I2 for", model))
  out
}

i2 <- do.call(rbind, Map(i2_rows, names(components), components,
                         MoreArgs = list(v_tilde = sampling_bm$v_tilde)))
rownames(i2) <- NULL

phylo_bm <- components$BM[["species_phylogenetic"]]
nonphylo_bm <- components$BM[["species_nonphylogenetic"]]
hphylo2 <- phylo_bm / (phylo_bm + nonphylo_bm)

sampling_summary <- data.frame(
  k = sampling_bm$k,
  p = sampling_bm$p,
  diagonal_V = sampling_bm$diagonal_V,
  trace_P = sampling_bm$trace_P_direct,
  trace_P_shortcut = sampling_bm$trace_P_shortcut,
  generalized_v_tilde_BM = sampling_bm$v_tilde,
  generalized_v_tilde_OU_Model_A = sampling_ou$v_tilde,
  same_sampling_V_and_X = same_vi && same_X,
  same_generalized_v_tilde = same_v_tilde,
  arithmetic_mean_vi = mean(fits$BM$vi),
  harmonic_mean_vi = length(fits$BM$vi) / sum(1 / fits$BM$vi),
  stringsAsFactors = FALSE
)

h2 <- data.frame(
  model = "BM",
  species_phylogenetic_variance = phylo_bm,
  species_nonphylogenetic_variance = nonphylo_bm,
  among_species_variance = phylo_bm + nonphylo_bm,
  Hphylo2 = hphylo2,
  Hphylo2_percent = 100 * hphylo2,
  interpretation = "Proportion of fitted among-species variance allocated to the phylogenetically structured component",
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c("same_effect_sizes", "same_sampling_variances", "same_fixed_effect_design",
            "diagonal_sampling_V", "direct_vs_shortcut_trace_P", "same_generalized_v_tilde",
            "BM_I2_partitions_total", "OU_Model_A_I2_partitions_total", "BM_Hphylo2_in_unit_interval"),
  value = c(same_yi, same_vi, same_X, sampling_bm$diagonal_V,
            sampling_bm$trace_difference, same_v_tilde,
            i2$I2_percent[i2$model == "BM" & i2$component == "total"] -
              sum(i2$I2_percent[i2$model == "BM" & i2$component != "total"]),
            i2$I2_percent[i2$model == "OU_Model_A" & i2$component == "total"] -
              sum(i2$I2_percent[i2$model == "OU_Model_A" & i2$component != "total"]),
            hphylo2),
  passed = c(same_yi, same_vi, same_X, sampling_bm$diagonal_V,
             abs(sampling_bm$trace_difference) < 1e-8, same_v_tilde,
             abs(i2$I2_percent[i2$model == "BM" & i2$component == "total"] -
                   sum(i2$I2_percent[i2$model == "BM" & i2$component != "total"])) < 1e-10,
             abs(i2$I2_percent[i2$model == "OU_Model_A" & i2$component == "total"] -
                   sum(i2$I2_percent[i2$model == "OU_Model_A" & i2$component != "total"])) < 1e-10,
             hphylo2 > 0 && hphylo2 < 1),
  stringsAsFactors = FALSE
)

write.csv(sampling_summary, file.path(out_dir, "sampling_summary.csv"), row.names = FALSE)
write.csv(i2, file.path(out_dir, "generalized_i2.csv"), row.names = FALSE)
write.csv(h2, file.path(out_dir, "bm_phylogenetic_heritability.csv"), row.names = FALSE)
write.csv(validation, file.path(out_dir, "validation.csv"), row.names = FALSE)
writeLines(c(
  "Generalized representative sampling variance follows equations (11)-(13) in the revised manuscript.",
  "BM and OU Model A use the same yi, diagonal sampling V, and fixed-effect X, so v_tilde is common.",
  "Each h_j is the diagonal marginal variance of its fitted random-effect covariance component.",
  "BM Hphylo2 is phylogenetic species variance divided by total fitted among-species variance; it is not variance explained by phylogeny.",
  "No model was refitted or altered by this script."
), file.path(out_dir, "readout.txt"))

print(sampling_summary, digits = 15, row.names = FALSE)
print(i2, digits = 15, row.names = FALSE)
print(h2, digits = 15, row.names = FALSE)
print(validation, digits = 15, row.names = FALSE)
cat("PHYLO_GENERALIZED_I2_AUDIT_COMPLETE\n")
