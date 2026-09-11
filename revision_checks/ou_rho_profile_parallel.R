#!/usr/bin/env Rscript

# Parallel fixed-rho REML profile for the Moura OU/exponential audit.
# This reads only the baseline fit object produced by ou_correctness_audit.R.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1L]]
}
baseline_path <- arg_value("--baseline", "revision_checks/ou_correctness_outputs/baseline_fit_objects.rds")
out_dir <- arg_value("--out-dir", "revision_checks/ou_correctness_outputs")
workers_requested <- as.integer(arg_value("--workers", "100"))
static_only <- "--static-only" %in% args
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

suppressPackageStartupMessages({ library(ape); library(metadat); library(metafor); library(nlme); library(parallel) })
assert(file.exists(baseline_path), paste("Missing baseline fit object:", baseline_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
profile_dir <- file.path(out_dir, "profile_points")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(baseline_path)
fit_raw <- baseline$fit_raw
rho_raw <- unname(fit_raw$rho)
assert(is.finite(rho_raw) && rho_raw > 0, "Invalid joint rho from baseline.")

# Reconstruct the exact analysis data and verify the saved tree/distance input.
dat <- dat.moura2021$dat
dat$species.id.phy <- dat$species.id
dat$effect.size.id <- factor(seq_len(nrow(dat)))
dat$const <- factor(1)
dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
tip_order <- baseline$tip_order
dat$species.id.phy <- factor(as.character(dat$species.id.phy), levels = tip_order)
assert(!anyNA(dat$species.id.phy) && identical(levels(dat$species.id.phy), tip_order),
       "Model grouping factor differs from saved baseline tip order.")
tree <- compute.brlen(dat.moura2021$tree)
D_raw <- cophenetic.phylo(tree)[tip_order, tip_order]
assert(isTRUE(all.equal(D_raw, baseline$D_raw, tolerance = 1e-12)),
       "Reconstructed direct patristic distance differs from baseline.")

# Matrix and scale checks use the saved joint fit but do not refit any model.
tree_height <- max(node.depth.edgelength(tree)[seq_along(tree$tip.label)])
A_bm <- baseline$A_bm
D_old <- 1 - A_bm
old_distance_error <- max(abs(D_old - D_raw / (2 * tree_height)))
literal_identity_minus_A <- diag(nrow(A_bm)) - A_bm
rho_old <- unname(baseline$fit_old$rho)
A_raw <- exp(-D_raw / rho_raw)
A_old_rescaled <- exp(-D_old / (rho_raw / (2 * tree_height)))
A_old_fitted <- exp(-D_old / rho_old)
martins <- corMartins(value = 1 / rho_raw, phy = tree, form = ~ species, fixed = TRUE)
martins <- Initialize(martins, data = data.frame(species = tip_order))
A_martins <- corMatrix(martins)[tip_order, tip_order]
matrix_checks <- data.frame(
  check = c("tree_is_ultrametric", "tree_height_constructed_branch_length_units",
    "taxa_missing_from_tree", "taxa_missing_from_data", "duplicated_tree_tips", "factor_order_matches_tree",
    "direct_patristic_symmetric_max_abs", "direct_patristic_diagonal_max_abs", "direct_patristic_missing_values",
    "raw_patristic_vs_old_normalised_max_abs", "old_distance_proportionality_constant",
    "literal_identity_minus_A_min_offdiagonal", "direct_vs_old_equivalent_rho_ratio",
    "direct_vs_old_fitted_logLik_abs_diff", "direct_vs_old_equivalent_matrix_max_abs", "direct_vs_old_fitted_matrix_max_abs",
    "corMartins_names_and_order_identical", "corMartins_direct_raw_max_abs", "corMartins_direct_raw_mean_abs"),
  value = c(as.numeric(is.ultrametric(tree)), tree_height, 0, 0, anyDuplicated(tip_order),
    as.numeric(identical(levels(dat$species.id.phy), tip_order)), max(abs(D_raw - t(D_raw))), max(abs(diag(D_raw))),
    sum(is.na(D_raw)), old_distance_error, 1 / (2 * tree_height),
    min(literal_identity_minus_A[row(literal_identity_minus_A) != col(literal_identity_minus_A)]),
    rho_raw / (2 * tree_height * rho_old), abs(as.numeric(logLik(fit_raw)) - as.numeric(logLik(baseline$fit_old))),
    max(abs(A_raw - A_old_rescaled)), max(abs(A_raw - A_old_fitted)),
    as.numeric(identical(rownames(A_martins), rownames(A_raw)) && identical(colnames(A_martins), colnames(A_raw))),
    max(abs(A_raw - A_martins)), mean(abs(A_raw - A_martins))))
metadata <- data.frame(
  field = c("R_version", "metafor_version", "ape_version", "metadat_version", "n_effect_sizes", "n_species",
    "tree_tips", "branch_length_method", "primary_distance", "SPEXP_parameterisation", "rho_units", "alpha_units"),
  value = c(R.version.string, as.character(packageVersion("metafor")), as.character(packageVersion("ape")),
    as.character(packageVersion("metadat")), nrow(dat), nlevels(dat$species.id.phy), length(tip_order),
    "ape::compute.brlen default Grafen topology-based scaling", "raw ape::cophenetic.phylo(tree)", "exp(-d/rho)",
    "constructed Grafen branch-length units", "per constructed Grafen branch-length unit"))
write.csv(matrix_checks, file.path(out_dir, "matrix_equivalence.csv"), row.names = FALSE)
write.csv(metadata, file.path(out_dir, "audit_metadata.csv"), row.names = FALSE)
if (static_only) {
  cat("OU_STATIC_MATRIX_CHECKS_COMPLETED\n")
  quit(status = 0L)
}

# A log-spaced grid provides dense coverage on both sides of the joint REML MLE.
# 0.02--50 times rho spans 3.4 orders of magnitude and contains rho exactly.
multipliers <- exp(seq(log(0.02), log(50), length.out = 101L))
rho_grid <- sort(unique(c(rho_raw, rho_raw * multipliers)))
workers <- max(1L, min(workers_requested, length(rho_grid), detectCores(logical = FALSE) - 2L))

fit_one <- function(rho_value) {
  stem <- formatC(rho_value, format = "f", digits = 12)
  point_path <- file.path(profile_dir, paste0("rho_", stem, ".rds"))
  if (file.exists(point_path)) return(readRDS(point_path))
  result <- tryCatch({
    fit <- rma.mv(yi, vi,
      random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id,
                    ~ species.id.phy | const),
      dist = list(species.id.phy = D_raw), struct = "SPEXP", rho = rho_value,
      control = list(sigma2.init = fit_raw$sigma2, tau2.init = fit_raw$tau2),
      data = dat, sparse = TRUE, method = "REML", test = "t")
    corr <- exp(-D_raw / rho_value)
    cov <- fit$tau2 * corr + diag(fit$sigma2[3L], nrow(corr))
    data.frame(rho = rho_value, alpha = 1 / rho_value,
      REML_logLik = as.numeric(logLik(fit)), AIC_nominal_fixed_rho = AIC(fit),
      pooled_mean = unname(coef(fit)[1L]), pooled_se = fit$se[1L],
      ci_lb = fit$ci.lb[1L], ci_ub = fit$ci.ub[1L],
      study_variance = fit$sigma2[1L], effect_size_variance = fit$sigma2[2L],
      species_nonphylogenetic_variance = fit$sigma2[3L], species_phylogenetic_variance = fit$tau2,
      correlation_offdiag_mean = mean(corr[upper.tri(corr)]),
      covariance_offdiag_mean = mean(cov[upper.tri(cov)]), error = NA_character_)
  }, error = function(e) {
    data.frame(rho = rho_value, alpha = 1 / rho_value, REML_logLik = NA_real_,
      AIC_nominal_fixed_rho = NA_real_, pooled_mean = NA_real_, pooled_se = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_, study_variance = NA_real_, effect_size_variance = NA_real_,
      species_nonphylogenetic_variance = NA_real_, species_phylogenetic_variance = NA_real_,
      correlation_offdiag_mean = NA_real_, covariance_offdiag_mean = NA_real_, error = conditionMessage(e))
  })
  saveRDS(result, point_path)
  result
}

profile <- do.call(rbind, mclapply(rho_grid, fit_one, mc.cores = workers, mc.preschedule = FALSE))
profile <- profile[order(profile$rho), , drop = FALSE]
assert(all(is.finite(profile$REML_logLik)), "At least one rho profile fit failed; inspect profile_points.")
profile$delta_REML_logLik_from_grid_max <- max(profile$REML_logLik) - profile$REML_logLik
profile$within_95pct_LR_grid <- profile$delta_REML_logLik_from_grid_max <= qchisq(0.95, 1) / 2
supported <- profile[profile$within_95pct_LR_grid, , drop = FALSE]
peak <- profile[which.max(profile$REML_logLik), , drop = FALSE]
bracketed <- min(supported$rho) > min(profile$rho) && max(supported$rho) < max(profile$rho)
summary <- data.frame(rho_joint_REML = rho_raw, alpha_joint_REML = 1 / rho_raw,
  profile_grid_peak_rho = peak$rho, profile_grid_peak_logLik = peak$REML_logLik,
  LR95_grid_lower_rho = min(supported$rho), LR95_grid_upper_rho = max(supported$rho),
  interval_bracketed_by_grid = bracketed, n_profile_points = nrow(profile), workers = workers,
  profile_shape = if (nrow(supported) > 0.25 * nrow(profile))
    "broad or weakly identified on evaluated grid" else "comparatively concentrated on evaluated grid",
  profile_boundary_note = if (bracketed) "Grid falls below LR threshold on both sides."
    else "At least one LR bound extends beyond the evaluated grid.")
write.csv(profile, file.path(out_dir, "rho_profile.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "rho_profile_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "rho_profile_session_info.txt"))
cat("OU_RHO_PROFILE_PARALLEL_COMPLETED\n")
