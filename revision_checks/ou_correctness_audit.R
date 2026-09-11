#!/usr/bin/env Rscript

# Reproducible correctness audit for the Moura phylogenetic exponential/OU example.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "--run"
root <- normalizePath(".")
out_dir <- file.path(root, "revision_checks", "ou_correctness_outputs")
audit_path <- file.path(root, "revision_checks", "OU_CORRECTNESS_AUDIT_2026-09-11.md")
required_outputs <- c("audit_metadata.csv", "matrix_equivalence.csv", "model_summary.csv",
                      "variance_components.csv", "correlation_by_distance.csv",
                      "effective_species_covariance.csv", "rho_profile.csv",
                      "rho_profile_summary.csv", "rho_profile_likelihood_interval.csv",
                      "rho_profile_supported_sensitivity.csv", "rho_profile_session_info.txt")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

if (identical(mode, "--check-tutorial")) {
  qmd <- paste(readLines(file.path(root, "tutorial_v2.qmd"), warn = FALSE), collapse = "\n")
  required <- c("cophenetic.phylo(tree)", "exp(-d_{ij}/\\rho)", "332.756", "345.345",
                "constructed Grafen branch-length units", "confidence interval", "profile")
  missing <- required[!vapply(required, grepl, logical(1), x = qmd, fixed = TRUE)]
  assert(!length(missing), paste("Missing tutorial content:", paste(missing, collapse = "; ")))
  prohibited <- c("indicating that phylogenetic similarity in assortative mating strength decays rapidly",
                  "pulled towards their own lineage-specific optima", "limited deep phylogenetic inertia")
  found <- prohibited[vapply(prohibited, grepl, logical(1), x = qmd, fixed = TRUE)]
  assert(!length(found), paste("Retained prohibited interpretation:", paste(found, collapse = "; ")))
  cat("OU_TUTORIAL_CHECKS_PASSED\n")
  quit(status = 0L)
}

if (identical(mode, "--check-executing-render")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result <- system2("quarto", c("render", "tutorial_v2.qmd", "--to", "html"), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  writeLines(c(paste("exit_status:", status), result), file.path(out_dir, "executing_render.log"))
  assert(status != 0L, "Executing render unexpectedly succeeded.")
  assert(any(grepl("moura2021_BM_meta_reg.rds", result, fixed = TRUE)),
         "Executing render did not reach the documented missing RDS boundary.")
  cat("EXPECTED_EXECUTING_RENDER_BLOCK_CONFIRMED\n")
  quit(status = 0L)
}

if (identical(mode, "--check")) {
  assert(file.exists(audit_path), "Missing OU audit Markdown.")
  absent <- required_outputs[!file.exists(file.path(out_dir, required_outputs))]
  assert(!length(absent), paste("Missing audit outputs:", paste(absent, collapse = ", ")))
  matrices <- read.csv(file.path(out_dir, "matrix_equivalence.csv"), check.names = FALSE)
  models <- read.csv(file.path(out_dir, "model_summary.csv"), check.names = FALSE)
  profile <- read.csv(file.path(out_dir, "rho_profile.csv"), check.names = FALSE)
  interval <- read.csv(file.path(out_dir, "rho_profile_likelihood_interval.csv"), check.names = FALSE)
  value <- function(name) matrices$value[matrices$check == name][[1L]]
  assert(value("raw_patristic_vs_old_normalised_max_abs") < 1e-10, "Incorrect distance equivalence.")
  assert(value("corMartins_direct_raw_max_abs") < 1e-10, "corMartins mismatch.")
  assert(models$estimated_parameters[models$model == "BM"] == 5L, "Incorrect BM AIC parameter count.")
  assert(models$estimated_parameters[models$model == "OU_joint_raw"] == 6L, "Incorrect joint OU AIC parameter count.")
  assert(nrow(profile) >= 10L && all(is.finite(profile$REML_logLik)), "Incomplete rho profile.")
  assert(isTRUE(interval$optimum_interior[[1L]]) && interval$lower_rho_log_interpolated[[1L]] > 0,
         "Invalid likelihood-profile interval.")
  cat("OU_CORRECTNESS_AUDIT_CHECKS_PASSED\n")
  quit(status = 0L)
}

assert(mode %in% c("--run", "--baseline"), paste("Unknown mode:", mode))
suppressPackageStartupMessages({ library(ape); library(metadat); library(metafor); library(nlme) })
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Exact tutorial data path and constructed Grafen tree.
dat <- dat.moura2021$dat
dat$species.id.phy <- dat$species.id
dat$effect.size.id <- factor(seq_len(nrow(dat)))
dat$const <- factor(1)
dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
tree <- compute.brlen(dat.moura2021$tree)
assert(is.ultrametric(tree), "The constructed tree is not ultrametric.")
tree_height <- max(node.depth.edgelength(tree)[seq_along(tree$tip.label)])
A_bm <- vcv(tree, corr = TRUE)
tip_order <- rownames(A_bm)
assert(identical(tip_order, colnames(A_bm)), "BM matrix names differ by margin.")

# Make grouping-factor and matrix ordering explicit.
observed_taxa <- as.character(dat$species.id.phy)
missing_from_tree <- setdiff(unique(observed_taxa), tip_order)
missing_from_data <- setdiff(tip_order, unique(observed_taxa))
assert(!anyDuplicated(tip_order), "Duplicated tree tips.")
assert(!length(missing_from_tree) && !length(missing_from_data), "Tree/data taxa differ.")
dat$species.id.phy <- factor(observed_taxa, levels = tip_order)
assert(identical(levels(dat$species.id.phy), tip_order) && !anyNA(dat$species.id.phy),
       "Grouping factor does not match tree order.")

D_raw <- cophenetic.phylo(tree)[tip_order, tip_order]
assert(identical(rownames(D_raw), tip_order) && identical(colnames(D_raw), tip_order), "Distance order mismatch.")
assert(isTRUE(all.equal(D_raw, t(D_raw), tolerance = 0)) && all(diag(D_raw) == 0) && !anyNA(D_raw),
       "Invalid direct patristic matrix.")
D_old <- 1 - A_bm # J - A: J is all ones, not an identity matrix.
old_distance_error <- max(abs(D_old - D_raw / (2 * tree_height)))
assert(old_distance_error < 1e-10, "Old distance is not d/(2h).")
literal_identity_minus_A <- diag(nrow(A_bm)) - A_bm

fit_bm <- rma.mv(yi, vi,
  random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
  R = list(species.id.phy = A_bm), data = dat, sparse = TRUE, method = "REML", test = "t")

# Model A: joint estimation on raw patristic/Grafen distance units.
fit_raw <- rma.mv(yi, vi,
  random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ species.id.phy | const),
  dist = list(species.id.phy = D_raw), struct = "SPEXP", control = list(rho.init = 0.04),
  data = dat, sparse = TRUE, method = "REML", test = "t")
rho_raw <- unname(fit_raw$rho)
alpha_raw <- 1 / rho_raw

# Old normalised implementation, fitted only to establish the special equivalence.
fit_old <- rma.mv(yi, vi,
  random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ species.id.phy | const),
  dist = list(species.id.phy = D_old), struct = "SPEXP", control = list(rho.init = 0.02),
  data = dat, sparse = TRUE, method = "REML", test = "t")
rho_old <- unname(fit_old$rho)

A_raw <- exp(-D_raw / rho_raw)
A_old_rescaled <- exp(-D_old / (rho_raw / (2 * tree_height)))
A_old_fitted <- exp(-D_old / rho_old)
martins <- corMartins(value = alpha_raw, phy = tree, form = ~ species, fixed = TRUE)
martins <- Initialize(martins, data = data.frame(species = tip_order))
A_martins <- corMatrix(martins)[tip_order, tip_order]
martins_max <- max(abs(A_raw - A_martins))
martins_mean <- mean(abs(A_raw - A_martins))
assert(martins_max < 1e-10, "Direct exponential does not match corMartins.")

# Model B: fixed matrix refit. Its nominal AIC excludes the prior rho estimate.
fit_fixed <- rma.mv(yi, vi,
  random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
  R = list(species.id.phy = A_raw), data = dat, sparse = TRUE, method = "REML", test = "t")

model_row <- function(name, fit, rho = NA_real_, alpha = NA_real_, note = "") {
  ll <- as.numeric(logLik(fit)); aic <- AIC(fit)
  data.frame(model = name, pooled_mean = unname(coef(fit)[1L]), pooled_se = fit$se[1L],
             ci_lb = fit$ci.lb[1L], ci_ub = fit$ci.ub[1L], REML_logLik = ll,
             estimated_parameters = as.integer(round((aic + 2 * ll) / 2)), AIC = aic,
             rho = rho, alpha = alpha, note = note, stringsAsFactors = FALSE)
}
models <- rbind(
  model_row("BM", fit_bm, note = "vcv(tree, corr = TRUE)"),
  model_row("OU_joint_raw", fit_raw, rho_raw, alpha_raw, "Model A: raw direct patristic distance; rho jointly estimated"),
  model_row("OU_joint_old_normalised", fit_old, rho_old, 1 / rho_old, "Diagnostic only: J-A = d/(2h) for this tree"),
  model_row("OU_fixed_matrix", fit_fixed, rho_raw, alpha_raw, "Model B: matrix fixed after Model A; nominal AIC omits rho"))
models$delta_AIC_vs_BM <- models$AIC - models$AIC[models$model == "BM"]

variances <- data.frame(component = c("study", "effect_size", "species_nonphylogenetic", "species_phylogenetic"),
  BM = fit_bm$sigma2,
  OU_joint_raw = c(fit_raw$sigma2, fit_raw$tau2),
  OU_joint_old_normalised = c(fit_old$sigma2, fit_old$tau2),
  OU_fixed_matrix = fit_fixed$sigma2, stringsAsFactors = FALSE)

off_diag <- upper.tri(D_raw)
distances <- unique(c(0, quantile(D_raw[off_diag], c(0.05, 0.25, 0.5, 0.75, 0.95)), max(D_raw)))
correlation_by_distance <- data.frame(patristic_distance = as.numeric(distances),
  BM_correlation = pmax(0, 1 - distances / (2 * tree_height)), OU_correlation = exp(-distances / rho_raw))

cov_summary <- function(name, correlation, fit, phylo, nonphylo) {
  cov <- phylo * correlation + diag(nonphylo, nrow(correlation))
  off <- upper.tri(cov)
  data.frame(model = name, phylogenetic_variance = phylo, nonphylogenetic_species_variance = nonphylo,
    correlation_offdiag_mean = mean(correlation[off]), correlation_offdiag_median = median(correlation[off]),
    covariance_offdiag_mean = mean(cov[off]), covariance_offdiag_median = median(cov[off]),
    covariance_offdiag_q05 = unname(quantile(cov[off], 0.05)), covariance_offdiag_q95 = unname(quantile(cov[off], 0.95)),
    pooled_se = fit$se[1L])
}
covariances <- rbind(cov_summary("BM", A_bm, fit_bm, fit_bm$sigma2[4L], fit_bm$sigma2[3L]),
                     cov_summary("OU_joint_raw", A_raw, fit_raw, fit_raw$tau2, fit_raw$sigma2[3L]))

# Preserve the completed model fits before starting the expensive profile.
write.csv(models, file.path(out_dir, "model_summary.csv"), row.names = FALSE)
write.csv(variances, file.path(out_dir, "variance_components.csv"), row.names = FALSE)
write.csv(correlation_by_distance, file.path(out_dir, "correlation_by_distance.csv"), row.names = FALSE)
write.csv(covariances, file.path(out_dir, "effective_species_covariance.csv"), row.names = FALSE)
saveRDS(list(tree = tree, tip_order = tip_order, D_raw = D_raw, A_bm = A_bm,
             fit_bm = fit_bm, fit_raw = fit_raw, fit_old = fit_old, fit_fixed = fit_fixed),
        file.path(out_dir, "baseline_fit_objects.rds"))
if (identical(mode, "--baseline")) {
  cat("OU_CORRECTNESS_BASELINE_COMPLETED\n")
  quit(status = 0L)
}

# Fixed-rho profile: all remaining components are re-estimated for each point.
# Each point is cached as it completes, so an expensive fit is never lost.
multipliers <- c(0.20, 0.35, 0.50, 0.70, 0.85, 1, 1.15, 1.40, 2, 3, 5, 10, 20)
rho_grid <- sort(unique(rho_raw * multipliers))
profile_dir <- file.path(out_dir, "profile_points")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
profile_rows <- lapply(rho_grid, function(rho_value) {
  point_path <- file.path(profile_dir, paste0("rho_", formatC(rho_value, format = "f", digits = 10), ".csv"))
  if (file.exists(point_path)) return(read.csv(point_path, check.names = FALSE))
  fit <- rma.mv(yi, vi,
    random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ species.id.phy | const),
    dist = list(species.id.phy = D_raw), struct = "SPEXP", rho = rho_value,
    control = list(sigma2.init = fit_raw$sigma2, tau2.init = fit_raw$tau2),
    data = dat, sparse = TRUE, method = "REML", test = "t")
  corr <- exp(-D_raw / rho_value); cov <- fit$tau2 * corr + diag(fit$sigma2[3L], nrow(corr))
  row <- data.frame(rho = rho_value, alpha = 1 / rho_value, REML_logLik = as.numeric(logLik(fit)),
    AIC_nominal_fixed_rho = AIC(fit), pooled_mean = unname(coef(fit)[1L]), pooled_se = fit$se[1L],
    ci_lb = fit$ci.lb[1L], ci_ub = fit$ci.ub[1L], study_variance = fit$sigma2[1L],
    effect_size_variance = fit$sigma2[2L], species_nonphylogenetic_variance = fit$sigma2[3L],
    species_phylogenetic_variance = fit$tau2, correlation_offdiag_mean = mean(corr[upper.tri(corr)]),
    covariance_offdiag_mean = mean(cov[upper.tri(cov)]))
  write.csv(row, point_path, row.names = FALSE)
  row
})
profile <- do.call(rbind, profile_rows)
profile$delta_REML_logLik_from_grid_max <- max(profile$REML_logLik) - profile$REML_logLik
profile$within_95pct_LR_grid <- profile$delta_REML_logLik_from_grid_max <= qchisq(0.95, 1) / 2
supported <- profile[profile$within_95pct_LR_grid, , drop = FALSE]
peak <- profile[which.max(profile$REML_logLik), , drop = FALSE]
bracketed <- min(supported$rho) > min(profile$rho) && max(supported$rho) < max(profile$rho)
profile_summary <- data.frame(rho_joint_REML = rho_raw, alpha_joint_REML = alpha_raw,
  profile_grid_peak_rho = peak$rho, profile_grid_peak_logLik = peak$REML_logLik,
  LR95_grid_lower_rho = min(supported$rho), LR95_grid_upper_rho = max(supported$rho),
  interval_bracketed_by_grid = bracketed,
  profile_shape = if (nrow(supported) > 4L) "broad or weakly identified on evaluated grid" else "comparatively concentrated on evaluated grid",
  profile_boundary_note = if (bracketed) "Grid falls below LR threshold on both sides." else "At least one LR bound extends beyond the evaluated grid.")

matrix_checks <- data.frame(check = c("tree_is_ultrametric", "tree_height_constructed_branch_length_units",
  "taxa_missing_from_tree", "taxa_missing_from_data", "duplicated_tree_tips", "factor_order_matches_tree",
  "direct_patristic_symmetric_max_abs", "direct_patristic_diagonal_max_abs", "direct_patristic_missing_values",
  "raw_patristic_vs_old_normalised_max_abs", "old_distance_proportionality_constant",
  "literal_identity_minus_A_min_offdiagonal", "direct_vs_old_equivalent_rho_ratio",
  "direct_vs_old_fitted_logLik_abs_diff", "direct_vs_old_equivalent_matrix_max_abs", "direct_vs_old_fitted_matrix_max_abs",
  "corMartins_names_and_order_identical", "corMartins_direct_raw_max_abs", "corMartins_direct_raw_mean_abs"),
  value = c(as.numeric(is.ultrametric(tree)), tree_height, length(missing_from_tree), length(missing_from_data),
  anyDuplicated(tip_order), as.numeric(identical(levels(dat$species.id.phy), tip_order)), max(abs(D_raw - t(D_raw))),
  max(abs(diag(D_raw))), sum(is.na(D_raw)), old_distance_error, 1 / (2 * tree_height),
  min(literal_identity_minus_A[row(literal_identity_minus_A) != col(literal_identity_minus_A)]),
  rho_raw / (2 * tree_height * rho_old), abs(as.numeric(logLik(fit_raw)) - as.numeric(logLik(fit_old))),
  max(abs(A_raw - A_old_rescaled)), max(abs(A_raw - A_old_fitted)),
  as.numeric(identical(rownames(A_martins), rownames(A_raw)) && identical(colnames(A_martins), colnames(A_raw))),
  martins_max, martins_mean), stringsAsFactors = FALSE)

metadata <- data.frame(field = c("R_version", "metafor_version", "ape_version", "metadat_version", "n_effect_sizes",
  "n_species", "tree_tips", "branch_length_method", "primary_distance", "SPEXP_parameterisation", "rho_units", "alpha_units"),
  value = c(R.version.string, as.character(packageVersion("metafor")), as.character(packageVersion("ape")),
  as.character(packageVersion("metadat")), nrow(dat), nlevels(dat$species.id.phy), length(tip_order),
  "ape::compute.brlen default Grafen topology-based scaling", "raw ape::cophenetic.phylo(tree)", "exp(-d/rho)",
  "constructed Grafen branch-length units", "per constructed Grafen branch-length unit"))

write.csv(metadata, file.path(out_dir, "audit_metadata.csv"), row.names = FALSE)
write.csv(matrix_checks, file.path(out_dir, "matrix_equivalence.csv"), row.names = FALSE)
write.csv(models, file.path(out_dir, "model_summary.csv"), row.names = FALSE)
write.csv(variances, file.path(out_dir, "variance_components.csv"), row.names = FALSE)
write.csv(correlation_by_distance, file.path(out_dir, "correlation_by_distance.csv"), row.names = FALSE)
write.csv(covariances, file.path(out_dir, "effective_species_covariance.csv"), row.names = FALSE)
write.csv(profile, file.path(out_dir, "rho_profile.csv"), row.names = FALSE)
write.csv(profile_summary, file.path(out_dir, "rho_profile_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("OU_CORRECTNESS_AUDIT_RUN_COMPLETED\n")
