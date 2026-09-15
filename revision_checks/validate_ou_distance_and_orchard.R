#!/usr/bin/env Rscript

# Regression checks for the restricted BM-correlation identity and the released
# orchaRd 2.2.1 adapter demonstration used in tutorial_v2.qmd. This script
# deliberately reconstructs patristic distance from the tree, never from A.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "all"
write_plot <- "--write-plot" %in% args
if (!mode %in% c("--mode", "ou", "conversions", "all")) {
  stop("Use --mode ou, --mode conversions, or --mode all.", call. = FALSE)
}
if (identical(mode, "--mode")) {
  if (length(args) < 2L || !args[[2L]] %in% c("ou", "conversions", "all")) {
    stop("Use --mode ou, --mode conversions, or --mode all.", call. = FALSE)
  }
  mode <- args[[2L]]
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_ou_distance <- function() {
  dat <- metadat::dat.moura2021
  tree <- ape::compute.brlen(dat$tree)
  tip_order <- tree$tip.label
  A <- ape::vcv(tree, corr = TRUE)[tip_order, tip_order]
  D <- ape::cophenetic.phylo(tree)[tip_order, tip_order]
  h <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
  J <- matrix(1, nrow(A), ncol(A), dimnames = dimnames(A))
  I <- diag(nrow(A))
  off_diagonal <- upper.tri(A)

  assert_true(ape::is.ultrametric(tree), "The regression example must be ultrametric.")
  assert_true(is.finite(h) && h > 0, "Tree height must be positive and finite.")
  assert_true(isTRUE(all.equal(unname(diag(A)), rep(1, nrow(A)), tolerance = 1e-12)),
              "A must be a unit-diagonal BM correlation matrix.")
  assert_true(identical(rownames(A), rownames(D)) && identical(colnames(A), colnames(D)),
              "A and direct patristic D must have the same ordered tips.")

  identity_error <- max(abs((J - A) - D / (2 * h)))
  I_minus_A <- I - A
  min_i_minus_a <- min(I_minus_A[off_diagonal])
  min_direct_distance <- min(D[off_diagonal])

  assert_true(identity_error < 1e-10,
              "J - A did not equal direct D/(2h) for the normalized ultrametric example.")
  assert_true(min_i_minus_a < -1e-10,
              "Literal I - A did not expose a negative off-diagonal non-distance entry.")
  assert_true(min_direct_distance >= -1e-12,
              "Direct patristic distances must be non-negative off the diagonal.")

  message(sprintf("OU normalized-ultrametric max error: %.3e", identity_error))
  message(sprintf("OU literal-I-minus-A minimum off-diagonal: %.9f", min_i_minus_a))
  message("OU_NORMALIZED_ULTRAMETRIC_IDENTITY_PASSED")
  message("OU_LITERAL_I_MINUS_A_NEGATIVE_CONTROL_PASSED")
  message("OU_DISTANCE_REGRESSION_PASSED")
  invisible(list(identity_error = identity_error, min_i_minus_a = min_i_minus_a))
}

lim_data <- function() {
  data("dat.lim2014", package = "metadat", envir = environment())
  dat <- metadat::dat.lim2014$o_o_unadj
  tree <- ape::compute.brlen(metadat::dat.lim2014$o_o_unadj_tree)
  dat <- metafor::escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  dat$species <- factor(as.character(dat$species))
  dat$phy <- factor(as.character(dat$species), levels = sort(tree$tip.label))
  assert_true(!anyNA(dat$phy), "Every Lim species must match the tree.")
  dat$id <- factor(seq_len(nrow(dat)))
  dat$g <- factor("all")
  A <- ape::vcv(tree, corr = TRUE)
  A <- A[levels(dat$phy), levels(dat$phy)]
  V <- diag(dat$vi, nrow = nrow(dat), ncol = nrow(dat))
  rownames(V) <- colnames(V) <- levels(dat$id)
  list(data = dat, A = A, V = V)
}

intervals_from_model <- function(beta, vb) {
  se <- sqrt(diag(vb))
  cbind(estimate = as.numeric(beta),
        ci.lb = as.numeric(beta) - stats::qnorm(0.975) * se,
        ci.ub = as.numeric(beta) + stats::qnorm(0.975) * se)
}

check_conversions <- function(write_plot = FALSE) {
  library(glmmTMB)
  inputs <- lim_data()
  dat <- inputs$data
  A <- inputs$A
  V <- inputs$V
  fit_tmb <- glmmTMB::glmmTMB(
    yi ~ 1 + environment +
      equalto(0 + id | g, V) +
      (1 | species) +
      propto(0 + phy | g, A),
    data = dat, REML = TRUE
  )
  fit_meta <- metafor::rma.mv(
    yi, vi, mods = ~ environment,
    random = list(~ 1 | id, ~ 1 | phy, ~ 1 | species),
    R = list(phy = A), data = dat, method = "REML", sparse = TRUE
  )
  converted <- orchaRd::glmmTMB_to_rma(
    fit_tmb, yi = "yi", vi = "vi", data = dat, measure = "GEN", test = "z"
  )

  tmb_beta <- glmmTMB::fixef(fit_tmb)$cond
  converted_beta <- as.numeric(converted$b)
  names(converted_beta) <- rownames(converted$b)
  tmb_ci <- intervals_from_model(tmb_beta, stats::vcov(fit_tmb)$cond)
  converted_ci <- intervals_from_model(converted$b[, 1L], converted$vb)
  meta_ci <- intervals_from_model(fit_meta$b[, 1L], fit_meta$vb)

  assert_true(identical(unname(converted_beta), unname(as.numeric(tmb_beta))),
              "Converted fixed effects differ from the original glmmTMB fixed effects.")
  assert_true(isTRUE(all.equal(unname(tmb_ci), unname(converted_ci), tolerance = 1e-12)),
              "Converted Wald confidence intervals differ from original glmmTMB intervals.")
  yi_difference <- max(abs(unname(converted$yi) - unname(dat$yi)))
  message(sprintf("Lim yi retention maximum absolute difference: %.17g", yi_difference))
  assert_true(yi_difference < 1e-15,
              "Converted yi differs from the supplied Lim effect-size vector.")
  vi_difference <- max(abs(unname(converted$vi) - unname(dat$vi)))
  message(sprintf("Lim vi retention maximum absolute difference: %.17g", vi_difference))
  assert_true(vi_difference < 1e-15,
              "Converted vi differs from the supplied Lim sampling-variance vector.")
  assert_true(identical(as.character(converted$data$id), as.character(dat$id)),
              "Converted data rows or id grouping are not retained in input order.")
  assert_true(identical(rownames(converted$X), rownames(stats::model.matrix(~ environment, dat))),
              "Converted fixed-effect model-matrix rows are not retained in input order.")
  assert_true(max(abs(tmb_ci - meta_ci)) < 0.005,
              "Matched Lim metafor and glmmTMB coefficient intervals diverged beyond tolerance.")

  orchard_tmb <- orchaRd::mod_results(converted, mod = "environment", group = "id")
  orchard_meta <- orchaRd::mod_results(fit_meta, mod = "environment", group = "id")
  assert_true(identical(names(orchard_tmb), names(orchard_meta)),
              "Converted and metafor orchard summaries have incompatible columns.")
  orchard_table_tmb <- orchard_tmb$mod_table
  orchard_table_meta <- orchard_meta$mod_table
  message("Lim orchard summary columns: ", paste(names(orchard_table_tmb), collapse = ", "))
  pi_columns <- c("lowerPR", "upperPR")
  assert_true(all(pi_columns %in% names(orchard_table_tmb)),
              "orchaRd summary did not expose lowerPR and upperPR prediction-interval endpoints.")
  assert_true(max(abs(as.matrix(orchard_table_tmb[, pi_columns, drop = FALSE]) -
                      as.matrix(orchard_table_meta[, pi_columns, drop = FALSE]))) < 0.005,
              "Converted and matched metafor orchard prediction intervals diverged beyond tolerance.")

  if (isTRUE(write_plot)) {
    plot <- orchaRd::orchard_plot(
      orchard_tmb, mod = "environment", group = "id", xlab = "Effect size (Fisher's Z)",
      angle = 45, g = FALSE
    ) + ggplot2::theme_classic()
    ggplot2::ggsave("figs/tutorial/lim_glmmtmb_orchard.png", plot,
                    width = 7, height = 4.8, units = "in", dpi = 300)
    message("GLMMTMB_ORCHARD_PLOT_WRITTEN")
  }

  message(sprintf("Lim converted-versus-metafor maximum coefficient/CI difference: %.6f",
                  max(abs(tmb_ci - meta_ci))))
  message(sprintf("Lim converted-versus-metafor maximum orchard-PI difference: %.6f",
                  max(abs(as.matrix(orchard_table_tmb[, pi_columns, drop = FALSE]) -
                          as.matrix(orchard_table_meta[, pi_columns, drop = FALSE])))))
  message("GLMMTMB_ORCHARD_CONVERSION_VALIDATION_PASSED")
  invisible(list(tmb = fit_tmb, metafor = fit_meta, converted = converted,
                 orchard_tmb = orchard_tmb, orchard_meta = orchard_meta))
}

if (mode %in% c("ou", "all")) check_ou_distance()
if (mode %in% c("conversions", "all")) check_conversions(write_plot)
