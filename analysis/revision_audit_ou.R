#!/usr/bin/env Rscript

# Independent audit of the OU example in tutorial_v2.qmd.
#
# Checks:
# 1. 1 - A is the normalised patristic distance for the ultrametric tree.
# 2. The fitted exponential matrix equals ape::corMartins() after rescaling.
# 3. BM and exponential-distance results are reproduced in a clean session.
# 4. AIC for the exponential model includes estimation of the range parameter.

suppressPackageStartupMessages({
  library(ape)
  library(metadat)
  library(metafor)
  library(nlme)
})

dat <- dat.moura2021$dat
dat$species.id.phy <- dat$species.id
dat$effect.size.id <- factor(seq_len(nrow(dat)))
dat$const <- factor(1)
dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)

tree <- compute.brlen(dat.moura2021$tree)
stopifnot(is.ultrametric(tree))

A_bm <- vcv(tree, corr = TRUE)
tip_order <- rownames(A_bm)
patristic <- cophenetic.phylo(tree)[tip_order, tip_order]
tree_height <- max(node.depth.edgelength(tree)[seq_along(tree$tip.label)])

# J - A, where J is an all-ones matrix. This is not I - A if I denotes
# an identity matrix.
D_from_bm <- 1 - A_bm
D_from_tree <- patristic / (2 * tree_height)
distance_error <- max(abs(D_from_bm - D_from_tree))
stopifnot(distance_error < 1e-10)

cat("VERSIONS\n")
cat("R:", R.version.string, "\n")
cat("metafor:", as.character(packageVersion("metafor")), "\n")
cat("metadat:", as.character(packageVersion("metadat")), "\n")
cat("ape:", as.character(packageVersion("ape")), "\n\n")

cat("DISTANCE CHECK\n")
cat("tree_height:", format(tree_height, digits = 12), "\n")
cat("max_abs_error_1_minus_A_vs_patristic:", format(distance_error, scientific = TRUE), "\n\n")

fit_exp <- rma.mv(
  yi,
  vi,
  random = list(
    ~ 1 | study.id,
    ~ 1 | effect.size.id,
    ~ 1 | species.id,
    ~ species.id.phy | const
  ),
  dist = list(species.id.phy = D_from_bm),
  struct = "SPEXP",
  control = list(rho.init = 1),
  data = dat,
  sparse = TRUE,
  method = "REML",
  test = "t"
)

rho_scaled <- unname(fit_exp$rho)
A_exp <- exp(-D_from_bm / rho_scaled)

# corMartins() uses exp(-alpha * patristic_distance).
alpha_patristic <- 1 / (rho_scaled * 2 * tree_height)
martins <- corMartins(
  value = alpha_patristic,
  phy = tree,
  form = ~ species,
  fixed = TRUE
)
martins <- Initialize(martins, data = data.frame(species = tip_order))
A_martins <- corMatrix(martins)
A_martins <- A_martins[tip_order, tip_order]
ou_matrix_error <- max(abs(A_exp - A_martins))
stopifnot(ou_matrix_error < 1e-10)

exp_result <- list(
  estimate = unname(coef(fit_exp)[1]),
  se = fit_exp$se[1],
  ci_lb = fit_exp$ci.lb[1],
  ci_ub = fit_exp$ci.ub[1],
  logLik = as.numeric(logLik(fit_exp)),
  AIC = AIC(fit_exp),
  sigma2 = fit_exp$sigma2,
  tau2 = fit_exp$tau2
)

cat("OU MATRIX CHECK\n")
cat("metafor_rho_on_normalised_distance:", format(rho_scaled, digits = 12), "\n")
cat("equivalent_corMartins_alpha:", format(alpha_patristic, digits = 12), "\n")
cat("max_abs_error_exp_vs_corMartins:", format(ou_matrix_error, scientific = TRUE), "\n\n")

rm(fit_exp, A_exp, A_martins, martins)
gc()

fit_bm <- rma.mv(
  yi,
  vi,
  random = list(
    ~ 1 | study.id,
    ~ 1 | effect.size.id,
    ~ 1 | species.id,
    ~ 1 | species.id.phy
  ),
  R = list(species.id.phy = A_bm),
  data = dat,
  sparse = TRUE,
  method = "REML",
  test = "t"
)

bm_result <- list(
  estimate = unname(coef(fit_bm)[1]),
  se = fit_bm$se[1],
  ci_lb = fit_bm$ci.lb[1],
  ci_ub = fit_bm$ci.ub[1],
  logLik = as.numeric(logLik(fit_bm)),
  AIC = AIC(fit_bm),
  sigma2 = fit_bm$sigma2
)

cat("MODEL RESULTS\n")
results <- data.frame(
  model = c("BM", "exponential_distance"),
  estimate = c(bm_result$estimate, exp_result$estimate),
  se = c(bm_result$se, exp_result$se),
  ci_lb = c(bm_result$ci_lb, exp_result$ci_lb),
  ci_ub = c(bm_result$ci_ub, exp_result$ci_ub),
  logLik = c(bm_result$logLik, exp_result$logLik),
  AIC = c(bm_result$AIC, exp_result$AIC)
)
print(results, row.names = FALSE, digits = 8)

cat("\nVARIANCE COMPONENTS\n")
variance <- data.frame(
  component = c("study", "effect", "species_nonphylo", "species_phylo"),
  BM = bm_result$sigma2,
  exponential_distance = c(exp_result$sigma2, exp_result$tau2)
)
print(variance, row.names = FALSE, digits = 8)

cat("\nAIC NOTE\n")
cat("The exponential-distance AIC above counts the fitted range parameter.\n")
cat("A two-stage fixed-matrix refit must not be compared with BM by AIC unless\n")
cat("the first-stage range parameter is included in the parameter count.\n")
