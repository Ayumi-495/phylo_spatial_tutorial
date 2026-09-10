# Generalized I2 for the retained/reportable full Grau-Andres SPGAU candidates.
# Reads saved fits only; no model is refitted.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
out_dir <- file.path(project_root, "revision_checks", "gaussian_global_outputs")
target_dir <- file.path(out_dir, "targeted_free_refits")

v_tilde_from_vi_X <- function(vi, X) {
  X <- as.matrix(X)
  stopifnot(length(vi) == nrow(X), all(is.finite(vi)), all(vi > 0),
            nrow(X) > ncol(X), qr(X)$rank == ncol(X))
  w <- 1 / vi
  WX <- X * w
  trace_P <- sum(w) - sum(diag(solve(crossprod(X, WX), crossprod(WX))))
  list(v_tilde = (nrow(X) - ncol(X)) / trace_P, trace_P = trace_P)
}

i2_rows <- function(model, solution, fit, components, v_tilde) {
  total <- sum(components)
  denom <- total + v_tilde
  out <- data.frame(
    model = model, solution = solution,
    component = c("total", names(components)),
    marginal_variance = c(total, unname(components)),
    generalized_v_tilde = v_tilde,
    I2_percent = 100 * c(total, unname(components)) / denom,
    logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]),
    AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
    rho_km = as.numeric(fit$rho[1]),
    stringsAsFactors = FALSE
  )
  stopifnot(abs(out$I2_percent[out$component == "total"] -
                  sum(out$I2_percent[out$component != "total"])) < 1e-10)
  out
}

spatial <- readRDS(file.path(out_dir, "spatial_only_spgau.rds"))
free_paths <- file.path(target_dir, paste0("free_",
  c("fixed200", "intermediate800", "primary3091"), ".rds"))
if (!all(file.exists(free_paths))) stop("Missing targeted free fits.")
free_fits <- lapply(free_paths, readRDS)
names(free_fits) <- c("fixed200_start", "intermediate800_start", "primary3091_start")

all_fits <- c(list(spatial), free_fits)
stopifnot(all(vapply(all_fits, function(x) inherits(x, "rma.mv"), logical(1))),
          all(vapply(all_fits, function(x) x$k == 2361L, logical(1))),
          all(vapply(all_fits, function(x) isTRUE(all.equal(x$vi, spatial$vi)), logical(1))),
          all(vapply(all_fits, function(x) isTRUE(all.equal(x$X, spatial$X)), logical(1))))

sampling <- v_tilde_from_vi_X(spatial$vi, spatial$X)
stopifnot(isTRUE(all.equal(sampling$v_tilde, 0.111123796928,
                          tolerance = 1e-11)))

results <- i2_rows(
  "spatial_only_spgau", "primary", spatial,
  c(effect_size = spatial$sigma2[1], spatial = spatial$tau2[1]),
  sampling$v_tilde
)
for (nm in names(free_fits)) {
  fit <- free_fits[[nm]]
  results <- rbind(results, i2_rows(
    "combined_spgau", nm, fit,
    c(effect_size = fit$sigma2[1], study = fit$sigma2[2], spatial = fit$tau2[1]),
    sampling$v_tilde
  ))
}

sampling_summary <- data.frame(
  k = spatial$k, p = spatial$p,
  generalized_v_tilde = sampling$v_tilde,
  trace_P = sampling$trace_P,
  verified_target = 0.111123796928,
  difference_from_verified_target = sampling$v_tilde - 0.111123796928,
  stringsAsFactors = FALSE
)
write.csv(sampling_summary, file.path(out_dir, "grau_spgau_i2_sampling_variance.csv"),
          row.names = FALSE)
write.csv(results, file.path(out_dir, "grau_spgau_i2_results.csv"), row.names = FALSE)
writeLines(c(
  "I2 uses the actual saved vi and intercept-only X.",
  "Spatial I2 is a marginal variance allocation, not variance explained by distance, rho, or pairwise correlation.",
  "Combined I2 is retained separately for each targeted free-refit solution so optimizer/ridge dependence is visible.",
  "No model was refitted by this script."
), file.path(out_dir, "grau_spgau_i2_readout.txt"))

print(sampling_summary, digits = 12)
print(results, digits = 12, row.names = FALSE)
cat("GRAU_SPGAU_I2_AUDIT_COMPLETE\n")
