# Generalized I2 for the settled Scholer metafor models.
# Reads saved fits only; no model is refitted.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
out_dir <- file.path(project_root, "revision_checks", "scholer_spatial_audit_outputs")

generalized_v_tilde <- function(V, X, tolerance = 1e-12) {
  X <- as.matrix(X)
  k <- nrow(X)
  p <- ncol(X)
  stopifnot(k > p, qr(X)$rank == p)
  if (is.null(dim(V))) {
    v <- as.numeric(V)
    stopifnot(length(v) == k, all(is.finite(v)), all(v > 0))
    diagonal <- TRUE
  } else {
    V <- as.matrix(V)
    stopifnot(identical(dim(V), c(k, k)),
              isTRUE(all.equal(V, t(V), tolerance = tolerance)),
              all(is.finite(V)))
    offdiag <- V
    diag(offdiag) <- 0
    diagonal <- max(abs(offdiag)) <= tolerance
    if (diagonal) {
      v <- diag(V)
      stopifnot(all(v > 0))
    }
  }
  if (diagonal) {
    w <- 1 / v
    WX <- X * w
    trace_P <- sum(w) - sum(diag(solve(crossprod(X, WX), crossprod(WX))))
  } else {
    W <- chol2inv(chol(V))
    WX <- W %*% X
    trace_P <- sum(diag(W)) -
      sum(diag(solve(crossprod(X, WX), crossprod(WX))))
  }
  stopifnot(is.finite(trace_P), trace_P > 0)
  list(v_tilde = (k - p) / trace_P, trace_P = trace_P,
       k = k, p = p, diagonal_V = diagonal)
}

i2_rows <- function(model, components, v_tilde) {
  total <- sum(components)
  denominator <- total + v_tilde
  out <- data.frame(
    model = model,
    component = c("total", names(components)),
    marginal_variance = c(total, unname(components)),
    generalized_v_tilde = v_tilde,
    I2_percent = 100 * c(total, unname(components)) / denominator,
    stringsAsFactors = FALSE
  )
  stopifnot(abs(out$I2_percent[out$component == "total"] -
                  sum(out$I2_percent[out$component != "total"])) < 1e-10)
  out
}

fits <- lapply(c("unstructured_only", "spatial_only", "combined"),
               function(nm) readRDS(file.path(out_dir, paste0(nm, ".rds"))))
names(fits) <- c("unstructured_only", "spatial_only", "combined")
stopifnot(all(vapply(fits, function(x) inherits(x, "rma.mv"), logical(1))),
          all(vapply(fits, function(x) x$k == 949L, logical(1))),
          all(vapply(fits, function(x) x$p == 1L, logical(1))),
          all(vapply(fits, function(x) isTRUE(all.equal(x$vi, fits[[1]]$vi)), logical(1))),
          all(vapply(fits, function(x) isTRUE(all.equal(x$X, fits[[1]]$X)), logical(1))))

sampling <- generalized_v_tilde(fits[[1]]$vi, fits[[1]]$X)
v <- fits[[1]]$vi
w <- 1 / v
ratio_typical <- (length(v) - 1) * sum(w) / (sum(w)^2 - sum(w^2))
stopifnot(isTRUE(all.equal(ratio_typical, sampling$v_tilde, tolerance = 1e-12)))

sampling_summary <- data.frame(
  k = sampling$k, p = sampling$p, diagonal_V = sampling$diagonal_V,
  arithmetic_mean_vi = mean(v),
  simple_harmonic_mean_vi = 1 / mean(1 / v),
  generalized_v_tilde = sampling$v_tilde,
  trace_P = sampling$trace_P,
  ratio_minus_generalized = ratio_typical - sampling$v_tilde,
  stringsAsFactors = FALSE
)

results <- rbind(
  i2_rows("unstructured_only",
          c(effect_size = fits$unstructured_only$sigma2[1],
            study = fits$unstructured_only$sigma2[2]), sampling$v_tilde),
  i2_rows("spatial_only",
          c(effect_size = fits$spatial_only$sigma2[1],
            spatial = fits$spatial_only$tau2[1]), sampling$v_tilde),
  i2_rows("combined",
          c(effect_size = fits$combined$sigma2[1],
            study = fits$combined$sigma2[2],
            spatial = fits$combined$tau2[1]), sampling$v_tilde)
)

write.csv(sampling_summary, file.path(out_dir, "scholer_i2_sampling_variance.csv"),
          row.names = FALSE)
write.csv(results, file.path(out_dir, "scholer_i2_results.csv"), row.names = FALSE)
writeLines(c(
  "I2 uses v_tilde=(k-p)/trace(P), P=W-WX(X'WX)^-1X'W, W=V^-1.",
  "Actual saved model vi and X were used; all three fits have identical vi and X.",
  "Spatial I2 is a marginal variance allocation, not variance explained by distance, a range, or pairwise correlation.",
  "No model was refitted."
), file.path(out_dir, "scholer_i2_readout.txt"))

print(sampling_summary, digits = 12)
print(results, digits = 12, row.names = FALSE)
cat("SCHOLER_I2_AUDIT_COMPLETE\n")
