# Generalized I2 for the saved published_cleaned SPEXP primary fits.
# No model is fitted or updated by this script.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
r18_dir <- file.path(project_root, "revision_checks", "reviewer18_influential_effects_outputs")
out_dir <- file.path(project_root, "revision_checks", "reviewer18_cleaned_primary_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(file.path(r18_dir, "published_cleaned_prepared.rds"))
fit_paths <- file.path(r18_dir, paste0(c("unstructured_only", "spatial_only", "combined"), ".rds"))
fits <- lapply(fit_paths, readRDS)
names(fits) <- c("unstructured_only", "spatial_only", "combined")

stopifnot(nrow(prepared$dat) == 2355L,
          all(vapply(fits, inherits, logical(1), what = "rma.mv")),
          all(vapply(fits, function(x) x$k == 2355L, logical(1))),
          all(vapply(fits, function(x) isTRUE(all.equal(as.numeric(x$vi),
                                                       prepared$dat$var_Hedges)), logical(1))),
          all(vapply(fits, function(x) isTRUE(all.equal(x$X, fits[[1L]]$X)), logical(1))))

v_tilde_from_vi_X <- function(vi, X) {
  X <- as.matrix(X)
  stopifnot(length(vi) == nrow(X), all(is.finite(vi)), all(vi > 0),
            nrow(X) > ncol(X), qr(X)$rank == ncol(X))
  w <- 1 / vi
  XtWX <- crossprod(X, X * w)
  XtW2X <- crossprod(X, X * (w^2))
  trace_P <- sum(w) - sum(diag(solve(XtWX, XtW2X)))
  list(k = nrow(X), p = ncol(X), trace_P = trace_P,
       v_tilde = (nrow(X) - ncol(X)) / trace_P,
       XtWX = XtWX, XtW2X = XtW2X)
}

sampling <- v_tilde_from_vi_X(fits[[1L]]$vi, fits[[1L]]$X)
stopifnot(is.finite(sampling$v_tilde), sampling$v_tilde > 0,
          is.finite(sampling$trace_P), sampling$trace_P > 0)

component_map <- list(
  unstructured_only = c(effect_size = fits$unstructured_only$sigma2[1],
                        study = fits$unstructured_only$sigma2[2]),
  spatial_only = c(effect_size = fits$spatial_only$sigma2[1],
                   spatial = fits$spatial_only$tau2[1]),
  combined = c(effect_size = fits$combined$sigma2[1],
               study = fits$combined$sigma2[2],
               spatial = fits$combined$tau2[1])
)

i2_rows <- function(model, components, v_tilde) {
  total_variance <- sum(components)
  denominator <- total_variance + v_tilde
  out <- data.frame(
    dataset = "published_cleaned",
    model = model,
    component = c(names(components), "total"),
    variance = c(unname(components), total_variance),
    generalized_v_tilde = v_tilde,
    denominator = denominator,
    I2_percent = 100 * c(unname(components), total_variance) / denominator,
    stringsAsFactors = FALSE
  )
  stopifnot(abs(out$I2_percent[out$component == "total"] -
                  sum(out$I2_percent[out$component != "total"])) < 1e-10)
  out
}

i2 <- do.call(rbind, Map(i2_rows, names(component_map), component_map,
                         MoreArgs = list(v_tilde = sampling$v_tilde)))
rownames(i2) <- NULL

sampling_summary <- data.frame(
  dataset = "published_cleaned",
  k = sampling$k,
  p = sampling$p,
  diagonal_V = TRUE,
  trace_P = sampling$trace_P,
  generalized_v_tilde = sampling$v_tilde,
  arithmetic_mean_vi = mean(fits[[1L]]$vi),
  harmonic_mean_vi = length(fits[[1L]]$vi) / sum(1 / fits[[1L]]$vi),
  stringsAsFactors = FALSE
)

write.csv(sampling_summary, file.path(out_dir, "cleaned_sampling_variance_summary.csv"),
          row.names = FALSE)
write.csv(i2, file.path(out_dir, "cleaned_generalized_i2.csv"), row.names = FALSE)
writeLines(c(
  "Generalized v_tilde uses the saved cleaned models' actual vi and X.",
  "trace(P) is evaluated as trace(W) - trace((X'WX)^(-1) X'W^2X), without materializing dense W or P.",
  "Component I2 values partition marginal heterogeneity variance; spatial I2 is not variance explained by distance or rho.",
  "No model was fitted or updated."
), file.path(out_dir, "cleaned_i2_readout.txt"))

print(sampling_summary, digits = 12, row.names = FALSE)
print(i2, digits = 12, row.names = FALSE)
cat("CLEANED_I2_COMPLETE\n")
