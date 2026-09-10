args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript reviewer18_compile_comparison.R <full_results_csv> <cleaned_output_dir>")
}

full_path <- args[[1L]]
out_dir <- args[[2L]]
clean_path <- file.path(out_dir, "published_cleaned_primary_results.csv")
stopifnot(file.exists(full_path), file.exists(clean_path), dir.exists(out_dir))

full <- read.csv(full_path, stringsAsFactors = FALSE)
clean <- read.csv(clean_path, stringsAsFactors = FALSE)
expected_models <- c("unstructured_only", "spatial_only", "combined")
stopifnot(setequal(full$model, expected_models), setequal(clean$model, expected_models))

full_h <- data.frame(
  dataset = "all_spatially_usable",
  n_effects = 2361L,
  n_studies = 393L,
  n_sites = 383L,
  model = full$model,
  mean = full$mean,
  ci_lb = full$ci_lb,
  ci_ub = full$ci_ub,
  effect_variance = full$sigma2_effect,
  study_variance = full$sigma2_study,
  spatial_variance = full$tau2_spatial,
  rho_km = full$rho_km,
  logLik_REML = full$logLik_REML,
  AIC_REML = full$AIC_REML,
  convergence_status = full$convergence_status,
  warnings = full$warnings,
  stringsAsFactors = FALSE
)

clean_h <- clean[c(
  "model", "n_effects", "n_studies", "n_sites", "mean", "ci_lb", "ci_ub",
  "effect_variance", "study_variance", "spatial_variance", "rho_km",
  "logLik_REML", "AIC_REML", "convergence_status", "warnings"
)]
clean_h$dataset <- "published_cleaned"
clean_h <- clean_h[c(
  "dataset", "n_effects", "n_studies", "n_sites", "model", "mean", "ci_lb",
  "ci_ub", "effect_variance", "study_variance", "spatial_variance", "rho_km",
  "logLik_REML", "AIC_REML", "convergence_status", "warnings"
)]

comparison <- rbind(full_h, clean_h)
comparison$model <- factor(comparison$model, levels = expected_models)
comparison <- comparison[order(comparison$dataset, comparison$model), , drop = FALSE]
comparison$model <- as.character(comparison$model)

comparison$delta_AIC_within_dataset <- ave(
  comparison$AIC_REML, comparison$dataset, FUN = function(x) x - min(x)
)

write.csv(comparison,
          file.path(out_dir, "full_vs_published_cleaned_comparison.csv"),
          row.names = FALSE)

wide_full <- full_h[match(expected_models, full_h$model), ]
wide_clean <- clean_h[match(expected_models, clean_h$model), ]
deltas <- data.frame(
  model = expected_models,
  delta_mean_clean_minus_full = wide_clean$mean - wide_full$mean,
  delta_effect_variance_clean_minus_full = wide_clean$effect_variance - wide_full$effect_variance,
  delta_study_variance_clean_minus_full = wide_clean$study_variance - wide_full$study_variance,
  delta_spatial_variance_clean_minus_full = wide_clean$spatial_variance - wide_full$spatial_variance,
  delta_rho_km_clean_minus_full = wide_clean$rho_km - wide_full$rho_km,
  stringsAsFactors = FALSE
)
write.csv(deltas, file.path(out_dir, "cleaned_minus_full_deltas.csv"), row.names = FALSE)

cat("COMPARISON_COMPILED\n")
