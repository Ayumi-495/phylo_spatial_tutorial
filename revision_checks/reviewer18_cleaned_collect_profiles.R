# Compile incrementally saved cleaned-data profile points. No models are fit.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
out_dir <- file.path(project_root, "revision_checks", "reviewer18_cleaned_primary_outputs")
profile_dir <- file.path(out_dir, "profiles")
r18_dir <- file.path(project_root, "revision_checks", "reviewer18_influential_effects_outputs")

grid <- read.csv(file.path(out_dir, "profile_grid.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(grid) == 43L,
          all(grid$model %in% c("spatial_only", "combined")),
          all(grid$component %in% c("tau2", "rho")),
          all(is.finite(grid$value)), all(grid$value >= 0))

profile_files <- list.files(profile_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(profile_files) == nrow(grid))
points <- do.call(rbind, lapply(profile_files, function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(nrow(x) == 1L)
  x
}))
row.names(points) <- NULL

key <- function(x) paste(x$model, x$component, sprintf("%.12g", x$value), sep = "|")
stopifnot(!anyDuplicated(key(grid)), !anyDuplicated(key(points)),
          setequal(key(grid), key(points)), all(points$status != "error"),
          all(is.finite(points$logLik_REML)), all(is.finite(points$AIC_REML)))

primary <- read.csv(file.path(r18_dir, "published_cleaned_primary_results.csv"),
                    stringsAsFactors = FALSE)
primary <- primary[primary$model %in% c("spatial_only", "combined"), ]
stopifnot(nrow(primary) == 2L)

points$free_logLik_REML <- primary$logLik_REML[match(points$model, primary$model)]
points$delta_logLik_from_free <- points$logLik_REML - points$free_logLik_REML
profile_group <- interaction(points$model, points$component)
points$profile_group_max <- ave(points$logLik_REML, profile_group, FUN = max)
points$delta_logLik_from_grid_max <- points$logLik_REML - points$profile_group_max
points <- points[order(points$model, points$component, points$value), ]
write.csv(points, file.path(out_dir, "cleaned_profile_results_compiled.csv"), row.names = FALSE)

primary_parameter <- function(model, component) {
  row <- primary[primary$model == model, ]
  if (component == "tau2") row$spatial_variance else row$rho_km
}

summarise_profile <- function(model, component) {
  x <- points[points$model == model & points$component == component, ]
  free_ll <- unique(x$free_logLik_REML)
  stopifnot(length(free_ll) == 1L)
  primary_value <- primary_parameter(model, component)
  primary_idx <- which.min(abs(x$value - primary_value))
  zero_idx <- if (component == "tau2") which(x$value == 0) else integer()
  near95 <- x$value[x$logLik_REML >= free_ll - 1.920729]
  near_half <- x$value[x$logLik_REML >= free_ll - 0.5]
  data.frame(
    model = model,
    component = component,
    primary_value = primary_value,
    primary_free_logLik_REML = free_ll,
    fixed_primary_grid_value = x$value[primary_idx],
    fixed_primary_logLik_REML = x$logLik_REML[primary_idx],
    fixed_primary_minus_free_logLik = x$logLik_REML[primary_idx] - free_ll,
    grid_max_value = x$value[which.max(x$logLik_REML)],
    grid_max_logLik_REML = max(x$logLik_REML),
    grid_max_minus_free_logLik = max(x$logLik_REML) - free_ll,
    zero_logLik_loss = if (length(zero_idx)) free_ll - x$logLik_REML[zero_idx] else NA_real_,
    near95_min = if (length(near95)) min(near95) else NA_real_,
    near95_max = if (length(near95)) max(near95) else NA_real_,
    near0.5_min = if (length(near_half)) min(near_half) else NA_real_,
    near0.5_max = if (length(near_half)) max(near_half) else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, list(
  summarise_profile("spatial_only", "tau2"),
  summarise_profile("spatial_only", "rho"),
  summarise_profile("combined", "tau2"),
  summarise_profile("combined", "rho")
))
write.csv(summary, file.path(out_dir, "cleaned_profile_summary.csv"), row.names = FALSE)

multistart_dir <- file.path(out_dir, "targeted_multistart")
multistart_files <- list.files(multistart_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(multistart_files) == 9L)
multistart <- do.call(rbind, lapply(multistart_files, read.csv,
                                   stringsAsFactors = FALSE))
row.names(multistart) <- NULL
base_tau <- points[points$model == "spatial_only" & points$component == "tau2",
                   c("value", "logLik_REML")]
multistart$single_start_logLik_REML <- base_tau$logLik_REML[
  match(multistart$fixed_tau2, base_tau$value)]
multistart$improvement_over_single_start <-
  multistart$logLik_REML - multistart$single_start_logLik_REML
multistart <- multistart[order(multistart$fixed_tau2, multistart$rho_start), ]
write.csv(multistart, file.path(out_dir, "targeted_multistart_compiled.csv"),
          row.names = FALSE)

multistart_summary <- do.call(rbind, lapply(split(multistart, multistart$fixed_tau2),
                                            function(x) data.frame(
  fixed_tau2 = x$fixed_tau2[1],
  best_logLik_REML = max(x$logLik_REML),
  worst_logLik_REML = min(x$logLik_REML),
  max_start_dependence_logLik = max(x$logLik_REML) - min(x$logLik_REML),
  max_improvement_over_single_start = max(x$improvement_over_single_start),
  best_final_rho_km = x$final_rho_km[which.max(x$logLik_REML)],
  stringsAsFactors = FALSE
)))
write.csv(multistart_summary, file.path(out_dir, "targeted_multistart_summary.csv"),
          row.names = FALSE)

# Harmonize the saved full- and cleaned-data generalized I2 partitions.
full_i2 <- read.csv(file.path(project_root, "revision_checks", "i2_definition_audit_outputs",
                              "grau_i2.csv"), stringsAsFactors = FALSE)
clean_i2 <- read.csv(file.path(out_dir, "cleaned_generalized_i2.csv"),
                     stringsAsFactors = FALSE)
full_i2 <- data.frame(
  dataset = "all_spatially_usable",
  model = full_i2$model,
  component = full_i2$component,
  variance = full_i2$variance,
  generalized_v_tilde = full_i2$generalized_v_tilde,
  I2_percent = full_i2$I2_percent,
  stringsAsFactors = FALSE
)
clean_i2 <- clean_i2[c("dataset", "model", "component", "variance",
                       "generalized_v_tilde", "I2_percent")]
i2_comparison <- rbind(full_i2, clean_i2)
i2_comparison <- i2_comparison[order(i2_comparison$model, i2_comparison$component,
                                     i2_comparison$dataset), ]
write.csv(i2_comparison, file.path(out_dir, "cleaned_vs_full_generalized_i2.csv"),
          row.names = FALSE)

full_clean_models <- read.csv(file.path(r18_dir, "full_vs_published_cleaned_comparison.csv"),
                              stringsAsFactors = FALSE)
model_order <- c("unstructured_only", "spatial_only", "combined")
rank_table <- do.call(rbind, lapply(split(full_clean_models, full_clean_models$dataset),
                                    function(x) {
  x$AIC_rank <- rank(x$AIC_REML, ties.method = "min")
  x[c("dataset", "model", "AIC_rank", "delta_AIC_within_dataset")]
}))
rank_table <- rank_table[order(match(rank_table$model, model_order), rank_table$dataset), ]
write.csv(rank_table, file.path(out_dir, "cleaned_vs_full_aic_ranks.csv"), row.names = FALSE)

cat("CLEANED_PROFILES_COMPILED\n")
print(summary, digits = 12, row.names = FALSE)
