# Compile incrementally saved Scholer spatial profile points into readout tables.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scholer_profile_summary.R <audit_output_dir>")
out_dir <- args[[1L]]
primary <- read.csv(file.path(out_dir, "scholer_primary_model_results.csv"), stringsAsFactors = FALSE)
files <- list.files(file.path(out_dir, "profiles"), pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No profile-point CSV files found.")
points <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
points <- points[order(points$model, points$parameter, points$fixed_value), , drop = FALSE]
points$delta_logLik_from_grid_max <- NA_real_

summary_rows <- list()
i <- 0L
for (model in unique(points$model)) {
  for (parameter in unique(points$parameter[points$model == model])) {
    ix <- points$model == model & points$parameter == parameter
    subset <- points[ix, , drop = FALSE]
    max_ll <- max(subset$logLik_REML, na.rm = TRUE)
    points$delta_logLik_from_grid_max[ix] <- max_ll - points$logLik_REML[ix]
    best <- subset[which.max(subset$logLik_REML), , drop = FALSE]
    primary_row <- primary[primary$model == model, , drop = FALSE]
    primary_value <- if (parameter == "tau2") primary_row$spatial_variance else primary_row$rho_km
    zero_ll <- if (parameter == "tau2") subset$logLik_REML[subset$fixed_value == 0] else NA_real_
    i <- i + 1L
    summary_rows[[i]] <- data.frame(
      model = model, parameter = parameter, n_profile_points = nrow(subset),
      primary_estimate = primary_value, primary_logLik_REML = primary_row$logLik_REML,
      grid_max_value = best$fixed_value, grid_max_logLik_REML = best$logLik_REML,
      primary_minus_grid_max = max_ll - primary_row$logLik_REML,
      tau2_zero_logLik_REML = zero_ll,
      tau2_zero_delta_logLik = if (parameter == "tau2") max_ll - zero_ll else NA_real_,
      grid_values_with_delta_le_1_92 = paste(range(subset$fixed_value[points$delta_logLik_from_grid_max[ix] <= 1.92]), collapse = " to "),
      stringsAsFactors = FALSE
    )
  }
}
summary <- do.call(rbind, summary_rows)
write.csv(points, file.path(out_dir, "scholer_profile_points_compiled.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "scholer_profile_summary.csv"), row.names = FALSE)

readout <- c(
  "Profile points were fit with the relevant spatial parameter fixed and all other free parameters re-optimized.",
  "Each point was saved independently before compilation. Delta log-likelihood is relative to the maximum among the stated bounded grid points.",
  "A grid range with delta <= 1.92 is descriptive only, not a formal confidence interval."
)
for (j in seq_len(nrow(summary))) {
  row <- summary[j, ]
  readout <- c(readout, paste(
    row$model, row$parameter,
    "primary=", signif(row$primary_estimate, 7),
    "grid maximum=", signif(row$grid_max_value, 7),
    "primary-minus-grid-max delta=", signif(row$primary_minus_grid_max, 6),
    "grid delta<=1.92 range=", row$grid_values_with_delta_le_1_92,
    if (row$parameter == "tau2") paste("tau2=0 delta=", signif(row$tau2_zero_delta_logLik, 6)) else ""
  ))
}
writeLines(readout, file.path(out_dir, "scholer_profile_readout.txt"))
