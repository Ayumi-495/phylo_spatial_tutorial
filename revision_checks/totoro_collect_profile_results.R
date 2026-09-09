# Combine incrementally saved profile-point records without refitting any model.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript totoro_collect_profile_results.R <profile_dir> <output_csv>")
}
profile_dir <- args[[1]]
output_csv <- args[[2]]

files <- list.files(profile_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No profile-point CSV files found.")
rows <- lapply(files, read.csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("model", "component", "value", "logLik_REML", "AIC_REML", "status", "elapsed_seconds", "warnings")
if (!all(vapply(rows, function(x) identical(names(x), required), logical(1)))) {
  stop("Profile-point CSV columns do not match the expected schema.")
}
profiles <- do.call(rbind, rows)
profiles <- profiles[order(profiles$model, profiles$component, profiles$value), ]
profiles$delta_logLik_from_profile_max <- NA_real_
for (i in interaction(profiles$model, profiles$component, drop = TRUE)) {
  index <- which(interaction(profiles$model, profiles$component, drop = TRUE) == i)
  valid <- is.finite(profiles$logLik_REML[index])
  if (any(valid)) {
    target <- index[valid]
    profiles$delta_logLik_from_profile_max[target] <-
      profiles$logLik_REML[target] - max(profiles$logLik_REML[target])
  }
}
write.csv(profiles, output_csv, row.names = FALSE)
