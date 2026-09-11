#!/usr/bin/env Rscript

# Summarise the completed fixed-rho likelihood profile without refitting models.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "revision_checks/ou_correctness_outputs"
profile_path <- file.path(out_dir, "rho_profile.csv")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
assert(file.exists(profile_path), paste("Missing profile:", profile_path))

profile <- read.csv(profile_path, check.names = FALSE)
profile <- profile[order(profile$rho), , drop = FALSE]
assert(all(is.finite(profile$rho)) && all(is.finite(profile$REML_logLik)), "Invalid profile values.")
peak_index <- which.max(profile$REML_logLik)
cutoff <- qchisq(0.95, 1) / 2
supported <- profile$delta_REML_logLik_from_grid_max <= cutoff
assert(any(supported), "No profile points meet the likelihood-ratio criterion.")
lower_fail <- max(which(!supported & seq_len(nrow(profile)) < peak_index))
upper_fail <- min(which(!supported & seq_len(nrow(profile)) > peak_index))
assert(is.finite(lower_fail) && is.finite(upper_fail), "Profile interval is not bracketed by the grid.")
interpolate_crossing <- function(i, j) {
  exp(approx(profile$delta_REML_logLik_from_grid_max[c(i, j)],
             log(profile$rho[c(i, j)]), xout = cutoff)$y)
}
interval <- data.frame(
  likelihood_ratio_cutoff = cutoff,
  peak_rho = profile$rho[peak_index],
  peak_alpha = profile$alpha[peak_index],
  lower_rho_grid_supported = min(profile$rho[supported]),
  upper_rho_grid_supported = max(profile$rho[supported]),
  lower_rho_log_interpolated = interpolate_crossing(lower_fail, lower_fail + 1L),
  upper_rho_log_interpolated = interpolate_crossing(upper_fail - 1L, upper_fail),
  optimum_interior = TRUE,
  stringsAsFactors = FALSE)
interval$lower_alpha_log_interpolated <- 1 / interval$upper_rho_log_interpolated
interval$upper_alpha_log_interpolated <- 1 / interval$lower_rho_log_interpolated

# Retain all likelihood-supported points for a compact pooled-mean/SE sensitivity table.
sensitivity <- profile[supported, c("rho", "alpha", "REML_logLik", "pooled_mean", "pooled_se",
  "ci_lb", "ci_ub", "study_variance", "effect_size_variance",
  "species_nonphylogenetic_variance", "species_phylogenetic_variance",
  "correlation_offdiag_mean", "covariance_offdiag_mean",
  "delta_REML_logLik_from_grid_max")]
write.csv(interval, file.path(out_dir, "rho_profile_likelihood_interval.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(out_dir, "rho_profile_supported_sensitivity.csv"), row.names = FALSE)
cat("OU_PROFILE_SUMMARY_COMPLETED\n")
