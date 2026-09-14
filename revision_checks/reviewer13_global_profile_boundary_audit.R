#!/usr/bin/env Rscript

# Recover the boundary status of every saved global profile shown in the
# supplementary figure. This script reads existing confint.rma objects only;
# it does not refit or overwrite any final model.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript revision_checks/reviewer13_global_profile_boundary_audit.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
legacy_dir <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs", "variance_profile_ci")
profile_dir <- file.path(root, "revision_checks", "reviewer13_figure3_profile_outputs")
out_path <- file.path(profile_dir, "global_profile_boundary_audit.csv")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

variance_targets <- tribble(
  ~model, ~parameter, ~parameter_index, ~component,
  "unstructured_only", "sigma2", 1L, "Effect-size variance",
  "unstructured_only", "sigma2", 2L, "Study variance",
  "spatial_only", "sigma2", 1L, "Effect-size variance",
  "spatial_only", "tau2", 1L, "Spatial variance",
  "combined", "sigma2", 1L, "Effect-size variance",
  "combined", "sigma2", 2L, "Study variance",
  "combined", "tau2", 1L, "Spatial variance"
)

variance_rows <- lapply(seq_len(nrow(variance_targets)), function(i) {
  target <- variance_targets[i, ]
  stem <- sprintf("%s_%s_%d", target$model, target$parameter, target$parameter_index)
  raw_path <- file.path(legacy_dir, paste0(stem, ".rds"))
  csv_path <- file.path(legacy_dir, paste0(stem, ".csv"))
  assert(file.exists(raw_path) && file.exists(csv_path), paste("Missing saved profile", stem))
  ci <- readRDS(raw_path)
  assert(inherits(ci, "confint.rma"), paste("Unexpected profile object", stem))
  random <- as.data.frame(ci$random)[1L, , drop = FALSE]
  saved <- read_csv(csv_path, show_col_types = FALSE)
  estimate <- unname(random[["estimate"]])
  lower <- unname(random[["ci.lb"]])
  upper <- unname(random[["ci.ub"]])
  assert(nrow(saved) == 1L && isTRUE(all.equal(saved$estimate, estimate, tolerance = 1e-12)) &&
           isTRUE(all.equal(saved$ci_lb, lower, tolerance = 1e-12)) &&
           isTRUE(all.equal(saved$ci_ub, upper, tolerance = 1e-12)),
         paste("CSV/RDS mismatch", stem))
  lower_sign <- ci$lb.sign
  upper_sign <- ci$ub.sign
  search_upper <- max(ifelse(estimate <= .Machine$double.eps^0.5, 10, max(10, estimate * 100)), 0)
  tibble(
    model = target$model, component = target$component, parameter = target$parameter,
    parameter_index = target$parameter_index, estimate = estimate,
    profile_lower = lower, profile_upper = upper,
    lower_sign = lower_sign, upper_sign = upper_sign,
    lower_cutoff_reached = identical(lower_sign, ""),
    upper_cutoff_reached = identical(upper_sign, ""),
    lower_is_search_boundary = !identical(lower_sign, ""),
    upper_is_search_boundary = !identical(upper_sign, ""),
    finite_two_sided = identical(lower_sign, "") && identical(upper_sign, ""),
    search_lower = 0, search_upper = search_upper,
    source_profile = sub(paste0(root, "/"), "", raw_path)
  )
}) |> bind_rows()

rho_rows <- read_csv(file.path(profile_dir, "figure3_profile_likelihood_summary.csv"), show_col_types = FALSE) |>
  filter(target %in% c("global_spatial_only_rho", "global_combined_rho")) |>
  transmute(
    model = if_else(target == "global_spatial_only_rho", "spatial_only", "combined"),
    component = "Spatial range rho", parameter = "rho", parameter_index = 1L,
    estimate, profile_lower = ci_lb, profile_upper = ci_ub,
    lower_sign = coalesce(lower_bound_sign, ""), upper_sign = coalesce(upper_bound_sign, ""),
    lower_cutoff_reached = !lower_reached_search_boundary,
    upper_cutoff_reached = !upper_reached_search_boundary,
    lower_is_search_boundary = lower_reached_search_boundary,
    upper_is_search_boundary = upper_reached_search_boundary,
    finite_two_sided, search_lower = search_min, search_upper = search_max,
    source_profile = paste0("reviewer13_figure3_profile_outputs/", target, "_profile.rds")
  )

audit <- bind_rows(variance_rows, rho_rows) |>
  mutate(
    profile_status = case_when(
      finite_two_sided ~ "finite two-sided profile interval",
      lower_is_search_boundary & upper_is_search_boundary ~ "both sides unresolved over searched range",
      lower_is_search_boundary ~ "lower side unresolved over searched range",
      upper_is_search_boundary ~ "upper side unresolved over searched range",
      TRUE ~ "inspect profile status"
    )
  ) |>
  arrange(match(model, c("unstructured_only", "spatial_only", "combined")),
          match(parameter, c("sigma2", "tau2", "rho")), parameter_index)

assert(nrow(audit) == 9L, "Expected seven variance and two rho profiles.")
assert(sum(audit$finite_two_sided) == 6L,
       "Expected six finite two-sided profiles: five variance profiles plus spatial-only tau2.")
write_csv(audit, out_path)
cat("GLOBAL_PROFILE_BOUNDARY_AUDIT_PASSED\n")
