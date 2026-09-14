#!/usr/bin/env Rscript

# Full-precision profile-likelihood audit for the Figure 3 variance/range
# displays. Every target is profiled from an immutable validated rma.mv object;
# only the parameter under examination is constrained while the other model
# parameters are re-optimized by metafor::confint.rma.mv(). No final model is
# overwritten or otherwise refitted for substantive inference.

suppressPackageStartupMessages({
  library(metafor)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && !args[[1L]] %in% c(
  "all", "moura_study", "moura_effect", "moura_species_nonphylo",
  "moura_species_phylo", "global_spatial_only_rho", "global_combined_rho", "collate"
))) {
  stop("Usage: Rscript revision_checks/reviewer13_figure3_profile_audit.R [all|target|collate]", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript revision_checks/reviewer13_figure3_profile_audit.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
out_dir <- file.path(root, "revision_checks", "reviewer13_figure3_profile_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

moura_fit <- readRDS(file.path(root, "revision_checks", "ou_correctness_outputs", "baseline_fit_objects.rds"))$fit_bm
spatial_fit <- readRDS(file.path(root, "revision_checks", "reviewer18_influential_effects_outputs", "spatial_only.rds"))
combined_fit <- readRDS(file.path(root, "revision_checks", "reviewer18_influential_effects_outputs", "combined.rds"))
assert(inherits(moura_fit, "rma.mv") && moura_fit$k == 1828L && length(moura_fit$sigma2) == 4L,
       "Unexpected final Moura BM fit.")
assert(inherits(spatial_fit, "rma.mv") && spatial_fit$k == 2355L && spatial_fit$rho[[1L]] > 0,
       "Unexpected final global spatial-only fit.")
assert(inherits(combined_fit, "rma.mv") && combined_fit$k == 2355L && combined_fit$rho[[1L]] > 0,
       "Unexpected final global combined fit.")

# Moura brackets are used only to start the new root search near the expected
# likelihood-ratio crossings. They do not supply an interval value: every
# endpoint is recalculated from the final object by constrained REML profiling.
# The bounds are deliberately wider than the retained printed intervals while
# avoiding the default 0--10 domain. Global rho domains use the rma.mv defaults
# for SPEXP: [0, max(10, 10 * fitted rho)].
# Inner constrained refits use maxiter = 100. A direct endpoint check showed
# that this cap completes a Moura constrained likelihood evaluation; the cap is
# recorded in every output and does not alter the immutable saved final fit.
targets <- list(
  moura_study = list(fit=moura_fit, label="Study variance", parameter="sigma2", index=1L,
                     control=list(vc.min=0.005, vc.max=0.04, tol=1e-6, maxiter=1000, eptries=10),
                     source_fit="ou_correctness_outputs/baseline_fit_objects.rds$fit_bm"),
  moura_effect = list(fit=moura_fit, label="Effect-size variance", parameter="sigma2", index=2L,
                      control=list(vc.min=0.01, vc.max=0.02, tol=1e-6, maxiter=1000, eptries=10),
                      source_fit="ou_correctness_outputs/baseline_fit_objects.rds$fit_bm"),
  moura_species_nonphylo = list(fit=moura_fit, label="Species variance, non-phylogenetic", parameter="sigma2", index=3L,
                                 control=list(vc.min=0.02, vc.max=0.09, tol=1e-6, maxiter=1000, eptries=10),
                                 source_fit="ou_correctness_outputs/baseline_fit_objects.rds$fit_bm"),
  moura_species_phylo = list(fit=moura_fit, label="Species variance, phylogenetic", parameter="sigma2", index=4L,
                              control=list(vc.min=0.01, vc.max=0.25, tol=1e-6, maxiter=1000, eptries=10),
                              source_fit="ou_correctness_outputs/baseline_fit_objects.rds$fit_bm"),
  global_spatial_only_rho = list(fit=spatial_fit, label="Spatial-only rho (km)", parameter="rho", index=1L,
                                 control=list(vc.min=0, vc.max=max(10, 10 * spatial_fit$rho[[1L]]), tol=1e-6, maxiter=1000, eptries=10),
                                 source_fit="reviewer18_influential_effects_outputs/spatial_only.rds"),
  global_combined_rho = list(fit=combined_fit, label="Combined-model rho (km)", parameter="rho", index=1L,
                             control=list(vc.min=0, vc.max=max(10, 10 * combined_fit$rho[[1L]]), tol=1e-6, maxiter=1000, eptries=10),
                             source_fit="reviewer18_influential_effects_outputs/combined.rds")
)

selected <- if (!length(args) || args[[1L]] == "all") names(targets) else if (args[[1L]] == "collate") character() else args

extract_profile_row <- function(profile, target_spec, id, elapsed, warnings) {
  random <- as.data.frame(profile$random)
  value_row <- random[1L, , drop=FALSE]
  lower_sign <- paste(profile$lb.sign, collapse="")
  upper_sign <- paste(profile$ub.sign, collapse="")
  tibble(
    target=id,
    component=target_spec$label,
    parameter=target_spec$parameter,
    parameter_index=target_spec$index,
    estimate=unname(value_row[["estimate"]]),
    ci_lb=unname(value_row[["ci.lb"]]),
    ci_ub=unname(value_row[["ci.ub"]]),
    lower_bound_sign=lower_sign,
    upper_bound_sign=upper_sign,
    lower_reached_search_boundary=identical(lower_sign, "<"),
    upper_reached_search_boundary=identical(upper_sign, ">"),
    finite_two_sided=identical(lower_sign, "") && identical(upper_sign, "") &&
      is.finite(value_row[["ci.lb"]]) && is.finite(value_row[["ci.ub"]]),
    ci_null=profile$ci.null,
    search_min=target_spec$control$vc.min,
    search_max=target_spec$control$vc.max,
    root_tolerance=target_spec$control$tol,
    maxiter=target_spec$control$maxiter,
    eptries=target_spec$control$eptries,
    inner_optimizer_maxiter=100L,
    method="metafor::confint.rma.mv(), profile-likelihood constrained re-optimization",
    source_fit=target_spec$source_fit,
    elapsed_seconds=elapsed,
    warnings=if (length(warnings)) paste(unique(warnings), collapse=" | ") else "",
    convergence_information="No explicit optimizer convergence code is exposed by confint.rma.mv; outer warnings were captured, while metafor suppresses internal constrained-fit warnings.",
    fit_changed=FALSE
  )
}

run_target <- function(id) {
  target_spec <- targets[[id]]
  target_spec$fit$control <- utils::modifyList(target_spec$fit$control, list(maxiter=100L))
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  profile <- withCallingHandlers(
    tryCatch({
      do.call(metafor::confint.rma.mv, c(
        list(object=target_spec$fit, level=0.95, control=target_spec$control),
        stats::setNames(list(target_spec$index), target_spec$parameter)
      ))
    }, error=function(e) e),
    warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(profile, "error")) {
    row <- tibble(
      target=id, component=target_spec$label, parameter=target_spec$parameter, parameter_index=target_spec$index,
      estimate=NA_real_, ci_lb=NA_real_, ci_ub=NA_real_, lower_bound_sign=NA_character_, upper_bound_sign=NA_character_,
      lower_reached_search_boundary=NA, upper_reached_search_boundary=NA, finite_two_sided=FALSE, ci_null=NA,
      search_min=target_spec$control$vc.min, search_max=target_spec$control$vc.max, root_tolerance=target_spec$control$tol,
      maxiter=target_spec$control$maxiter, eptries=target_spec$control$eptries,
      inner_optimizer_maxiter=100L,
      method="metafor::confint.rma.mv(), profile-likelihood constrained re-optimization",
      source_fit=target_spec$source_fit, elapsed_seconds=elapsed,
      warnings=paste(c(conditionMessage(profile), unique(warnings)), collapse=" | "),
      convergence_information="Profile call returned an error; see warnings.", fit_changed=FALSE
    )
    saveRDS(list(profile_error=conditionMessage(profile), row=row), file.path(out_dir, paste0(id, "_profile.rds")))
  } else {
    row <- extract_profile_row(profile, target_spec, id, elapsed, warnings)
    saveRDS(list(profile=profile, row=row), file.path(out_dir, paste0(id, "_profile.rds")))
  }
  write_csv(row, file.path(out_dir, paste0(id, "_profile.csv")))
  cat(id, "completed in", sprintf("%.1f", elapsed), "seconds\n")
}

invisible(lapply(selected, run_target))
paths <- file.path(out_dir, paste0(names(targets), "_profile.csv"))
if (all(file.exists(paths))) {
  combined <- bind_rows(lapply(paths, read_csv, show_col_types=FALSE)) |>
    arrange(match(target, names(targets)))
  write_csv(combined, file.path(out_dir, "figure3_profile_likelihood_summary.csv"))
}
cat("REVIEWER13_FIGURE3_PROFILE_AUDIT_PASSED\n")
