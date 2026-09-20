#!/usr/bin/env Rscript

# Regression checks for the public package-specific visualisations. The checks
# read saved interval artifacts and image files only; they never fit or profile
# a statistical model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--mode" ||
    !args[[2L]] %in% c("source", "intervals", "assets", "precompute", "all")) {
  stop("Usage: Rscript validate_visualisation_intervals.R --mode {source|intervals|assets|precompute|all}", call. = FALSE)
}
mode <- args[[2L]]
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
source_path <- file.path(root, "tutorial_v2.qmd")
source_text <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
interval_dir <- file.path(root, "revision_checks", "visualisation_interval_outputs")

assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
required_images <- c(
  "revision_checks/moura_bm_brms_precompute_outputs/moura_bm_brms_parameter_distributions.png",
  "figs/tutorial/tmb_eg1_1.png",
  "figs/tutorial/metafor_lim.png",
  "revision_checks/brms_consistency_outputs/lim_mr_adapt_delta_0.99/posterior_parameter_distributions.png",
  "figs/tutorial/lim_glmmtmb_orchard.png",
  "figs/tutorial/lim_glmmtmb_results.png",
  "figs/tutorial/spain_metafor_spatial_only.png",
  "figs/tutorial/spain_glmmtmb_spatial_only.png",
  "revision_checks/regional_cross_package_audit_outputs/brms_output/brms_parameter_distributions.png",
  "figs/tutorial/spain_cross_package_results.png"
)

check_source <- function() {
  assert(!grepl("posterior-distribution figure is shown", source_text, fixed = TRUE),
         "A public brms visualisation tab still redirects the reader instead of showing a figure.")
  assert(all(vapply(required_images, function(path) grepl(path, source_text, fixed = TRUE), logical(1))),
         "One or more required package-specific visualisation images are not referenced by tutorial_v2.qmd.")
  assert(grepl("Moura BM meta-analysis fitted with `brms`", source_text, fixed = TRUE) &&
           grepl("Lim BM meta-regression fitted with `brms`", source_text, fixed = TRUE),
         "One or more brms visualisation panels are missing their local posterior figure.")
  assert(grepl("Lim BM meta-regression fitted with `glmmTMB`", source_text, fixed = TRUE) &&
           grepl("Spain regional subset: spatial-only `glmmTMB` result", source_text, fixed = TRUE),
         "One or more glmmTMB variance-component figures are missing from the tutorial.")
  assert(grepl("95% profile-likelihood confidence intervals", source_text, fixed = TRUE) &&
           grepl("95% Wald confidence intervals", source_text, fixed = TRUE) &&
           grepl("95% credible intervals", source_text, fixed = TRUE),
         "The tutorial does not distinguish the three displayed interval types.")
  message("VISUALISATION_SOURCE_AUDIT_PASSED")
}

check_assets <- function() {
  paths <- file.path(root, required_images)
  assert(all(file.exists(paths)), "A required visualisation image is missing.")
  assert(all(file.info(paths)$size > 0L), "A required visualisation image is empty.")
  message("VISUALISATION_ASSET_AUDIT_PASSED")
}

read_interval <- function(name) {
  path <- file.path(interval_dir, name)
  assert(file.exists(path), paste("Missing interval artifact:", name))
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")
  assert(all(needed %in% names(dat)), paste("Malformed interval artifact:", name))
  assert(all(is.finite(as.matrix(dat[c("estimate", "ci_lb", "ci_ub")]))) &&
           all(dat$ci_lb <= dat$estimate) && all(dat$estimate <= dat$ci_ub),
         paste("Invalid interval values:", name))
  dat
}

check_precompute <- function() {
  moura <- read_interval("moura_frequentist_intervals.csv")
  lim <- read_interval("lim_frequentist_intervals.csv")
  spain <- read_interval("spain_metafor_intervals.csv")
  assert(all(c("metafor", "glmmTMB") %in% moura$package), "Moura interval output lacks a frequentist package.")
  assert(all(c("metafor", "glmmTMB") %in% lim$package), "Lim interval output lacks a frequentist package.")
  assert(all(c("IID effect-size variance", "Spatial variance", "Exponential range (km)") %in% spain$parameter),
         "Spain metafor output lacks one of its required profile intervals.")
  assert(file.exists(file.path(interval_dir, "moura_metafor_bm.rds")) &&
           file.exists(file.path(interval_dir, "moura_glmmTMB_bm.rds")) &&
           file.exists(file.path(interval_dir, "lim_metafor_bm_meta_regression.rds")) &&
           file.exists(file.path(interval_dir, "lim_glmmTMB_bm_meta_regression.rds")) &&
           file.exists(file.path(interval_dir, "spain_metafor_spatial_only.rds")),
         "One or more Totoro precomputed fit artifacts are absent.")
  assert(!file.exists(file.path(root, "Rdata", "tutorial_v2", "moura_glmmTMB_bm.rds")) &&
           !file.exists(file.path(root, "Rdata", "tutorial_v2", "lim_glmmTMB_bm_meta_regression.rds")),
         "A newly re-estimated frequentist fit was written into the tutorial RDS bundle.")
  message("VISUALISATION_PRECOMPUTE_AUDIT_PASSED")
}

check_intervals <- function() {
  check_precompute()
  spain_gt <- read_interval("spain_glmmTMB_wald_intervals.csv")
  assert(all(c("IID effect-size variance", "Spatial variance", "Exponential range (km)") %in% spain_gt$parameter),
         "Spain glmmTMB Wald interval artifact is incomplete.")
  assert(all(grepl("Wald|profile-likelihood|credible|confidence", c(
    read_interval("moura_frequentist_intervals.csv")$interval_method,
    read_interval("lim_frequentist_intervals.csv")$interval_method,
    read_interval("spain_metafor_intervals.csv")$interval_method,
    spain_gt$interval_method
  ))), "An interval artifact lacks a declared inferential method.")
  message("VISUALISATION_INTERVAL_AUDIT_PASSED")
}

if (mode %in% c("source", "all")) check_source()
if (mode %in% c("assets", "all")) check_assets()
if (mode %in% c("precompute", "all")) check_precompute()
if (mode %in% c("intervals", "all")) check_intervals()
