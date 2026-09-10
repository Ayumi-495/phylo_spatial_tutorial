# Final verification gates for the spatial numerical audit.
# Usage: Rscript revision_checks/validate_spatial_audit_final.R <stage>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one validation stage.")
stage <- args[[1L]]
allowed <- c("grau_gaussian", "scholer_primary", "scholer_profiles",
             "scholer_i2", "grau_gaussian_i2", "synthesis", "protected_sources")
stopifnot(stage %in% allowed)

root <- normalizePath(".")
if (basename(root) == "revision_checks") {
  root <- dirname(root)
  setwd(root)
}
checks <- file.path(root, "revision_checks")
gauss <- file.path(checks, "gaussian_global_outputs")
scholer <- file.path(checks, "scholer_spatial_audit_outputs")

close_enough <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol))
}

read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

if (stage == "grau_gaussian") {
  target <- file.path(gauss, "targeted_free_refits")
  compiled_path <- file.path(target, "targeted_free_refits_compiled.csv")
  report_path <- file.path(checks, "grau_global_gaussian_kernel_audit_2026-09-09.md")
  stopifnot(file.exists(compiled_path), file.exists(report_path),
            file.exists(file.path(gauss, "combined_spgau_profile_tau2.csv")),
            file.exists(file.path(gauss, "combined_spgau_profile_rho.csv")),
            file.exists(file.path(gauss, "fixed_profiles", "combined_spgau_tau2_0.csv")))
  tab <- read.csv(compiled_path, stringsAsFactors = FALSE)
  stopifnot(nrow(tab) == 3L,
            setequal(tab$start_label, c("fixed200", "intermediate800", "primary3091")),
            all(is.finite(tab$logLik_REML)), all(is.finite(tab$AIC_REML)),
            all(is.finite(tab$spatial_variance)), all(tab$spatial_variance >= 0),
            all(is.finite(tab$rho_km)), all(tab$rho_km > 0),
            !any(grepl("error|not_converged", tab$convergence_status, ignore.case = TRUE)),
            all(abs(tab$AIC_REML - (-2 * tab$logLik_REML + 10)) < 1e-5))
  for (label in tab$start_label) {
    fit <- readRDS(file.path(target, paste0("free_", label, ".rds")))
    row <- tab[tab$start_label == label, ]
    stopifnot(close_enough(fit$fit.stats["ll", "REML"], row$logLik_REML),
              close_enough(fit$sigma2[1], row$iid_effect_variance),
              close_enough(fit$sigma2[2], row$study_variance),
              close_enough(fit$tau2[1], row$spatial_variance),
              close_enough(fit$rho[1], row$rho_km))
  }
  zero <- read.csv(file.path(gauss, "fixed_profiles", "combined_spgau_tau2_0.csv"),
                   stringsAsFactors = FALSE)
  stopifnot(max(tab$logLik_REML) - zero$logLik_REML[1] < 1.92)
  report <- read_text(report_path)
  stopifnot(grepl("targeted", report, ignore.case = TRUE),
            grepl("weakly|flat|optimizer", report, ignore.case = TRUE),
            grepl("tau2|spatial variance", report, ignore.case = TRUE))
  cat("GRAU_GAUSSIAN_VALIDATED\n")
}

if (stage == "scholer_primary") {
  prepared <- readRDS(file.path(checks, "scholer_structure_audit_outputs", "scholer_prepared.rds"))
  dat <- prepared$dat
  D <- prepared$distance_km
  stopifnot(nrow(dat) == 949L, nlevels(dat$study_id) == 205L,
            nlevels(dat$site_id) == 454L,
            identical(rownames(D), levels(dat$site_id)),
            identical(colnames(D), levels(dat$site_id)))
  expected <- list(
    unstructured_only = list(s = c("effect_id", "study_id"), g = character(), struct = NULL),
    spatial_only = list(s = "effect_id", g = c("site_id", "const"), struct = "SPEXP"),
    combined = list(s = c("effect_id", "study_id"), g = c("site_id", "const"), struct = "SPEXP")
  )
  fits <- list()
  for (nm in names(expected)) {
    fit <- readRDS(file.path(scholer, paste0(nm, ".rds")))
    fits[[nm]] <- fit
    stopifnot(inherits(fit, "rma.mv"), fit$k == 949L, fit$p == 1L,
              identical(unname(fit$s.names), expected[[nm]]$s),
              identical(as.character(unname(unlist(fit$g.names))), expected[[nm]]$g),
              isTRUE(all.equal(fit$vi, dat$vi)), all(fit$X[, 1] == 1))
    if (!is.null(expected[[nm]]$struct)) stopifnot(all(fit$struct == expected[[nm]]$struct))
  }
  tab <- read.csv(file.path(scholer, "scholer_primary_model_results.csv"),
                  stringsAsFactors = FALSE)
  stopifnot(nrow(tab) == 3L, setequal(tab$model, names(expected)),
            all(is.finite(tab$pooled_mean)), all(is.finite(tab$AIC_REML)),
            !any(grepl("error|not_converged", tab$convergence_status, ignore.case = TRUE)))
  cat("SCHOLER_PRIMARY_VALIDATED\n")
}

if (stage == "scholer_profiles") {
  profile_dir <- file.path(scholer, "profiles")
  csvs <- list.files(profile_dir, pattern = "\\.csv$", full.names = TRUE)
  rdss <- list.files(profile_dir, pattern = "\\.rds$", full.names = TRUE)
  stopifnot(length(csvs) == 39L, length(rdss) == 39L)
  tab <- read.csv(file.path(scholer, "scholer_profile_points_compiled.csv"),
                  stringsAsFactors = FALSE)
  summary <- read.csv(file.path(scholer, "scholer_profile_summary.csv"),
                      stringsAsFactors = FALSE)
  stopifnot(nrow(tab) == 39L, nrow(summary) == 4L,
            !any(grepl("error|not_converged", tab$convergence_status, ignore.case = TRUE)))
  c_tau <- summary[summary$model == "combined" & summary$parameter == "tau2", ]
  s_tau <- summary[summary$model == "spatial_only" & summary$parameter == "tau2", ]
  c_rho <- summary[summary$model == "combined" & summary$parameter == "rho", ]
  s_rho <- summary[summary$model == "spatial_only" & summary$parameter == "rho", ]
  stopifnot(nrow(c_tau) == 1L, nrow(s_tau) == 1L, nrow(c_rho) == 1L, nrow(s_rho) == 1L,
            c_tau$tau2_zero_delta_logLik < 1.92,
            s_tau$tau2_zero_delta_logLik > 100,
            close_enough(c_rho$primary_estimate, c_rho$grid_max_value, tol = 1e-6),
            close_enough(s_rho$primary_estimate, s_rho$grid_max_value, tol = 1e-6),
            grepl("10 to 12000", c_rho$grid_values_with_delta_le_1_92, fixed = TRUE))
  cat("SCHOLER_PROFILES_VALIDATED\n")
}

if (stage == "scholer_i2") {
  sampling <- read.csv(file.path(scholer, "scholer_i2_sampling_variance.csv"))
  tab <- read.csv(file.path(scholer, "scholer_i2_results.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(sampling) == 1L, nrow(tab) == 10L,
            close_enough(sampling$generalized_v_tilde, 0.000725865052242, tol = 1e-11),
            abs(sampling$ratio_minus_generalized) < 1e-12,
            !close_enough(sampling$generalized_v_tilde, sampling$arithmetic_mean_vi, tol = 1e-3))
  for (nm in unique(tab$model)) {
    z <- tab[tab$model == nm, ]
    stopifnot(abs(z$I2_percent[z$component == "total"] -
                    sum(z$I2_percent[z$component != "total"])) < 1e-8)
  }
  cat("SCHOLER_I2_VALIDATED\n")
}

if (stage == "grau_gaussian_i2") {
  sampling <- read.csv(file.path(gauss, "grau_spgau_i2_sampling_variance.csv"))
  tab <- read.csv(file.path(gauss, "grau_spgau_i2_results.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(sampling) == 1L,
            close_enough(sampling$generalized_v_tilde, 0.111123796928, tol = 1e-11),
            abs(sampling$difference_from_verified_target) < 1e-10,
            nrow(tab) == 15L)
  keys <- unique(tab[c("model", "solution")])
  for (i in seq_len(nrow(keys))) {
    z <- tab[tab$model == keys$model[i] & tab$solution == keys$solution[i], ]
    stopifnot(abs(z$I2_percent[z$component == "total"] -
                    sum(z$I2_percent[z$component != "total"])) < 1e-8)
  }
  cat("GRAU_GAUSSIAN_I2_VALIDATED\n")
}

if (stage == "synthesis") {
  report_path <- file.path(checks, "spatial_audit_final_synthesis_2026-09-09.md")
  table_path <- file.path(checks, "spatial_audit_final_comparison.csv")
  stopifnot(file.exists(report_path), file.exists(table_path))
  tab <- read.csv(table_path, stringsAsFactors = FALSE)
  stopifnot(all(c("Grau-Andres", "Scholer") %in% unique(tab$dataset)),
            all(c("SPEXP", "SPGAU") %in% unique(tab$kernel)),
            all(c("unstructured_only", "spatial_only", "combined") %in% unique(tab$model)),
            all(is.finite(tab$pooled_mean)), all(is.finite(tab$AIC_REML)))
  text <- read_text(report_path)
  stopifnot(grepl("pooled mean", text, ignore.case = TRUE),
            grepl("weakly identified", text, ignore.case = TRUE),
            grepl("small estimated spatial variance", text, ignore.case = TRUE),
            grepl("I²", text, fixed = TRUE))
  cat("SPATIAL_SYNTHESIS_VALIDATED\n")
}

if (stage == "protected_sources") {
  changed <- system2("git", c("diff", "--name-only", "0e57293", "--"), stdout = TRUE)
  status <- system2("git", c("status", "--porcelain=v1", "--untracked-files=all"), stdout = TRUE)
  status_paths <- if (length(status)) substring(status, 4L) else character()
  all_changed <- unique(c(changed, status_paths))
  forbidden <- function(paths) grepl("(^|/)(tutorial_v2\\.qmd|.*response.*\\.(md|docx|txt)|.*manuscript.*\\.(md|qmd|docx|tex|pdf))$",
                                     paths, ignore.case = TRUE)
  # Positive control proves the absence checker recognizes a protected source.
  stopifnot(forbidden("tutorial_v2.qmd"), !any(forbidden(all_changed)))
  cat("PROTECTED_SOURCES_UNCHANGED\n")
}
