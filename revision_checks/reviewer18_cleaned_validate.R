args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Specify i2, profiles, comparison, scope, or report")
stage <- args[[1L]]
root <- normalizePath(".")
out <- file.path(root, "revision_checks", "reviewer18_cleaned_primary_outputs")
r18 <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")

if (stage == "i2") {
  prepared <- readRDS(file.path(r18, "published_cleaned_prepared.rds"))
  fit <- readRDS(file.path(r18, "unstructured_only.rds"))
  reported_sampling <- read.csv(file.path(out, "cleaned_sampling_variance_summary.csv"))
  reported_i2 <- read.csv(file.path(out, "cleaned_generalized_i2.csv"))
  vi <- prepared$dat$var_Hedges
  X <- as.matrix(fit$X)
  W <- diag(1 / vi)
  P <- W - W %*% X %*% solve(crossprod(X, W %*% X)) %*% t(X) %*% W
  direct_v_tilde <- (nrow(X) - ncol(X)) / sum(diag(P))
  stopifnot(nrow(reported_sampling) == 1L,
            abs(direct_v_tilde - reported_sampling$generalized_v_tilde) < 1e-12,
            abs(sum(diag(P)) - reported_sampling$trace_P) < 1e-8,
            setequal(unique(reported_i2$model),
                     c("unstructured_only", "spatial_only", "combined")))
  for (model in unique(reported_i2$model)) {
    x <- reported_i2[reported_i2$model == model, ]
    stopifnot(abs(x$I2_percent[x$component == "total"] -
                    sum(x$I2_percent[x$component != "total"])) < 1e-10)
  }
  cat("CLEANED_I2_VALIDATED\n")
} else if (stage == "profiles") {
  grid <- read.csv(file.path(out, "profile_grid.csv"), stringsAsFactors = FALSE)
  points <- read.csv(file.path(out, "cleaned_profile_results_compiled.csv"),
                     stringsAsFactors = FALSE)
  summary <- read.csv(file.path(out, "cleaned_profile_summary.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  multistart <- read.csv(file.path(out, "targeted_multistart_compiled.csv"),
                         stringsAsFactors = FALSE)
  key <- function(x) paste(x$model, x$component, sprintf("%.12g", x$value), sep = "|")
  spatial_free_logLik <- unique(points$free_logLik_REML[points$model == "spatial_only"])
  stopifnot(nrow(grid) == 43L, nrow(points) == 43L,
            !anyDuplicated(key(points)), setequal(key(grid), key(points)),
            all(points$status != "error"), all(is.finite(points$logLik_REML)),
            nrow(summary) == 4L,
            nrow(multistart) == 9L, all(multistart$status != "error"),
            max(multistart$improvement_over_single_start) < 0.10,
            length(spatial_free_logLik) == 1L,
            max(multistart$logLik_REML) < spatial_free_logLik + 0.10,
            all(abs(summary$fixed_primary_minus_free_logLik) < 0.05),
            all(summary$grid_max_minus_free_logLik < 0.10),
            summary$zero_logLik_loss[summary$model == "spatial_only" &
                                       summary$component == "tau2"] > 1.920729,
            summary$zero_logLik_loss[summary$model == "combined" &
                                       summary$component == "tau2"] < 1.920729)
  cat("CLEANED_PROFILES_VALIDATED\n")
} else if (stage == "comparison") {
  i2 <- read.csv(file.path(out, "cleaned_vs_full_generalized_i2.csv"),
                 stringsAsFactors = FALSE)
  ranks <- read.csv(file.path(out, "cleaned_vs_full_aic_ranks.csv"),
                    stringsAsFactors = FALSE)
  models <- read.csv(file.path(r18, "full_vs_published_cleaned_comparison.csv"),
                     stringsAsFactors = FALSE)
  rank_signature <- function(dataset) {
    x <- ranks[ranks$dataset == dataset, ]
    x$model[order(x$AIC_rank)]
  }
  stopifnot(nrow(i2) == 20L,
            setequal(unique(i2$dataset), c("all_spatially_usable", "published_cleaned")),
            nrow(ranks) == 6L,
            identical(rank_signature("all_spatially_usable"),
                      rank_signature("published_cleaned")),
            max(abs(models$mean[models$dataset == "published_cleaned"] -
                      models$mean[models$dataset == "all_spatially_usable"])) < 0.03)
  cat("CLEANED_COMPARISON_VALIDATED\n")
} else if (stage == "scope") {
  qmd_status <- system2("git", c("diff", "--quiet", "--", "tutorial_v2.qmd"))
  gaussian_hits <- list.files(out, pattern = "gau|spgau|gaussian", recursive = TRUE,
                              ignore.case = TRUE, full.names = TRUE)
  forbidden_fit_names <- file.path(out, paste0(c("unstructured_only", "spatial_only",
                                                 "combined"), ".rds"))
  stopifnot(qmd_status == 0L, length(gaussian_hits) == 0L,
            !any(file.exists(forbidden_fit_names)))
  cat("CLEANED_SCOPE_VALIDATED\n")
} else if (stage == "report") {
  report_path <- file.path(root, "revision_checks",
                           "reviewer18_cleaned_primary_2026-09-09.md")
  report <- readLines(report_path, warn = FALSE)
  required <- c("Evidence", "Interpretation", "Uncertainty", "Recommendation",
                "generalized I2", "spatial variance", "rho identifiability",
                "pooled mean", "AIC ranking", "variance allocation",
                "full-data sensitivity", "Gaussian")
  stopifnot(all(vapply(required, function(x) any(grepl(x, report, fixed = TRUE)),
                       logical(1))))
  cat("CLEANED_REPORT_VALIDATED\n")
} else stop("Unknown stage: ", stage)
