args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Specify provenance, datasets, fits, report, or source")
stage <- args[[1L]]
root <- normalizePath(".")
out <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")

if (stage == "provenance") {
  h <- read.csv(file.path(out, "dryad_file_hashes.csv"), stringsAsFactors = FALSE)
  e <- read.csv(file.path(out, "published_influential_effects.csv"), stringsAsFactors = FALSE)
  archived_code <- readLines(file.path(root, "data", "Roger_etal_2024",
                                       "global_fire_analysis_240325.R"), warn = FALSE)
  archived_text <- paste(archived_code, collapse = "\n")
  code_study_ids <- sub("-[0-9]+$", "", e$effect_uid)
  stopifnot(nrow(h) == 2L, all(h$exact_match), nrow(e) == 7L,
            sum(e$present_in_dryad_csv) == 6L,
            identical(e$effect_uid[!e$present_in_dryad_csv], "Launonen_1999-1"),
            all(e$response[e$present_in_dryad_csv] == e$response_expected[e$present_in_dryad_csv]),
            all(vapply(code_study_ids, grepl, logical(1), x = archived_text, fixed = TRUE)),
            grepl("cooksd.ab>0.015", archived_text, fixed = TRUE),
            grepl("cooksd.di>0.015", archived_text, fixed = TRUE),
            grepl("cooksd.fi>0.06", archived_text, fixed = TRUE),
            grepl("data = reg.fi", archived_text, fixed = TRUE),
            grepl("reg.fit <- reg[reg$response==\"fitness\",]", archived_text, fixed = TRUE))
  cat("PROVENANCE_VALIDATED\n")
} else if (stage == "datasets") {
  c <- read.csv(file.path(out, "dataset_counts.csv"), stringsAsFactors = FALSE)
  a <- read.csv(file.path(out, "all_spatially_usable.csv"), stringsAsFactors = FALSE)
  p <- read.csv(file.path(out, "published_cleaned.csv"), stringsAsFactors = FALSE)
  e <- read.csv(file.path(out, "published_influential_effects.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(a) == 2361L, nrow(p) == 2355L,
            nrow(a) - nrow(p) == sum(e$present_in_dryad_csv),
            nrow(c) == 2L, !any(p$effect_uid %in% e$effect_uid[e$present_in_dryad_csv]))
  cat("DATASETS_VALIDATED\n")
} else if (stage == "fits") {
  r <- read.csv(file.path(out, "published_cleaned_primary_results.csv"), stringsAsFactors = FALSE)
  stopifnot(setequal(r$model, c("unstructured_only", "spatial_only", "combined")),
            nrow(r) == 3L, all(r$n_effects == 2355L),
            all(is.finite(r$mean)), all(is.finite(r$ci_lb)), all(is.finite(r$ci_ub)),
            all(is.finite(r$effect_variance)), all(is.finite(r$logLik_REML)),
            all(is.finite(r$AIC_REML)), all(r$convergence_status != "error"),
            file.exists(file.path(out, "unstructured_only.rds")),
            file.exists(file.path(out, "spatial_only.rds")),
            file.exists(file.path(out, "combined.rds")))
  p <- readRDS(file.path(out, "published_cleaned_prepared.rds"))
  stopifnot(identical(rownames(p$distance_km), levels(p$dat$site_id)),
            identical(colnames(p$distance_km), levels(p$dat$site_id)))
  cat("FITS_VALIDATED\n")
} else if (stage == "report") {
  comparison <- read.csv(file.path(out, "full_vs_published_cleaned_comparison.csv"), stringsAsFactors = FALSE)
  report <- readLines(file.path(root, "revision_checks", "reviewer18_influential_effects_2026-09-09.md"), warn = FALSE)
  required <- c("pooled biological conclusion", "model ranking", "variance allocation",
                "spatial component", "Spain")
  stopifnot(nrow(comparison) == 6L,
            all(vapply(required, function(x) any(grepl(x, report, fixed = TRUE)), logical(1))))
  cat("REPORT_VALIDATED\n")
} else if (stage == "source") {
  status <- system2("git", c("diff", "--quiet", "--", "tutorial_v2.qmd"))
  stopifnot(status == 0L)
  cat("SOURCE_UNCHANGED\n")
} else stop("Unknown stage: ", stage)
