# Validate the revised spatial tutorial against finalized audit CSV outputs.
# This script does not refit any model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || args[[1L]] != "spatial_tutorial") {
  stop("Usage: Rscript revision_checks/validate_spatial_tutorial_revision.R spatial_tutorial")
}

root <- normalizePath(".")
if (basename(root) == "revision_checks") {
  root <- dirname(root)
  setwd(root)
}

qmd_path <- file.path(root, "tutorial_v2.qmd")
text <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
spatial <- sub("^[\\s\\S]*?# Spatial meta-analysis\\n", "", text)
spatial <- sub("\\n# Software and package versions[\\s\\S]*$", "", spatial)

expect_text <- function(pattern, label, fixed = TRUE) {
  if (!grepl(pattern, spatial, fixed = fixed)) stop("Missing: ", label)
}
expect_absent <- function(pattern, label, fixed = TRUE) {
  if (grepl(pattern, spatial, fixed = fixed)) stop("Unexpected: ", label)
}
fmt <- function(x, digits) formatC(x, format = "f", digits = digits)

clean <- read.csv(file.path(
  root, "revision_checks", "reviewer18_influential_effects_outputs",
  "published_cleaned_primary_results.csv"
), stringsAsFactors = FALSE)
profile <- read.csv(file.path(
  root, "revision_checks", "reviewer18_cleaned_primary_outputs",
  "cleaned_profile_summary.csv"
), stringsAsFactors = FALSE)
i2 <- read.csv(file.path(
  root, "revision_checks", "reviewer18_cleaned_primary_outputs",
  "cleaned_generalized_i2.csv"
), stringsAsFactors = FALSE)
gaussian <- read.csv(file.path(
  root, "revision_checks", "cleaned_gaussian_kernel_outputs",
  "combined_spgau_targeted_multistart.csv"
), stringsAsFactors = FALSE)
gaussian_zero <- read.csv(file.path(
  root, "revision_checks", "cleaned_gaussian_kernel_outputs",
  "cleaned_gaussian_identification_summary.csv"
), stringsAsFactors = FALSE)
scholer <- read.csv(file.path(
  root, "revision_checks", "scholer_spatial_audit_outputs",
  "scholer_primary_model_results.csv"
), stringsAsFactors = FALSE)

stopifnot(nrow(clean) == 3L, nrow(scholer) == 3L,
          all(clean$n_effects == 2355L), all(clean$n_studies == 390L),
          all(clean$n_sites == 380L))

expect_text("published_cleaned")
expect_text("2,355 effect sizes")
expect_text("390 studies")
expect_text("380 recorded-coordinate locations")
expect_text("all_spatially_usable")
expect_text("2,361 effect sizes")
expect_text("Sensitivity to influential-effect exclusions")
expect_text("WGS84 ellipsoidal geodesic distance")
expect_absent("EPSG:3857")
expect_absent("display: none")
expect_absent("~ effect_id | const")
expect_text("~ site_id | const")
expect_text("identifiability diagnostics, not formal confidence intervals")
expect_text("tested values with nearly equivalent likelihoods")
expect_text("restricted common target")
expect_text("not the preferred biological model")
expect_text("95% CI")
expect_text("95% CrI")
expect_text("Scholer et al. (2020)")
expect_text("949 effect sizes")
expect_text("205 references")
expect_text("454 recorded-coordinate locations")

for (nm in clean$model) {
  row <- clean[clean$model == nm, ]
  expect_text(fmt(row$mean, 3), paste(nm, "mean"))
  expect_text(fmt(row$AIC_REML, 3), paste(nm, "AIC"))
}
for (nm in scholer$model) {
  row <- scholer[scholer$model == nm, ]
  expect_text(fmt(row$pooled_mean, 3), paste("Scholer", nm, "mean"))
  expect_text(fmt(row$AIC_REML, 3), paste("Scholer", nm, "AIC"))
}

combined <- clean[clean$model == "combined", ]
combined_profile <- profile[profile$model == "combined" & profile$component == "tau2", ]
expect_text(fmt(combined$spatial_variance, 5), "cleaned combined spatial variance")
expect_text(fmt(combined_profile$zero_logLik_loss, 3), "cleaned zero-variance loss")

best_gaussian <- gaussian[gaussian$label == "short_312_km", ]
long_gaussian <- gaussian[gaussian$label == "long_3000_km", ]
expect_text(fmt(best_gaussian$rho_km, 2), "best Gaussian rho")
expect_text(fmt(long_gaussian$rho_km, 2), "long Gaussian rho")
expect_text(fmt(long_gaussian$delta_logLik_from_best, 3), "Gaussian log-likelihood spread")
expect_text(fmt(gaussian_zero$best_vs_tau2_zero_logLik_loss, 3), "Gaussian zero-variance loss")

i2_combined <- i2[i2$model == "combined", ]
expect_text(fmt(i2_combined$I2_percent[i2_combined$component == "spatial"], 3),
            "cleaned combined spatial I2")

if (length(gregexpr("\n\\| Combined \\| 0\\.657", spatial, perl = TRUE)[[1L]]) != 1L) {
  stop("Scholer combined row is duplicated or missing")
}

cat("SPATIAL_TUTORIAL_REVISION_VALIDATED\n")
