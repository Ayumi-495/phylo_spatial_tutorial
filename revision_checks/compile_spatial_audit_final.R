# Compile the finalized spatial-audit results without refitting any model.

checks <- "revision_checks"

grau_exp <- read.csv(file.path(
  checks, "totoro_spatial_audit_outputs", "primary_model_results.csv"
), stringsAsFactors = FALSE)
grau_exp_i2 <- read.csv(file.path(
  checks, "i2_definition_audit_outputs", "grau_i2.csv"
), stringsAsFactors = FALSE)
grau_gau_sp <- read.csv(file.path(
  checks, "gaussian_global_outputs", "spatial_only_spgau_result.csv"
), stringsAsFactors = FALSE)
grau_gau_ms <- read.csv(file.path(
  checks, "gaussian_global_outputs", "targeted_free_refits",
  "targeted_free_refits_compiled.csv"
), stringsAsFactors = FALSE)
grau_gau_i2 <- read.csv(file.path(
  checks, "gaussian_global_outputs", "grau_spgau_i2_results.csv"
), stringsAsFactors = FALSE)
scholer <- read.csv(file.path(
  checks, "scholer_spatial_audit_outputs", "scholer_primary_model_results.csv"
), stringsAsFactors = FALSE)
scholer_i2 <- read.csv(file.path(
  checks, "scholer_spatial_audit_outputs", "scholer_i2_results.csv"
), stringsAsFactors = FALSE)

stopifnot(nrow(grau_exp) == 3L, nrow(grau_gau_sp) == 1L,
          nrow(grau_gau_ms) == 3L, nrow(scholer) == 3L)

i2_value <- function(tab, model, component, solution = NULL) {
  keep <- tab$model == model & tab$component == component
  if (!is.null(solution)) keep <- keep & tab$solution == solution
  value <- tab$I2_percent[keep]
  stopifnot(length(value) == 1L)
  value
}

empty_row <- function() {
  data.frame(
    dataset = character(), kernel = character(), model = character(),
    pooled_mean = numeric(), ci_lb = numeric(), ci_ub = numeric(),
    effect_variance = numeric(), study_variance = numeric(),
    spatial_variance = numeric(), rho_km = numeric(),
    logLik_REML = numeric(), AIC_REML = numeric(),
    I2_total_percent = numeric(), I2_effect_percent = numeric(),
    I2_study_percent = numeric(), I2_spatial_percent = numeric(),
    spatial_identification = character(), optimizer_note = character(),
    stringsAsFactors = FALSE
  )
}

out <- empty_row()
for (nm in c("unstructured_only", "spatial_only", "combined")) {
  x <- grau_exp[grau_exp$model == nm, ]
  out <- rbind(out, data.frame(
    dataset = "Grau-Andres", kernel = if (nm == "unstructured_only") "none" else "SPEXP",
    model = nm, pooled_mean = x$mean, ci_lb = x$ci_lb, ci_ub = x$ci_ub,
    effect_variance = x$sigma2_effect,
    study_variance = if (nm == "spatial_only") NA_real_ else x$sigma2_study,
    spatial_variance = if (nm == "unstructured_only") NA_real_ else x$tau2_spatial,
    rho_km = if (nm == "unstructured_only") NA_real_ else x$rho_km,
    logLik_REML = x$logLik_REML, AIC_REML = x$AIC_REML,
    I2_total_percent = i2_value(grau_exp_i2, nm, "total"),
    I2_effect_percent = i2_value(grau_exp_i2, nm, "effect_size"),
    I2_study_percent = if (nm == "spatial_only") NA_real_ else
      i2_value(grau_exp_i2, nm, "study"),
    I2_spatial_percent = if (nm == "unstructured_only") NA_real_ else
      i2_value(grau_exp_i2, nm, "spatial"),
    spatial_identification = switch(nm,
      unstructured_only = "not applicable",
      spatial_only = "variance away from zero; very short rho poorly resolved",
      combined = "additional spatial component weakly identified; rho broad"
    ),
    optimizer_note = "verified primary fit",
    stringsAsFactors = FALSE
  ))
}

out <- rbind(out, data.frame(
  dataset = "Grau-Andres", kernel = "SPGAU", model = "spatial_only",
  pooled_mean = grau_gau_sp$mean, ci_lb = grau_gau_sp$ci_lb,
  ci_ub = grau_gau_sp$ci_ub,
  effect_variance = grau_gau_sp$iid_effect_variance,
  study_variance = NA_real_, spatial_variance = grau_gau_sp$spatial_variance,
  rho_km = grau_gau_sp$e_folding_range_km,
  logLik_REML = grau_gau_sp$logLik_REML, AIC_REML = grau_gau_sp$AIC_REML,
  I2_total_percent = i2_value(grau_gau_i2, "spatial_only_spgau", "total", "primary"),
  I2_effect_percent = i2_value(grau_gau_i2, "spatial_only_spgau", "effect_size", "primary"),
  I2_study_percent = NA_real_,
  I2_spatial_percent = i2_value(grau_gau_i2, "spatial_only_spgau", "spatial", "primary"),
  spatial_identification = "large variance; near-zero fitted rho is not broad-scale spatial structure",
  optimizer_note = "verified primary fit",
  stringsAsFactors = FALSE
))

best <- grau_gau_ms[which.max(grau_gau_ms$logLik_REML), ]
best_solution <- paste0(best$start_label, "_start")
out <- rbind(out, data.frame(
  dataset = "Grau-Andres", kernel = "SPGAU", model = "combined",
  pooled_mean = best$pooled_mean, ci_lb = best$ci_lb, ci_ub = best$ci_ub,
  effect_variance = best$iid_effect_variance, study_variance = best$study_variance,
  spatial_variance = best$spatial_variance, rho_km = best$rho_km,
  logLik_REML = best$logLik_REML, AIC_REML = best$AIC_REML,
  I2_total_percent = i2_value(grau_gau_i2, "combined_spgau", "total", best_solution),
  I2_effect_percent = i2_value(grau_gau_i2, "combined_spgau", "effect_size", best_solution),
  I2_study_percent = i2_value(grau_gau_i2, "combined_spgau", "study", best_solution),
  I2_spatial_percent = i2_value(grau_gau_i2, "combined_spgau", "spatial", best_solution),
  spatial_identification = "weakly identified additional component; rho optimizer-dependent",
  optimizer_note = paste(
    "best observed targeted solution; second stationary solution at rho",
    sprintf("%.1f km has delta logLik %.3f", grau_gau_ms$rho_km[grau_gau_ms$start_label == "primary3091"],
            grau_gau_ms$delta_logLik_from_best[grau_gau_ms$start_label == "primary3091"])
  ),
  stringsAsFactors = FALSE
))

for (nm in c("unstructured_only", "spatial_only", "combined")) {
  x <- scholer[scholer$model == nm, ]
  out <- rbind(out, data.frame(
    dataset = "Scholer", kernel = if (nm == "unstructured_only") "none" else "SPEXP",
    model = nm, pooled_mean = x$pooled_mean, ci_lb = x$ci_lb, ci_ub = x$ci_ub,
    effect_variance = x$iid_effect_variance,
    study_variance = if (nm == "spatial_only") NA_real_ else x$study_variance,
    spatial_variance = if (nm == "unstructured_only") NA_real_ else x$spatial_variance,
    rho_km = if (nm == "unstructured_only") NA_real_ else x$rho_km,
    logLik_REML = x$logLik_REML, AIC_REML = x$AIC_REML,
    I2_total_percent = i2_value(scholer_i2, nm, "total"),
    I2_effect_percent = i2_value(scholer_i2, nm, "effect_size"),
    I2_study_percent = if (nm == "spatial_only") NA_real_ else
      i2_value(scholer_i2, nm, "study"),
    I2_spatial_percent = if (nm == "unstructured_only") NA_real_ else
      i2_value(scholer_i2, nm, "spatial"),
    spatial_identification = switch(nm,
      unstructured_only = "not applicable",
      spatial_only = "variance and rho clearly identified within restricted model",
      combined = "small estimated variance and weakly identified; rho broad"
    ),
    optimizer_note = "primary solution agrees with profile-grid maximum",
    stringsAsFactors = FALSE
  ))
}

write.csv(out, file.path(checks, "spatial_audit_final_comparison.csv"), row.names = FALSE)
cat("SPATIAL_AUDIT_FINAL_TABLE_COMPILED\n")
