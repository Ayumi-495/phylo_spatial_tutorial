# Verify locally referenced paths in SPATIAL_AUDIT_HANDOFF_2026-09-09.md.
# This check does not fit models or edit tutorial/manuscript/response sources.

root <- normalizePath(".")
required_paths <- c(
  "tutorial_v2.qmd",
  "data/Roger_etal_2024/Roger_etal_2024.csv",
  "data/Roger_etal_2024/global_fire_analysis_240325.R",
  "revision_checks/SPATIAL_AUDIT_HANDOFF_2026-09-09.md",
  "revision_checks/reviewer18_influential_effects_2026-09-09.md",
  "revision_checks/reviewer18_influential_effects.R",
  "revision_checks/reviewer18_influential_effects_outputs/dryad_file_hashes.csv",
  "revision_checks/reviewer18_influential_effects_outputs/published_cleaned_primary_results.csv",
  "revision_checks/reviewer18_influential_effects_outputs/full_vs_published_cleaned_comparison.csv",
  "revision_checks/reviewer18_influential_effects_outputs/published_influential_effects.csv",
  "revision_checks/reviewer18_influential_effects_outputs/spain_exclusion_overlap.csv",
  "revision_checks/reviewer18_influential_effects_outputs/published_cleaned_prepared.rds",
  "revision_checks/reviewer18_influential_effects_outputs/published_cleaned_distance_km.csv",
  "revision_checks/reviewer18_influential_effects_outputs/unstructured_only.rds",
  "revision_checks/reviewer18_influential_effects_outputs/spatial_only.rds",
  "revision_checks/reviewer18_influential_effects_outputs/combined.rds",
  "revision_checks/reviewer18_cleaned_i2.R",
  "revision_checks/reviewer18_cleaned_primary_2026-09-09.md",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_generalized_i2.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_sampling_variance_summary.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/profile_grid.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_profile_results_compiled.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_profile_summary.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/targeted_multistart_compiled.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/targeted_multistart_summary.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_vs_full_generalized_i2.csv",
  "revision_checks/reviewer18_cleaned_primary_outputs/cleaned_vs_full_aic_ranks.csv",
  "revision_checks/reviewer18_cleaned_collect_profiles.R",
  "revision_checks/reviewer18_cleaned_validate.R",
  "revision_checks/GATES_reviewer18_cleaned_primary.md",
  "revision_checks/grau_global_gaussian_kernel_audit_2026-09-09.md",
  "revision_checks/grau_global_gaussian_kernel_audit.R",
  "revision_checks/gaussian_global_outputs/combined_spgau.rds",
  "revision_checks/gaussian_global_outputs/targeted_free_refits/targeted_free_refits_compiled.csv",
  "revision_checks/regional_cross_package_audit_2026-09-08.md",
  "revision_checks/regional_cross_package_audit.R",
  "revision_checks/regional_cross_package_audit_outputs/spain_prepared.rds",
  "revision_checks/regional_cross_package_audit_outputs/regional_cross_package_comparison.csv",
  "revision_checks/regional_cross_package_audit_outputs/brms_output/brms_spatial_only.rds",
  "revision_checks/i2_definition_audit_2026-09-09.md",
  "revision_checks/i2_definition_audit_outputs/spain_i2_metafor_glmmTMB.csv",
  "revision_checks/i2_definition_audit_outputs/spain_i2_brms_posterior.csv",
  "revision_checks/scholer_spatial_audit_2026-09-08.md",
  "revision_checks/scholer_structure_audit_outputs/scholer_prepared.rds",
  "revision_checks/scholer_spatial_audit_outputs/scholer_primary_model_results.csv",
  "revision_checks/scholer_spatial_audit_outputs/scholer_profile_summary.csv",
  "revision_checks/scholer_spatial_audit_outputs/scholer_i2_results.csv"
)

missing <- required_paths[!file.exists(file.path(root, required_paths))]
stopifnot(!length(missing))

report <- readLines(file.path(root, "revision_checks", "SPATIAL_AUDIT_HANDOFF_2026-09-09.md"),
                    warn = FALSE)
required_headings <- paste0("## ", LETTERS[1:15], ".")
stopifnot(all(vapply(required_headings,
                     function(x) any(startsWith(report, x)), logical(1))))
cat("SPATIAL_AUDIT_HANDOFF_PATHS_VALIDATED\n")
