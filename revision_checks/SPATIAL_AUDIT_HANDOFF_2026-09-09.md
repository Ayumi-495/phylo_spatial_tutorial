# Spatial audit handoff: authoritative record before tutorial revision

## A. Current status

- **Branch:** `revision-online-tutorial`.
- **HEAD before this audit-only checkpoint:** `0e572932d6cd6716c7d6eacd8a02f5fe61b46aa7` (`0e57293`, *Complete audited spatial tutorial revision*).
- **Revision worktree:** `/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/`.
- **Totoro staging/output path:** `/home/amizuno/phylo_spatial_tutorial_revision/revision_checks/`. The completed cleaned-primary profiles and companion report were recovered from there into the local revision worktree in this checkpoint.
- **Intentionally unchanged in this checkpoint:** `tutorial_v2.qmd`, the manuscript, and the response letter. This commit contains audit material and this handoff only.
- **OU:** no OU work was started as part of this audit. Do not infer any OU change from the spatial audit files.

The spatial numerical record is now settled for tutorial-writing purposes. No additional spatial model should be run unless a genuine implementation error is found or Ayumi explicitly asks for a new analysis.

## B. Source-of-truth hierarchy

When numerical values disagree, use this order:

1. validated saved model objects and compiled machine-readable CSV outputs;
2. final audit reports derived from those objects;
3. this handoff;
4. historical output embedded in `tutorial_v2.qmd`.

Historical tutorial output is not numerical authority. In particular, values printed before the Reviewer 18 cleaned-data audit, old global projected implementations, and earlier Scholer grouping/output must not be copied into a revision without checking the files below.

The `.rds` model objects named below are locally retained source artifacts and are intentionally git-ignored because of their size. This audit-only commit preserves their scripts, CSV summaries, reports, and validation evidence; it does not claim that a fresh clone contains every fitted object.

## C. Grau-Andres data provenance and dataset definitions

### Raw public data

The audited Dryad data file is `data/Roger_etal_2024/Roger_etal_2024.csv`; the archived analysis code is `data/Roger_etal_2024/global_fire_analysis_240325.R`. They match the Dryad archive for DOI `10.5061/dryad.0vt4b8h6j` exactly:

- data SHA-256: `435c87cac364acc3bfb61e978b83c6bed06a5b08a2336c15e5f7427ed8674b45`;
- code SHA-256: `9a7aa3cbd14faeb59a8dd6773f52524db2a240637a40616a5cd8d42db80d015c`.

The verified hashes are in `revision_checks/reviewer18_influential_effects_outputs/dryad_file_hashes.csv`. Use the local files above, not a separately downloaded derivative, when reproducing the published screening decision.

### `all_spatially_usable`

`all_spatially_usable` has **2,361 effects, 393 studies, and 383 recorded-coordinate locations**. It excludes only `Pellegrini_2021-1` and `Pellegrini_2021-2`, which lack usable coordinates. The saved CSV is `revision_checks/reviewer18_influential_effects_outputs/all_spatially_usable.csv`; reproducible construction and the full-data/global fit workflow are in `revision_checks/reviewer18_influential_effects.R` and the earlier spatial audit material.

### `published_cleaned`

`published_cleaned` has **2,355 effects, 390 studies, and 380 recorded-coordinate locations**. It starts from `all_spatially_usable` and removes the six recoverable influential records named in the published archived analysis:

`Ngugi_2022-2`, `Gagnon_2015-2`, `Moris_2017-1`, `Schwilk_1997-1`, `Silveira_2016-4`, and `Ansley_2015-1`.

The original code names seven exclusions, but `Launonen_1999-1` is already absent from the public 2,363-row CSV. It cannot be removed again, and no substitute record was invented. The full table, including response, Hedges' d, variance, imputation flag, and source field, is `revision_checks/reviewer18_influential_effects_outputs/published_influential_effects.csv`.

The two coordinate-missing records are different observations from the influential records. None of the recoverable influential records occurs in the Spain regional subset; see `revision_checks/reviewer18_influential_effects_outputs/spain_exclusion_overlap.csv`.

## D. Reviewer 18 audit files

- `revision_checks/reviewer18_influential_effects_2026-09-09.md`: narrative recovery of the published Cook's-distance exclusions, cleaned primary models, and full-data sensitivity interpretation. Use it for the reviewer-facing rationale.
- `revision_checks/reviewer18_influential_effects.R`: reproducible data-screening and cleaned-primary fitting script. It checks the Dryad hashes and asserts distance/site-factor alignment.
- `revision_checks/reviewer18_influential_effects_outputs/`: source-of-truth R18 data products and saved cleaned primary fits.
- `.../published_cleaned.csv` and `.../all_spatially_usable.csv`: the two datasets to distinguish in the tutorial.
- `.../published_cleaned_prepared.rds`, `.../published_cleaned_site_lookup.csv`, and `.../published_cleaned_distance_km.csv`: prepared cleaned data, sorted 380-site lookup, and WGS84 ellipsoidal geodesic distance matrix in kilometres.
- `.../unstructured_only.rds`, `.../spatial_only.rds`, and `.../combined.rds`: the three saved cleaned primary `metafor` fits.
- `.../published_cleaned_primary_results.csv`: authoritative numerical table for the cleaned primary fits.
- `.../full_vs_published_cleaned_comparison.csv` and `.../cleaned_minus_full_deltas.csv`: the full versus cleaned sensitivity comparison. Compare AIC only within a dataset.
- `revision_checks/reviewer18_cleaned_i2.R`: reads saved cleaned objects only and calculates generalized I2 from the actual cleaned `V` and `X`.
- `revision_checks/reviewer18_cleaned_primary_outputs/cleaned_generalized_i2.csv` and `.../cleaned_sampling_variance_summary.csv`: authoritative cleaned generalized-I2 components and its sampling/design denominator.
- `.../profile_grid.csv`, `.../profiles/`, and `.../cleaned_profile_results_compiled.csv`: the targeted 43-point cleaned spatial-parameter profile grid, incrementally saved point fits, and compiled result table.
- `.../cleaned_profile_summary.csv`: concise profile maxima, zero-variance losses, and near-equivalent grid ranges.
- `.../targeted_multistart/`, `.../targeted_multistart_compiled.csv`, and `.../targeted_multistart_summary.csv`: nine limited checks prompted by the spatial-only fixed-variance profile, not a broad multi-start search.
- `.../cleaned_vs_full_generalized_i2.csv` and `.../cleaned_vs_full_aic_ranks.csv`: machine-readable cleaned/full I2 and within-dataset AIC-rank comparisons.
- `revision_checks/reviewer18_cleaned_primary_2026-09-09.md`: companion narrative for the cleaned-primary I2/profile support.
- `revision_checks/reviewer18_cleaned_collect_profiles.R` and `revision_checks/reviewer18_cleaned_validate.R`: compilation and five validation stages; neither refits a model.
- `revision_checks/GATES_reviewer18_cleaned_primary.md`: the audit-checkpoint gate ledger. It records validation evidence after this commit.

## E. Final primary Grau-Andres SPEXP results

All models use `published_cleaned` (2,355 effects, 390 studies, 380 recorded-coordinate locations), Hedges' d, known diagonal `V = diag(var_Hedges)`, intercept-only fixed effect, REML with `test = "t"`, iid effect-size heterogeneity, and WGS84 ellipsoidal geodesic distances in kilometres. In spatial fits, `site_id` is a distinct grouping factor and the outer group `const` has one level, so covariance can occur across studies. No iid site intercept is included.

| Model | Mean (95% CI) | Effect variance | Study variance | Spatial variance | rho km | REML logLik | AIC | Status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Unstructured-only | -0.363934 (-0.484464, -0.243404) | 0.751560 | 1.138302 | -- | -- | -3969.187514 | 7944.375029 | completed; no fit warning |
| Spatial-only | -0.355840 (-0.478525, -0.233154) | 0.763024 | -- | 1.149481 | 0.136939 | -3980.247274 | 7968.494547 | completed; no fit warning |
| Combined | -0.362424 (-0.544630, -0.180219) | 0.751387 | 1.094252 | 0.047846 | 2058.906947 | -3968.254342 | 7946.508684 | completed; Matrix/S4 deprecation warning only |

The combined warning is a matrix-class/S4 deprecation warning, not an optimizer failure. The saved `rma.mv` objects do not provide an explicit optimizer status, so the reported status is deliberately conservative.

## F. Cleaned-primary identifiability results

- **Spatial-only spatial variance:** fixing `tau2 = 0` is **381.047** REML log-likelihood units below the maximum. The variance is separated from zero under this deliberately restricted model.
- **Spatial-only rho:** the free solution is **0.1369 km**, while the profile-grid maximum is near **0.05 km** and differs by only **0.0033** log-likelihood units. The short-range likelihood is poorly resolved; targeted grid values from approximately 0.005 to 0.5 km remain near-equivalent.
- **Combined spatial variance:** fixing `tau2 = 0` is only **0.933** log-likelihood units below the combined maximum. The additional spatial variance is weakly identified, even though its point estimate is numerically small.
- **Combined rho:** the profile is broad. Approximately **200 to 12,000 km** are within the targeted near-equivalent support region. The grid maximum near **385 km** and free estimate near **2,059 km** differ by only **0.017** log-likelihood units.
- **Targeted multi-start:** the largest improvement was **0.0196** log-likelihood units. It found no substantively competing optimum.

Interpretation boundary: a small point estimate for spatial variance and weak identification are different claims. Do not present any combined-model rho as a well-determined biological correlation distance.

## G. Cleaned-primary generalized I2

The generalized representative sampling variance is **`v_tilde = 0.110937236`**, calculated from the cleaned actual `V` and intercept-only `X` as `(k - p) / tr(P)`.

| Model | Effect I2 | Study I2 | Spatial I2 | Total I2 |
|---|---:|---:|---:|---:|
| Unstructured-only | 37.563% | 56.892% | -- | 94.455% |
| Spatial-only | 37.709% | -- | 56.808% | 94.517% |
| Combined | 37.486% | 54.592% | 2.387% | 94.465% |

Component I2 is **marginal variance allocation**. It is not variance explained by geographic distance, the strength of pairwise spatial correlation, a spatial range, or proof that a component is precisely identified. Use `reviewer18_cleaned_i2.R`, `cleaned_generalized_i2.csv`, and `cleaned_sampling_variance_summary.csv` for the calculation and numerator/denominator provenance.

## H. Full-data Reviewer 18 sensitivity analysis

`all_spatially_usable` remains the **2,361-effect sensitivity analysis**, not the primary worked example. Removing the six recoverable source-publication exclusions changes means only slightly and all confidence intervals remain below zero. The AIC ordering is unchanged within both datasets:

1. unstructured-only;
2. combined;
3. spatial-only.

Total I2 remains very similar. In the combined model, spatial I2 changes from approximately **2.93%** in the full data to **2.39%** in `published_cleaned`; spatial variance/rho identifiability conclusions do not change. The tutorial should frame this as a Reviewer 18 influential-effect sensitivity: it supports robustness of the pooled biological conclusion and covariance-model ranking, but not a resolved spatial range.

Do not compare absolute AIC values between the 2,355- and 2,361-effect datasets. AIC comparisons are valid only among models fitted to the same observations.

## I. Gaussian full-data sensitivity

The finalized Gaussian audit is `revision_checks/grau_global_gaussian_kernel_audit_2026-09-09.md`; its saved objects, profiles, targeted refits, WGS84 distance matrix, and generalized I2 are under `revision_checks/gaussian_global_outputs/`. The script is `revision_checks/grau_global_gaussian_kernel_audit.R`.

This is a **full 2,361-effect kernel/optimizer sensitivity analysis**, not a cleaned-primary analysis. Do not rerun Gaussian fits on `published_cleaned` without a separate decision.

- Gaussian spatial-only has large spatial variance and an extremely short 0.362-km e-folding range.
- The Gaussian combined fit has optimizer-dependent stationary solutions around **312 km** and **3,091 km**; their log-likelihood difference is only **0.384**.
- Relative to the best observed combined solution, zero spatial variance costs only **1.063** log-likelihood units.
- The Gaussian rho is practically weakly identified. Neither rho should receive substantive biological interpretation.

## J. Spain cross-package analysis

Spain is the regional cross-package implementation example, not the global primary model comparison. Its source report is `revision_checks/regional_cross_package_audit_2026-09-08.md`, and the staged implementation is `revision_checks/regional_cross_package_audit.R`.

- **Data:** 186 effects, 30 studies, 32 recorded-coordinate locations; none of the R18 influential exclusions is present.
- **Geometry:** Spain-specific WGS84 Lambert Conformal Conic projection (`+proj=lcc +lat_1=38 +lat_2=43 +lat_0=40.5 +lon_0=-3.5 +datum=WGS84 +units=m +no_defs`), expressed in kilometres. Maximum pairwise distortion versus WGS84 geodesic distance is approximately **0.103%**.
- **Target model:** matched exponential spatial-only covariance, with known diagonal sampling variance, intercept-only fixed effect, and iid effect-size heterogeneity. No study random intercept was included for this cross-package demonstration.
- **Inputs and geometry:** `revision_checks/regional_cross_package_audit_outputs/spain_data.rds`, `spain_prepared.rds`, `spain_site_lookup_projected_km.csv`, `spain_distance_lcc_km.csv`, `spain_distance_great_circle_km.csv`, and `spain_projection_distortion.csv`.
- **Frequentist fits:** `metafor_spatial_only.rds`/`metafor_spatial_only_result.csv` and `glmmTMB_spatial_only.rds`/`glmmTMB_spatial_only_result.csv`, with diagnostics in `glmmTMB_diagnostics.txt` and parameter conversion in `glmmTMB_theta.csv`.
- **Bayesian fit:** `brms_output/brms_spatial_only.rds`, `brms_output/brms_posterior_summary.csv`, `brms_output/brms_diagnostics.csv`, `brms_output/brms_pp_check_dens_overlay.png`, and `brms_output/brms_stancode.stan`.
- **Comparison:** `regional_cross_package_comparison.csv` and `regional_cross_package_comparison_notes.txt` are the compact authoritative comparison tables.

The matched exponential `metafor` and `glmmTMB` estimates agree to numerical precision. The `brms` fit uses 4 chains, seed 20260908, 10 Stan threads per chain, explicit priors, and passed the saved R-hat/ESS/divergence/treedepth/PPC checks. Use the authoritative I2 files `revision_checks/i2_definition_audit_outputs/spain_i2_metafor_glmmTMB.csv` and `spain_i2_brms_posterior.csv` with the definition report `revision_checks/i2_definition_audit_2026-09-09.md`.

Confidence intervals and credible intervals are not interchangeable. Their numerical similarity is a cross-package consistency check only, not equivalence of frequentist and Bayesian inference.

## K. Scholer Example 2

Scholer is the advanced cross-classified empirical example, not a second package-comparison exercise.

- **Prepared data and structure:** `revision_checks/scholer_structure_audit_outputs/scholer_prepared.rds`, `scholer_site_lookup.csv`, and `scholer_distance_great_circle_km.csv`.
- **Data:** 949 effects, 205 references/studies, and 454 recorded-coordinate locations, using WGS84 ellipsoidal geodesic distances in kilometres.
- **Fits and profiles:** `revision_checks/scholer_spatial_audit_outputs/unstructured_only.rds`, `spatial_only.rds`, `combined.rds`, `scholer_primary_model_results.csv`, `profiles/`, `scholer_profile_points_compiled.csv`, and `scholer_profile_summary.csv`.
- **Generalized I2:** `scholer_i2_audit.R`, `scholer_i2_sampling_variance.csv`, and `scholer_i2_results.csv`; the finalized narrative is `revision_checks/scholer_spatial_audit_2026-09-08.md`.

The central result is: spatial-only is strongly disfavoured by AIC and shifts the pooled mean; combined is close to unstructured-only; combined spatial variance is small and weakly identified; combined rho is effectively unresolved across a very broad range; combined spatial I2 is approximately **2.94%**.

Do not reuse the historical tutorial code that paired the 454 by 454 site distance matrix with `effect_id`. The correct spatial grouping is `site_id`, distinct from both `effect_id` and `study_id`.

## L. Files and content that must not be reused

- Historical numerical output in `tutorial_v2.qmd` when it disagrees with the audit objects/CSVs.
- The pre-R18 full-data Grau analysis presented as though it were the primary analysis. `published_cleaned` is primary; full data is sensitivity.
- Historical global projected `brms`/`glmmTMB` implementations. The global primary uses great-circle distances in `metafor`; the projected three-package comparison is Spain only.
- Arithmetic-mean or otherwise inconsistent I2 values. Use the generalized `v_tilde` outputs above.
- Historical Scholer model/output that mismatches the spatial matrix and grouping factor.
- Any isolated rho point estimate without its variance, model structure, and profile-identifiability qualification.

## M. Intended tutorial-writing structure

This is a handoff recommendation only. It does not authorise tutorial edits in this audit checkpoint.

1. Explain spatial data structure and geographic distance.
2. Present the global Grau-Andres `published_cleaned` primary SPEXP analysis.
3. Present the Reviewer 18 full-data influential-effect sensitivity.
4. Explain generalized I2 as marginal variance allocation.
5. Present the full-data Gaussian kernel/optimizer sensitivity.
6. Present the Spain regional cross-package implementation.
7. Present Scholer as an advanced cross-classified spatial meta-analysis.
8. End with a short reporting checklist: data definition, distance/unit, model structure, variance/rho profiles, and an explicit distinction between pooled robustness and spatial-structure uncertainty.

## N. Remaining limitations and unresolved issues

- All completed Totoro cleaned-primary profile points, compiled CSVs, targeted multi-start outputs, and the companion report were recovered locally in this checkpoint. No requested cleaned-primary output remains Totoro-only.
- A no-execute Quarto render previously confirmed the intended spatial sections are visible. An executing full render remains blocked upstream by the pre-existing missing `Rdata/tutorial_v2/moura2021_BM_meta_reg.rds`; do not solve that by altering non-spatial content.
- The remaining work is editorial, not numerical: decide the amount of detail and reader guidance to put around the identifiability results, then update `tutorial_v2.qmd` from these source files. Do not report a biological spatial range from the combined models.

## O. Checklist before any tutorial edit

- Confirm branch and HEAD.
- Read this handoff before using historical QMD output.
- Inspect the named reports and CSVs/models, especially `published_cleaned_primary_results.csv`, `cleaned_profile_summary.csv`, and `cleaned_generalized_i2.csv`.
- Confirm `published_cleaned` remains the primary global dataset and `all_spatially_usable` remains the R18 sensitivity analysis.
- Confirm Spain is unaffected by R18.
- Confirm that no additional spatial fitting is required.
- Preserve the distinctions among primary analysis, influential-effect sensitivity, kernel/optimizer sensitivity, and regional cross-package implementation.
