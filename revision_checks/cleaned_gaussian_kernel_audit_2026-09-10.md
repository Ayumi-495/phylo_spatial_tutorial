# Cleaned-data Gaussian kernel sensitivity audit (2026-09-10)

## Scope

This is the requested fixed-dataset kernel check on the validated Grau-Andres `published_cleaned` object. It reuses the saved prepared input, the WGS84 ellipsoidal geodesic distance matrix in kilometres, and the previously validated effect/study/site hierarchy. It does not reconstruct the R18 screen, refit SPEXP models, run a Gaussian profile grid, or edit tutorial/manuscript/response sources.

The input validation confirms **2,355 effect sizes**, **390 studies**, **380 recorded-coordinate locations**, one constant spatial outer group, and identical distance-matrix row/column names and `site_id` levels. The saved cleaned SPEXP and Gaussian fits use the identical response, diagonal known sampling variances, and intercept-only fixed-effect design.

`SPGAU` uses `Cor(d) = exp(-d^2/rho^2)`, so `rho` is the e-folding distance in kilometres.

## Fits

| Model / solution | Mean (95% CI) | Effect variance | Study variance | Spatial variance | rho (km) | REML logLik | AIC | Status | Elapsed seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| Spatial-only | -0.3577 [-0.4810, -0.2344] | 0.76343 | — | 1.15270 | 0.346 | -3980.1863 | 7968.3726 | completed; no optimizer warning | ~108 |
| Combined, short start | -0.3705 [-0.5060, -0.2350] | 0.75154 | 1.03729 | 0.08403 | 306.47 | -3967.9028 | 7945.8055 | completed; no optimizer warning | ~1,209 |
| Combined, intermediate start | -0.3705 [-0.5060, -0.2350] | 0.75154 | 1.03729 | 0.08402 | 306.48 | -3967.9028 | 7945.8055 | completed; no optimizer warning | ~987 |
| Combined, ~3,000-km start | -0.3525 [-0.5587, -0.1463] | 0.75139 | 1.10393 | 0.05016 | 3,673.79 | -3968.1882 | 7946.3764 | completed; no optimizer warning | ~328 |

Elapsed times are reconstructed to the nearest few seconds from the sequential Totoro output timestamps because the original process ended after the free fits and before it could write its compact summary. The completed fitted objects were preserved and no fit was repeated.

## Targeted identifiability check

The short and intermediate starts converge to the same approximately **306-km** stationary solution. The ~3,000-km start converges instead to a **3,674-km** solution. The latter is only **0.2854** REML log-likelihood units below the best solution (Delta AIC = 0.5708). Thus, the point estimate for the combined Gaussian range is optimizer-dependent.

For `tau2 = 0`, the Gaussian combined specification is mathematically the already saved cleaned unstructured model: the spatial range is nonidentified and no separate numerical fit is warranted. Its REML log-likelihood is -3969.1875, only **1.2847** units below the best Gaussian combined solution. The additional Gaussian spatial component is therefore weakly identified.

No Gaussian profile grid was run: the requested limited starts already reveal competing near-equal range solutions, and the exact zero-variance restriction supplies the needed boundary comparison.

## Generalized I2

Using the already validated cleaned-data `v_tilde = 0.110937235977794`:

- Spatial-only: effect-size I2 = 37.66%, spatial I2 = 56.87%, total I2 = 94.53%.
- Best combined (short) solution: effect-size I2 = 37.88%, study I2 = 52.29%, spatial I2 = 4.24%, total I2 = 94.41%.

These are variance-allocation summaries conditional on the selected solution, not evidence for a uniquely estimated spatial range.

## Interpretation relative to prior audits

The cleaned Gaussian sensitivity supports, rather than materially changes, the full-data Gaussian interpretation. With the dataset held fixed to the primary `published_cleaned` observations:

1. The pooled biological conclusion remains negative across SPEXP and SPGAU models.
2. Spatial-only again assigns substantial heterogeneity to an extremely short-range component; that does not demonstrate broad-scale autocorrelation because the model omits study-level unstructured heterogeneity.
3. Once study-level heterogeneity is included, Gaussian spatial variance is small and weakly identified, and the range has competing short/intermediate and long-range solutions with near-equal likelihood.

The full 2,361-effect Gaussian audit remains useful historical/R18-sensitivity evidence, but this cleaned-data result is the preferred kernel-sensitivity source for future tutorial writing.

## Source artifacts

- Fit script: `revision_checks/cleaned_gaussian_kernel_audit.R`
- No-refit record finalizer: `revision_checks/cleaned_gaussian_kernel_finalize.R`
- Input proof: `revision_checks/cleaned_gaussian_kernel_outputs/cleaned_gaussian_input_validation.csv`
- Fit records: `spatial_only_spgau_result.csv`, `combined_spgau_targeted_multistart.csv`, and `combined_spgau_tau2_zero_result.csv` in that same output directory
- Identification and I2: `cleaned_gaussian_identification_summary.csv`, `cleaned_gaussian_generalized_i2.csv`, and `cleaned_gaussian_validation_summary.csv`
- Saved model objects are local, gitignored `.rds` files in `revision_checks/cleaned_gaussian_kernel_outputs/`.
