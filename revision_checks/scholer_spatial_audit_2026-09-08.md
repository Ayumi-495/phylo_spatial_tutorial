# Scholer spatial `metafor` audit

This audit uses the full, audited Scholer dataset: 949 effect sizes, 205 references, and 454 recorded-coordinate locations. No QMD, manuscript, OU, `brms`, or `glmmTMB` content was changed.

The prepared data use `vi = se^2`, an intercept-only fixed effect, an iid effect-size random term, and WGS84 ellipsoidal great-circle distances in kilometres from `geosphere::distGeo()`. The distance matrix was built over the 454 sorted recorded-coordinate `site_id` levels. Its row and column names were asserted against those factor levels immediately before every fit. `effect_id`, `study_id`, and `site_id` are distinct grouping structures; the spatial outer group `const` has one level, so spatial covariance can occur across studies. No iid location intercept was added.

All models used `method = "REML"`, `test = "t"`, and `sparse = TRUE`.

## Primary fits

| Model | Pooled mean (95% CI) | iid effect variance | Study variance | Spatial variance | rho (km) | REML logLik | AIC | Elapsed | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Unstructured-only | 0.669 (0.558, 0.779) | 0.231 | 0.476 | -- | -- | -825.239 | 1656.479 | 0.9 s | completed; no warning |
| Spatial-only | 0.458 (0.350, 0.567) | 0.280 | -- | 0.323 | 168.7 | -928.304 | 1864.607 | 32.5 s | completed; Matrix class-coercion warning |
| Combined | 0.657 (0.536, 0.778) | 0.225 | 0.449 | 0.0204 | 535.6 | -823.545 | 1657.090 | 39.8 s | completed; Matrix class-coercion warning |

`rma.mv()` returned successfully for all three models and emitted no optimizer/convergence warning. The Matrix class-coercion warning in the spatial fits was retained in the result records; it is not an optimizer warning.

The spatial-only model is strongly disfavoured by AIC: it is 208.129 AIC units above the unstructured-only model. The combined model is only 0.611 AIC units above the unstructured-only model. Thus, once study-level heterogeneity is admitted, adding the spatial component does not improve this AIC comparison.

## Incrementally saved profile grid

The two spatial models were profiled by holding either spatial variance or rho fixed, re-optimizing the other parameters, and saving each profile fit independently. Thirty-nine points were run on Totoro with eight concurrent, one-thread R workers. There were no optimizer warnings or errors. The profile grid is a bounded diagnostic grid, not a formal confidence interval.

| Model / parameter | Primary estimate | Grid maximum | Key profile result |
|---|---:|---:|---|
| Spatial-only tau2 | 0.3227 | 0.3227 | tau2 = 0 was 126.498 log-likelihood units below the grid maximum. |
| Spatial-only rho | 168.7 km | 168.7 km | The coarse grid had a distinct maximum at the primary estimate; neighbouring 100 and 300 km points were 2.14 and 2.60 units lower. |
| Combined tau2 | 0.0204 | 0.0200 | tau2 = 0 was only 1.694 log-likelihood units below the maximum; values 0–0.16 were all within 1.92 units on this grid. |
| Combined rho | 535.6 km | 535.6 km | Every tested value from 10 to 12,000 km was within 1.92 log-likelihood units of the maximum. |

The primary combined solution agrees with the profile-grid maximum, and no higher competing local optimum appeared. However, the combined profile shows a broad spatial-variance/range ridge: at very short rho, the estimated spatial variance collapses essentially to zero and the likelihood approaches the unstructured-only model; at much longer rho, higher spatial variance gives similarly close likelihoods. Therefore, the combined rho point estimate of about 536 km is weakly identified and must not be interpreted as a well-defined correlation range.

The spatial-only model has an identifiable spatial term under its deliberately restricted covariance structure, but its much worse AIC and different pooled mean show that it is absorbing heterogeneity that the combined model allocates primarily to the study-level term. It should not be used to claim robust broad-scale spatial autocorrelation.

No broad multi-start analysis was run because the primary/profile comparison did not show a competing likelihood maximum or starting-value dependence.

## Generalized I² audit (added 2026-09-09)

The corrected Scholer fits and profiles above were already present in spatial-completion commit `0e57293`. They were revalidated against the saved prepared data after a later handoff incorrectly described them as not yet fitted. The saved fit objects contain 949 observations, the same `vi` and intercept-only `X`, and the expected distinct `effect_id`, `study_id`, and `site_id` structures. No Scholer model was refitted in this I² pass.

Using the actual saved sampling variances and design matrix,

\[
\widetilde v = \frac{k-p}{\operatorname{tr}(P)}=0.000725865052242,
\quad P=W-WX(X'WX)^{-1}X'W,\quad W=V^{-1}.
\]

For comparison, the arithmetic mean sampling variance is 0.00785229548999 and the simple harmonic mean is 0.000721567012154. The large difference from the arithmetic mean explains why older Scholer I² values based on “mean sampling variance” must not be reused.

| Model | Total I² (%) | Effect-size I² (%) | Study I² (%) | Spatial I² (%) |
|---|---:|---:|---:|---:|
| Unstructured-only | 99.8974 | 32.5854 | 67.3120 | -- |
| Spatial-only | 99.8798 | 46.4345 | -- | 53.4453 |
| Combined | 99.8956 | 32.4334 | 64.5256 | 2.9366 |

These are plug-in variance-partition summaries, not evidence that the corresponding components are precisely identified. In particular:

- the spatial-only spatial I² is large and its variance/rho profiles have clear maxima under that deliberately restricted model, but the model is strongly disfavoured by AIC;
- the combined spatial point variance is small (0.0204), and it is also weakly identified: fixing it to zero costs only 1.694 log-likelihood units and every tested rho from 10 to 12,000 km remains within 1.92 units of the profile maximum;
- the combined spatial I² of 2.94% is the share of typical marginal variance assigned to that fitted spatial random component. It is not a percentage explained by geographic distance, a spatial range, or a pairwise correlation measure.

Exact I² code and outputs are in `scholer_i2_audit.R`, `scholer_i2_sampling_variance.csv`, and `scholer_i2_results.csv`.

## Saved evidence

- `scholer_primary_model_results.csv`: primary estimates, AIC, warnings, and elapsed times.
- `unstructured_only.rds`, `spatial_only.rds`, and `combined.rds`: completed primary fits.
- `profiles/`: one RDS and one CSV for each completed profile point.
- `scholer_profile_points_compiled.csv` and `scholer_profile_summary.csv`: compiled likelihood grid.
- `scholer_model_metadata.txt`: data, distance, and common fitting settings.
- `scholer_i2_sampling_variance.csv`, `scholer_i2_results.csv`, and `scholer_i2_readout.txt`: generalized I² denominator and component calculations.
