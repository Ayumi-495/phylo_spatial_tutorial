# Full-dataset spatial-model audit

Status: completed on Totoro; no tutorial QMD numerical result or interpretation has been changed.

## Specification verified before fitting

- Data: 2,361 effect sizes, 393 study IDs, and 383 unique recorded-coordinate `site_id`s.
- Sampling errors: `V = var_Hedges` (diagonal; independent-sampling-error assumption retained).
- Fixed effects: intercept only in every model.
- Shared effect-size heterogeneity: `~ 1 | effect_id` in every model.
- Spatial locations: WGS84 ellipsoidal great-circle distances in kilometres, generated with `geosphere::distm(..., fun = geosphere::distGeo) / 1000`.
- Matrix checks: distance-matrix row and column names were explicitly asserted identical to the ordered `site_id` factor levels before every primary fit and every distributed profile point.
- Spatial outer group: a one-level `const` factor, allowing spatial covariance across studies.
- No iid location intercept was added.

## Primary REML fits

| Model | Pooled mean (95% CI) | Effect-size variance | Study variance | Spatial variance | rho (km) | logLik | AIC | Elapsed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Unstructured-only | -0.345 (-0.470, -0.221) | 0.775 | 1.242 | -- | -- | -4018.489 | 8042.979 | 2.1 s |
| Spatial-only | -0.334 (-0.460, -0.208) | 0.793 | -- | 1.231 | 0.050 | -4034.475 | 8076.949 | 408.8 s |
| Combined | -0.349 (-0.491, -0.206) | 0.775 | 1.167 | 0.0619 | 384.96 | -4017.802 | 8045.604 | 876.7 s |

`rma.mv` did not expose an explicit optimiser-status code in its saved fit object. All three primary fits completed without an optimiser warning. The combined fit emitted a Matrix class-coercion warning; it is retained in the raw record and is not an optimiser convergence warning.

## Profile-likelihood audit

Forty-eight independent fixed-parameter REML refits completed and were saved individually (zero errors): 10 spatial-only `tau2` points, 12 spatial-only `rho` points, 11 combined `tau2` points, and 15 combined `rho` points.

- Spatial-only `tau2` is clearly away from zero: fixing it to zero lowers the REML log-likelihood by 378.25 units relative to the profiled maximum. Its `rho` is, however, weakly resolved at very short distances: 0.03, 0.05, 0.08, and 0.15 km differ by at most 0.04 log-likelihood units; the grid maximum is 0.15 km. Thus the primary 0.050 km estimate is not a hard numerical boundary, but it is not precisely estimated.
- In the combined model, the spatial component is weakly identified. The profile maximum is `tau2 = 0.06` and matches the primary solution. Yet `tau2 = 0` is only 0.687 log-likelihood units lower, and values from 0.02 through 1.6 remain within 0.60 units of the maximum.
- The combined `rho` profile also matches the primary solution at 385 km, but is broad: 200, 385, 700, 1,500, and 3,000 km are all within 0.11 log-likelihood units of the profile maximum. At 0.03--30 km, the refit collapses to the unstructured-only solution; 6,000 and 12,000 km are also only 0.31 and 0.48 units below the maximum.

## Interpretation for the next decision gate

The unstructured-only model has the lowest AIC. The spatial-only model is strongly disfavoured by AIC. The combined model gives a negligible log-likelihood increase over unstructured-only (0.687) but costs two additional parameters (AIC +2.625). The primary combined solution is reproducible through the profiles, so there is no evidence here for a distinct competing local optimum that warrants a multi-start analysis. However, its spatial variance and range are too weakly identified to support a strong spatial-structure interpretation without a later, explicitly approved sensitivity analysis.

Raw models, primary result record, distance matrix, profile-point records, and the compiled profile table are in `totoro_spatial_audit_outputs/`.
