# Spain regional cross-package audit

This audit is separate from `tutorial_v2.qmd`. No QMD, manuscript, or response-letter edits were made. The full Spain-labelled subset was retained: 186 effect sizes, 30 studies, and 32 recorded-coordinate locations. The geographically isolated singleton was not removed.

## Implementation check before fitting

The current local versions were `metafor 5.0.1`, `glmmTMB 1.1.15`, `sf 1.1.1`, and `geosphere 1.6.8`.

`glmmTMB::equalto()` was verified with a 20-effect-size preflight using the actual `var_Hedges` values. With `dispformula = ~0`, the known diagonal sampling-variance model agreed with the corresponding `metafor` model:

| Check | metafor | glmmTMB |
|---|---:|---:|
| Pooled mean | 0.8708312 | 0.8708312 |
| REML log-likelihood | -30.3410719 | -30.3410718 |
| Hessian | -- | positive definite |

The `equalto()` matrix rows and columns were explicitly matched to the effect-ID factor levels. The current glmmTMB implementation therefore removes the previous sampling-variance limitation; the remaining package difference is spatial geometry and kernel parameterisation.

## Regional geometry

Coordinates were transformed from WGS84 longitude/latitude to:

```
+proj=lcc +lat_1=38 +lat_2=43 +lat_0=40.5 +lon_0=-3.5 +datum=WGS84 +units=m +no_defs
```

Projected x/y coordinates were divided by 1,000 and used as kilometres. The 32-site Euclidean distances were compared with WGS84 ellipsoidal great-circle distances from `geosphere::distm(..., fun = geosphere::distGeo) / 1000`:

- maximum absolute pairwise distortion: 0.1031%
- 95th percentile absolute distortion: 0.0883%
- maximum great-circle separation: 998.33 km

The distance and coordinate orders were asserted against the site/effect factor levels before fitting.

## Target model

All fitted models used the same 186 effect sizes, intercept-only fixed effect, and independent known diagonal sampling variances. The target covariance was:

1. known sampling-error covariance `V = diag(var_Hedges)`;
2. iid effect-size heterogeneity;
3. exponential spatial covariance over the 32 recorded coordinate locations.

No study-level random intercept was included in this regional comparison.

## Primary fits

For likelihood comparability, the `metafor` fit used `control = list(REMLf = FALSE)`. This is still REML; it uses the likelihood convention that matches the `glmmTMB` `equalto()` implementation. The variance and range estimates are unchanged from the default `metafor` REML constant.

| Package | Mean (95% interval) | iid effect variance | Spatial variance | Spatial SD | rho (km) | REML logLik | AIC | Diagnostics |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `metafor` | -0.0973 (-0.489, 0.295) | 0.2137 | 0.6257 | 0.7910 | 28.96 | -241.320700 | 490.6414 | completed; Matrix class-coercion warning |
| `glmmTMB` | -0.0973 (-0.491, 0.296) | 0.2137 | 0.6257 | 0.7910 | 28.96 | -241.320700 | 490.6414 | convergence code 0; `pdHess = TRUE`; `diagnose()` flags Wald caution |

The estimates and REML/AIC values agree to numerical precision. The small difference in the displayed mean interval reflects `metafor`'s t-based interval (`test = "t"`) versus glmmTMB's Wald z interval. The glmmTMB `diagnose()` message concerns an unusually large absolute z statistic for the Gaussian dispersion intercept; it is a caution about the Wald approximation, not an optimizer or positive-definite-Hessian failure.

For glmmTMB, the exponential structure is parameterised as

```
correlation(d) = exp(-exp(-theta[2]) * d)
```

so the common e-folding range is `rho = exp(theta[2])` km. The fitted values were `theta[1] = -0.23447` (log spatial SD) and `theta[2] = 3.36601` (log rho).

## `metafor` profile check

Profile likelihoods were saved for the spatial variance and rho. The spatial variance was not on the zero boundary:

- estimate `tau2 = 0.6257`
- profile 95% interval: `0.2558–1.5691`

The range was estimable but substantially less precise:

- estimate `rho = 28.96 km`
- profile 95% interval: approximately `2.90–164.17 km`

Thus the Spain spatial-only model is not degenerate, but it should be presented as a regional demonstration with a moderately broad range estimate rather than as strong evidence for a precisely determined spatial scale.

## Files

- `regional_cross_package_audit.R`: staged audit script
- `spain_prepared.rds`, projected site lookup, and both distance matrices
- `metafor_spatial_only.rds`, `metafor_spatial_only_result.csv`
- `metafor_profile_tau2.csv`, `metafor_profile_rho.csv`, and profile CIs
- `glmmTMB_spatial_only.rds`, `glmmTMB_spatial_only_result.csv`, and `glmmTMB_diagnostics.txt`

The brms fit and diagnostics are described below; it was run on Totoro using four chains and 10 threads per chain (40 requested CPU cores).

## `brms` regional fit (Totoro)

The saved regional `brms` fit used the same 186 effect sizes, 32 projected coordinate locations, intercept-only fixed effect, and projected x/y coordinates in kilometres. Its formula was:

```
d_Hedges | se(sqrt(var_Hedges), sigma = TRUE) ~
  1 + gp(x_km, y_km, cov = "exponential", scale = FALSE)
```

In the generated Stan code, `sigma = TRUE` gives the observation model `Normal(mu, sqrt(se_i^2 + sigma^2))`. Thus `sigma` is the iid effect-size SD and is not duplicated by an additional effect-ID random effect. The exponential GP covariance is `sdgp^2 * exp(-distance / lscale)`; with `scale = FALSE` and coordinates in km, `lscale` is the e-folding range in km. Repeated coordinate pairs are grouped into the same GP input, giving 32 unique latent spatial locations.

The run used `brms 2.22.0`, CmdStan `2.36.0`, four chains, four concurrent chain processes, ten Stan threads per chain (40 requested CPU cores), seed `20260908`, 3,000 iterations per chain with 1,500 warmup iterations, `adapt_delta = 0.95`, and `max_treedepth = 12`. The completed fit was saved before diagnostics. All four chains completed; max R-hat was 1.0023, minimum bulk ESS 1,121, minimum tail ESS 1,708, divergences were 0, and maximum treedepth was 7/12. A density-overlay posterior predictive check was saved as `brms_pp_check_dens_overlay.png`.

Posterior medians and 95% CrIs were: pooled mean `-0.100` (`-0.533, 0.388`); iid effect variance `0.220` (`0.132, 0.352`); spatial variance `0.661` (`0.251, 1.671`); and rho `36.6 km` (`8.1, 155.3 km`). The transformed variance summaries are calculated draw-by-draw, rather than by squaring interval endpoints.

## Compact cross-package comparison

| Package | Mean (interval) | iid effect variance | spatial variance (SD) | rho (km) | Likelihood / diagnostics |
|---|---:|---:|---:|---:|---|
| `metafor` | -0.097 (-0.489, 0.295) | 0.214 | 0.626 (0.791) | 29.0 | REML logLik -241.321; AIC 490.641; Matrix class-coercion warning |
| `glmmTMB` | -0.097 (-0.491, 0.296) | 0.214 | 0.626 (0.791) | 29.0 | REML logLik -241.321; AIC 490.641; convergence 0; pdHess TRUE; `diagnose()` Wald caution |
| `brms` | -0.100 (-0.533, 0.388) | 0.220 (0.132–0.352) | 0.661 (0.251–1.671; SD median 0.813) | 36.6 (8.1–155.3) | max R-hat 1.0023; ESS bulk/tail 1,121/1,708; 0 divergences; treedepth 7/12; PPC saved |

The `metafor` and `glmmTMB` rows are REML fits using the same known diagonal sampling VCV. `metafor` used `REMLf = FALSE` solely to match the `equalto()` likelihood convention; the variance and rho estimates are unchanged. The `brms` row is Bayesian posterior inference with explicit priors and is not assigned a REML log-likelihood or AIC. Its interval is a 95% credible interval, whereas the frequentist intervals are confidence intervals. Similar intervals are therefore only a consistency check, not an equivalence of inferential interpretation.

The close `metafor`/`glmmTMB` agreement verifies the `equalto()` implementation for this regional geometry. The `brms` result is numerically similar but differs slightly because it uses posterior inference and explicit priors. The glmmTMB `diagnose()` flag concerns a Wald approximation for a dispersion component; optimizer convergence and the positive-definite Hessian were satisfactory.

No QMD, manuscript, or response-letter text has been changed in this audit. All fit objects, generated Stan code, diagnostics, posterior-predictive draws, and the compact CSV comparison are kept under `revision_checks/regional_cross_package_audit_outputs/`.
