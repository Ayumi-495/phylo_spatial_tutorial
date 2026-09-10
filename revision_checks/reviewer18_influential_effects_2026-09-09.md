# Reviewer Comment 18: influential-effect and spatial-model audit

## Scope and result

This audit recovered the influential-effect exclusions specified in the archived
Grau-Andres et al. analysis, created `published_cleaned` and
`all_spatially_usable` datasets, and fitted the three pre-specified global
`metafor` models to `published_cleaned`. It did not edit `tutorial_v2.qmd`, the
manuscript, or the response letter.

The substantive result is robust: removing the recoverable source-publication
exclusions shifts the pooled means only slightly in the negative direction, all
95% CIs remain below zero, and the AIC ranking remains unstructured-only,
combined, then spatial-only. The additional spatial component in the combined
model remains weakly identified.

## Provenance and recovery of the seven named exclusions

The local archived files are byte-for-byte matches to the files currently
served by Dryad for DOI 10.5061/dryad.0vt4b8h6j:

- data SHA-256: `435c87cac364acc3bfb61e978b83c6bed06a5b08a2336c15e5f7427ed8674b45`
- code SHA-256: `9a7aa3cbd14faeb59a8dd6773f52524db2a240637a40616a5cd8d42db80d015c`

The archived R code identifies three abundance effects and two diversity effects
using Cook's-distance thresholds of 0.015, and two fitness effects using a
threshold of 0.06. No new absolute-effect-size cutoff was introduced.

There is one important archive discrepancy. The code names seven exclusions,
but only six are present in the archived 2,363-row CSV. `Launonen_1999-1` is
named in the code but is already absent from the CSV. Its row-level values
therefore cannot be recovered from this archive, and no substitute observation
was invented. The archived code also defines the fitness data as `reg.fit` but
uses the undefined object `reg.fi` in the preliminary fitness model. Thus, the
fitness Cook's-distance calculation cannot be rerun from the script unchanged;
the two intended fitness exclusions are nevertheless explicit in the subsequent
exclusion statements.

| Effect ID | Study ID | Response | Hedges' d | Sampling variance | SD imputed | Source field | Among two missing-coordinate records? |
|---|---|---:|---:|---:|---|---|---|
| `Ngugi_2022-2` | `Ngugi_2022` | abundance | 5.412367 | 0.517968 | no | Fig5 | no |
| `Gagnon_2015-2` | `Gagnon_2015` | abundance | -9.375232 | 0.799125 | no | Fig5 | no |
| `Moris_2017-1` | `Moris_2017` | abundance | 4.447277 | 0.106241 | no | Fig3a | no |
| `Schwilk_1997-1` | `Schwilk_1997` | diversity | 3.359365 | 1.607111 | no | Fig4 | no |
| `Silveira_2016-4` | `Silveira_2016` | diversity | -2.432012 | 0.217417 | no | Fig2 | no |
| `Launonen_1999-1` | `Launonen_1999` | fitness | unavailable | unavailable | unavailable | unavailable | no: absent from the archived CSV, rather than one of those two records |
| `Ansley_2015-1` | `Ansley_2015` | fitness | 7.435583 | 1.862226 | no | Fig4 | no |

The paper describes the imputation flag in terms of missing standard deviations;
all six recoverable records have `imputed = no`. The two records excluded from
the spatial data because coordinates are missing are `Pellegrini_2021-1` and
`Pellegrini_2021-2`; neither is among the seven code-named influential effects.

## Are these simply the largest absolute effects?

No. Among the 2,361 spatially usable records, the six recoverable exclusions
rank 17 (`Gagnon_2015-2`), 30 (`Ansley_2015-1`), 54 (`Ngugi_2022-2`), 85
(`Moris_2017-1`), 149 (`Schwilk_1997-1`), and 248 (`Silveira_2016-4`) by
descending absolute Hedges' d. None is among the seven largest absolute effects.
`Launonen_1999-1` cannot be ranked because it is absent from the archived CSV.
This confirms that the published screening was not a simple magnitude cutoff:
Cook's distance also reflects model leverage and precision. It does not imply
that a numerically extreme effect is a data error.

## Audited datasets

| Dataset | Effect sizes | Studies | Recorded-coordinate locations |
|---|---:|---:|---:|
| `all_spatially_usable` | 2,361 | 393 | 383 |
| `published_cleaned` | 2,355 | 390 | 380 |

Operationally, `published_cleaned` is the 2,361-row spatial dataset minus the
six code-named exclusions that are present in the archived CSV. The seventh
named effect was already absent before spatial filtering.

## Cleaned-data global spatial fits

All models used Hedges' d, known diagonal sampling variances, an intercept-only
fixed effect, REML with `test = "t"`, and iid effect-size heterogeneity. Spatial
models used a 380 by 380 WGS84 ellipsoidal geodesic distance matrix in kilometres,
with row and column names asserted identical to sorted `site_id` factor levels.
The spatial outer grouping factor was constant across all studies, and no iid
site intercept was added.

| Model | Mean [95% CI] | Effect variance | Study variance | Spatial variance | rho (km) | REML logLik | REML AIC | Status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Unstructured-only | -0.364 [-0.484, -0.243] | 0.752 | 1.138 | -- | -- | -3969.188 | 7944.375 | completed; no fit warning |
| Spatial-only | -0.356 [-0.479, -0.233] | 0.763 | -- | 1.149 | 0.137 | -3980.247 | 7968.495 | completed; no fit warning |
| Combined | -0.362 [-0.545, -0.180] | 0.751 | 1.094 | 0.04785 | 2058.907 | -3968.254 | 7946.509 | completed; Matrix S4 deprecation warning only |

The recorded combined-model warning (`Setting class(x) to multiple strings ...`)
is a Matrix/S4 deprecation warning, not an optimizer-convergence warning. The
`rma.mv` objects do not expose an explicit optimizer status in the inspected
fields, so status is reported conservatively as completed without an explicit
optimizer status rather than as guaranteed convergence.

## Full-data sensitivity comparison

| Dataset | Model | Mean | Effect var. | Study var. | Spatial var. | rho (km) | AIC | Delta AIC within dataset |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| All spatially usable | Unstructured-only | -0.345 | 0.775 | 1.242 | -- | -- | 8042.979 | 0.000 |
| All spatially usable | Spatial-only | -0.334 | 0.793 | -- | 1.231 | 0.050 | 8076.949 | 33.970 |
| All spatially usable | Combined | -0.349 | 0.775 | 1.167 | 0.06188 | 384.958 | 8045.604 | 2.625 |
| Published-cleaned | Unstructured-only | -0.364 | 0.752 | 1.138 | -- | -- | 7944.375 | 0.000 |
| Published-cleaned | Spatial-only | -0.356 | 0.763 | -- | 1.149 | 0.137 | 7968.495 | 24.120 |
| Published-cleaned | Combined | -0.362 | 0.751 | 1.094 | 0.04785 | 2058.907 | 7946.509 | 2.134 |

AIC and log-likelihood values are compared only among models fitted to the same
dataset; their absolute values should not be compared across datasets with
different observations.

- **pooled biological conclusion:** unchanged. Cleaned means are 0.013 to 0.022
  more negative, and every CI remains wholly below zero.
- **model ranking:** unchanged. Unstructured-only has the lowest AIC, combined
  is second, and spatial-only is clearly worst under both datasets.
- **variance allocation:** modestly lower after screening. Effect-size variance
  falls by about 0.024 to 0.030; study variance falls by 0.103 in the
  unstructured model and 0.072 in the combined model; spatial variance falls by
  0.082 in the spatial-only model and 0.014 in the combined model.
- **spatial component:** the spatial-only model still allocates substantial
  variance to an extremely short-range term, while the combined model allocates
  little additional variance to space once study heterogeneity is present.

The cleaned combined point estimate for rho moves from 385 km to 2,059 km, but
this is not evidence for a newly resolved scale. It lies within the approximately
200--3,000 km near-flat region in the existing full-data profile. Moreover,
fixing spatial variance to zero reduces the combined model to the cleaned
unstructured-only model, whose log-likelihood is only 0.933 units below the
cleaned combined maximum. Consequently zero spatial variance remains close to
the maximum and rho cannot be strongly identified jointly. This directly
preserves the current identifiability interpretation, so no expensive profile
grid or targeted refit was needed for this sensitivity audit.

## Spain regional subset

None of the six recoverable exclusions is in Spain, and the seventh named effect
is absent from the Dryad CSV. The 186-observation Spain subset is therefore
unchanged, so the existing cross-package audit is unaffected and no regional
models need to be rerun for Reviewer Comment 18.

## Saved evidence

- Reproducible audit: `revision_checks/reviewer18_influential_effects.R`
- Comparison compiler: `revision_checks/reviewer18_compile_comparison.R`
- Independent checks: `revision_checks/reviewer18_validate_outputs.R`
- Machine-readable outputs and fitted objects:
  `revision_checks/reviewer18_influential_effects_outputs/`

