# I² definition and calculation audit

Date: 2026-09-09  
Scope: read-only audit of the manuscript/tutorial definitions and calculation from already saved full Grau-Andrés and Spain exponential-model fits. No model was refitted. `tutorial_v2.qmd`, the manuscript, and the response letter were not edited. OU and Scholer I² were deliberately excluded.

## Conclusion

Use the generalized design- and sampling-covariance-based representative sampling variance

\[
\widetilde v = \frac{k-p}{\operatorname{tr}(P)},\qquad
P=W-WX(X'WX)^{-1}X'W,\qquad W=V^{-1}.
\]

For the models audited here, define \(H=\sum_j h_j\), where each \(h_j\) is the fitted *marginal variance* contributed by a random-effect component to an observation. Then report

\[
I^2_{\mathrm{total}}=100\frac{H}{H+\widetilde v},\qquad
I^2_j=100\frac{h_j}{H+\widetilde v}.
\]

This is the natural extension of the multilevel formulation described in the [`metafor` I² guidance](https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate). That guidance derives the general \(P\)-matrix formula, sums multilevel heterogeneity variances, and explicitly extends \(W\) to \(V^{-1}\) when the known sampling-error covariance matrix is non-diagonal.

This recommendation is defensible for the tutorial's iid random intercepts and normalized phylogenetic/spatial correlation structures because each component has a unit diagonal. It is not a universal rule for every correlated random effect; the qualifications below must accompany it.

## 1. Current definitions and inconsistencies

| Location | Current calculation | Representative sampling variance | Combining heterogeneity | Component I² |
|---|---|---|---|---|
| Manuscript, Eq. 11 (manuscript p. 16) | \(100\sum_j\sigma_j^2/(\sum_j\sigma_j^2+\bar v)\) | Described as “mean sampling-error variance”; as written this suggests an arithmetic mean, but the exact operation is not defined | Sum of all components | Yes, reported with the same total-variance denominator |
| Manuscript spatial text (manuscript p. 21) | Reports total 83.6%, effect 32.7%, spatial 50.8% | Again explicitly says “mean sampling variance” | Effect plus spatial | Yes |
| Current `metafor` tutorial | `orchaRd::i2_ml(phylo_eg1_meta_ma_BM)` with no `method` argument | `orchaRd` default `method="ratio"`: \((k-1)\sum w_i / [ (\sum w_i)^2-\sum w_i^2]\), \(w_i=1/v_i\) | Sums `model$sigma2` | Yes, each `sigma2` uses the common denominator |
| `orchaRd::i2_ml(method="matrix")` | Generalized formula above | \((k-p)/\operatorname{tr}(P)\), using `model$V` and `model.matrix(model)` | Sums `model$sigma2` | Yes |
| Current `brms` tutorial | No I² calculation; it is explicitly deferred | None currently | None currently | None currently |
| Removed pre-checkpoint `brms` code at commit `8ad126b` | Intended draw-wise \(I^2\) ratios | Simple harmonic mean, `1 / mean(1 / vi)` | Intended `var_total_hetero` plus component numerators | Intended, but not runnable: `var_total_hetero`, component draws, and `summ_I2()` were undefined |

Important details:

- The `orchaRd` ratio quantity is neither the arithmetic mean nor the simple harmonic mean. For an intercept-only model with diagonal \(V\), it is algebraically identical to the generalized \((k-p)/\operatorname{tr}(P)\) formula. The audit verified this equality to machine precision for both Grau-Andrés and Spain.
- With moderators, `orchaRd`'s default ratio method still uses \(k-1\), whereas the generalized method uses \(k-p\) and the full design matrix. The matrix definition is therefore the coherent default for a tutorial meant to cover meta-regression.
- `orchaRd` 2.2.1 stops when any `model$tau2 > 0`. Spatial `rma.mv` variances from `SPEXP` are stored in `tau2`, so `i2_ml()` cannot directly compute the requested spatial I². A transparent local calculation must explicitly include the spatial `tau2` as a marginal variance component.
- The current phylogenetic `orchaRd` result is possible because its phylogenetic variance is represented as a `sigma2` random-intercept component with a supplied correlation matrix. The default ratio denominator is valid for that intercept-only, diagonal-\(V\) fit, but its biological interpretation in the tutorial is too strong.
- The manuscript's 97.3% total and 13.3/10.0/38.6/35.5% components reproduce the current tutorial's `orchaRd` output after rounding, yet the manuscript calls its denominator a mean sampling variance. The generating call instead uses the Higgins-Thompson typical variance. The equation's prose and the implemented calculation therefore disagree.
- The manuscript's statement that all three packages give an “identical” I² is not reproducible from the current tutorial: the current `brms` section has no calculation, and the removed code used a different (simple harmonic-mean) denominator. Similar rounded values would not establish a common definition.
- No visible current spatial tutorial code calculates I². The manuscript's older spatial values and “mean sampling variance” wording must not be reused for the corrected spatial fits.

The simple harmonic-mean alternative has appeared in the literature, but the [`metafor` methodological note](https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate) identifies it as an alternative rather than the standard generalized implementation. The generalized definition also has the decisive advantage of accommodating moderators and non-diagonal known sampling covariance under one formula.

## 2. Meaning for correlated phylogenetic and spatial effects

Suppose component \(j\) has

\[
u_j\sim N(0,\tau_j^2R_j),\qquad \operatorname{diag}(R_j)=1.
\]

Then, for observation \(i\),

\[
\operatorname{Var}(u_{ji})=\tau_j^2,
\]

while for two observations \(i\ne l\),

\[
\operatorname{Cov}(u_{ji},u_{jl})=\tau_j^2R_{j,il}.
\]

The scalar \(\tau_j^2\) can therefore enter the marginal variance sum \(H\). Its component I² is the percentage of the typical total marginal variance of an observed effect allocated to that fitted random-effect component. It is **not**:

- the percentage of heterogeneity “explained by geographic distance” or phylogeny;
- a measure of the correlation between two distinct locations/species;
- a measure of correlation range or decay;
- evidence that the structured component is well identified.

For `metafor` spatial models, the documentation states that covariance is \(\tau^2\) times the spatial correlation; for `SPEXP`, \(R_{il}=\exp(-d_{il}/\rho)\). Thus \(d_{ii}=0\) gives a unit diagonal, while off-diagonal covariance depends jointly on \(\tau^2\), \(\rho\), and distance ([`rma.mv` spatial covariance documentation](https://wviechtb.github.io/metafor/reference/rma.mv.html)). A high spatial component I² can coexist with negligible covariance between distinct sites if \(\rho\) is extremely short.

The same marginal interpretation applies to a phylogenetic component only when the supplied matrix is a correlation matrix or is otherwise scaled so its diagonal is one. If a phylogenetic covariance matrix has nonconstant diagonal elements, or if a random-slope design makes \(\operatorname{diag}(ZGZ')\) observation-dependent, summing scalar variance parameters is not enough. One must state an averaging rule for the component diagonals (for example, the mean of `diag(Z G Z')`) or report observation-specific quantities. The current tutorial's Brownian correlation matrix and the exponential spatial kernels have unit diagonals, so this complication does not affect the audited models.

I² also remains conditional on sampling precision. A high I² does not by itself mean that absolute heterogeneity is scientifically large; the same heterogeneity variance yields different I² under different sampling variances.

## 3. Non-diagonal known sampling V

The proposed definition naturally accommodates correlated sampling errors:

1. Use the known sampling variance-covariance matrix \(V\), in exactly the same row order as the response and \(X\).
2. Set \(W=V^{-1}\), not `diag(1 / diag(V))`.
3. Compute \(P=W-WX(X'WX)^{-1}X'W\).
4. Compute \(\widetilde v=(k-p)/\operatorname{tr}(P)\).

This matches the non-independent-sampling-errors extension in the `metafor` guidance. \(V\) must be symmetric positive definite for the ordinary inverse; singular or nearly singular sampling covariance requires a separately justified numerical/statistical treatment, not silent diagonalization. The resulting \(\widetilde v\) is a design- and weight-based representative variance, not the arithmetic mean of `diag(V)`.

The audit code verifies the trace calculation against an explicitly constructed \(P\) for a positive-definite non-diagonal test matrix; the absolute difference was \(5.55\times10^{-17}\).

## 4. Full global Grau-Andrés I²

All three fits use the same 2,361 observations, diagonal sampling \(V\), and intercept-only \(X\). From the actual saved `vi` and `X`:

- arithmetic mean of `vi`: 0.398691308748
- simple harmonic mean of `vi`: 0.111037631577
- recommended generalized \(\widetilde v\): **0.111123796928**
- \(\operatorname{tr}(P)\): 21237.575256081

The values below are plug-in I² estimates from the saved REML variance estimates; no uncertainty interval for the frequentist I² was calculated.

| Model | Total I² (%) | Effect-size I² (%) | Study I² (%) | Spatial I² (%) |
|---|---:|---:|---:|---:|
| Unstructured-only | 94.7786 | 36.4329 | 58.3457 | — |
| Spatial-only | 94.7965 | 37.1529 | — | 57.6437 |
| Combined | 94.7457 | 36.6612 | 55.1584 | 2.9261 |

Interpretation:

- Total I² is essentially unchanged across the three covariance specifications, consistent with similar total fitted marginal heterogeneity.
- In the spatial-only model, 57.64% is allocated to the spatially structured random effect's *marginal variance*. It is not evidence for strong broad-scale spatial autocorrelation. The fitted \(\rho\approx0.050\) km is extremely short and poorly resolved, so the implied covariance between distinct locations falls off almost immediately.
- In the combined model, the point-estimate spatial contribution is 2.93%. This does not justify “no spatial effect”: the earlier profile audit showed that the additional spatial variance and range are weakly identified. I² inherits that uncertainty and must be presented beside the variance/range profiles, not as a definitive partition.

## 5. Spain cross-package I²

All packages use the same 186 Spain observations and known sampling variances. The common intercept-only representative variance is:

- arithmetic mean of `vi`: 0.339165600844
- simple harmonic mean of `vi`: 0.131741099353
- recommended generalized \(\widetilde v\): **0.132532287191**
- \(\operatorname{tr}(P)\): 1395.886269835

### Frequentist matched fits

| Package | Total I² (%) | Effect-size I² (%) | Spatial I² (%) |
|---|---:|---:|---:|
| `metafor` | 86.363533 | 21.988383 | 64.375150 |
| `glmmTMB` | 86.363527 | 21.988388 | 64.375138 |

The variance components agree within \(4.7\times10^{-7}\), and applying the same \(\widetilde v\) gives the same I² to numerical precision. This agreement is expected from this deliberately matched `equalto()`/Gaussian-residual implementation; it should not be generalized to arbitrary `glmmTMB` meta-analytic specifications.

### `brms` posterior, draw-by-draw

The same fixed \(\widetilde v=0.132532287191\) was used in every posterior draw. Each draw used `sigma^2` as iid effect-size heterogeneity and `sdgp^2` as spatial marginal variance before calculating total and component I².

| Quantity | Posterior median (%) | 95% CrI (%) |
|---|---:|---:|
| Total I² | 86.9989 | 78.5307, 93.4688 |
| Effect-size I² | 21.5329 | 9.5424, 42.5766 |
| Spatial I² | 65.2178 | 38.1214, 83.2313 |

At every draw, the two component I² values sum exactly to total I². Their separately reported posterior medians and interval endpoints need not add because quantiles are nonlinear summaries of different posterior distributions. Confidence-based plug-in estimates and Bayesian posterior intervals also have different inferential meanings; numerical similarity is only a cross-package consistency check.

## 6. Reproducibility and files

- Exact code: [`i2_definition_audit.R`](i2_definition_audit.R)
- Full Grau results: [`i2_definition_audit_outputs/grau_i2.csv`](i2_definition_audit_outputs/grau_i2.csv)
- Full Grau sampling denominator: [`i2_definition_audit_outputs/grau_sampling_variance_summary.csv`](i2_definition_audit_outputs/grau_sampling_variance_summary.csv)
- Spain `metafor`/`glmmTMB`: [`i2_definition_audit_outputs/spain_i2_metafor_glmmTMB.csv`](i2_definition_audit_outputs/spain_i2_metafor_glmmTMB.csv)
- Spain `brms`: [`i2_definition_audit_outputs/spain_i2_brms_posterior.csv`](i2_definition_audit_outputs/spain_i2_brms_posterior.csv)
- Algebra checks: [`i2_definition_audit_outputs/formula_validation_checks.csv`](i2_definition_audit_outputs/formula_validation_checks.csv)
- Package record: [`i2_definition_audit_outputs/session_info.txt`](i2_definition_audit_outputs/session_info.txt)

The script reads saved objects only and asserts common data/design structures, the Grau cross-model `vi`/`X` match, Spain `metafor`/`glmmTMB` variance-component agreement, draw-wise Bayesian additivity, equivalence of the ratio and generalized formulas for the two intercept-only diagonal-\(V\) analyses, and agreement of direct versus trace-shortcut calculations for a non-diagonal \(V\).

## Sources checked

- Viechtbauer, W. [`I² for Multilevel and Multivariate Models`](https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate). This provides the standard typical variance, generalized \(P\)-matrix definition, multilevel component partition, and non-diagonal-\(V\) extension.
- Viechtbauer, W. [`rma.mv` documentation](https://wviechtb.github.io/metafor/reference/rma.mv.html). This defines the spatial covariance as spatial correlation times \(\tau^2\), documents `SPEXP`, and distinguishes known sampling \(V\) from the fitted marginal covariance.
- Locally installed `orchaRd` 2.2.1 source for `i2_ml()`, `ratio_i2()`, and `matrix_i2()`.
- Current revision source `tutorial_v2.qmd` and its pre-mechanical-correction state at commit `8ad126b`.
- Revision manuscript PDF: `/Users/ayumi/Downloads/Tutorial___Phylo_spatial_meta_analysis_2 (6).pdf`.
