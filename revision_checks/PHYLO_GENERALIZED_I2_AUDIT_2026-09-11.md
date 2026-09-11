# Moura phylogenetic generalized I2 audit (2026-09-11)

## Scope

This audit recalculates generalized representative sampling variance, total I2,
and component-specific I2 for the current Moura BM model and audited OU Model A.
It reads the saved `rma.mv` fits from
`revision_checks/ou_correctness_outputs/baseline_fit_objects.rds`; neither model
is refitted. The manuscript and tutorial source are not edited in this audit.

## Definition

Equations 11--13 in the revised manuscript define

\[
\widetilde v = \frac{k-p}{\operatorname{tr}(P)},\quad
P = W-WX(X^\top WX)^{-1}X^\top W,\quad W=V^{-1},
\]

\[
I^2_{\rm total}=\frac{\sum_j h_j}{\sum_j h_j+\widetilde v}\times100,
\qquad
I^2_j=\frac{h_j}{\sum_j h_j+\widetilde v}\times100.
\]

Here each `h_j` is the diagonal marginal variance contributed by a fitted
heterogeneity component. Both phylogenetic correlation matrices have diagonal
one, so the fitted BM phylogenetic `sigma2` and OU Model A phylogenetic `tau2`
are their respective phylogenetic `h_j` values.

## Common representative sampling variance

Both saved fits have the same 1,828 effect sizes, same diagonal sampling matrix
`V = diag(vi)`, and the same intercept-only fixed-effect design (`p = 1`). This
was verified directly from `yi`, `vi`, and `X` in the two fitted objects.

| Quantity | Value |
|---|---:|
| `k` | 1828 |
| `p` | 1 |
| `tr(P)` | 469570.852267723 |
| `v_tilde` (BM) | 0.00389078664311631 |
| `v_tilde` (OU Model A) | 0.00389078664311631 |
| Arithmetic mean `vi` (not used as `v_tilde`) | 0.0490059378755240 |
| Harmonic mean `vi` (not used as `v_tilde`) | 0.00347069277154184 |

The directly formed equation-11 trace and its algebraically equivalent
diagonal-`V` shortcut differed by `9.90e-10`, well below the audit tolerance.

## I2 results

| Model | Component | `h_j` | Denominator | I2 (%) |
|---|---|---:|---:|---:|
| BM | Study | 0.0191584508703000 | 0.144384789320162 | 13.2690229770794 |
| BM | Effect size | 0.0144501380509337 | 0.144384789320162 | 10.0080750326765 |
| BM | Species, non-phylogenetic | 0.0556617542807465 | 0.144384789320162 | 38.5509820964041 |
| BM | Species, phylogenetic | 0.0512236594750657 | 0.144384789320162 | 35.4771854544049 |
| BM | Total | 0.140494002677046 | 0.144384789320162 | 97.3052655605648 |
| OU Model A | Study | 0.0157332991878589 | 0.136982192054056 | 11.4856529538163 |
| OU Model A | Effect size | 0.0144193530012275 | 0.136982192054056 | 10.5264434632038 |
| OU Model A | Species, non-phylogenetic | 0.00000000642492680617 | 0.136982192054056 | 0.00000469033726927 |
| OU Model A | Species, phylogenetic | 0.102938746796926 | 0.136982192054056 | 75.1475394380495 |
| OU Model A | Total | 0.133091405410939 | 0.136982192054056 | 97.1596405454068 |

For each model, the component I2 values sum to total I2 to numerical precision.
The common `v_tilde` is a property of the unchanged sampling-error and
fixed-effect specification; the different denominators reflect the two fitted
sets of heterogeneity variances.

## BM phylogenetic heritability

For the BM fit,

\[
H^2_{\rm phylo} = \frac{0.0512236594750657}
{0.0512236594750657 + 0.0556617542807465}
= 0.479239006288454 \; (47.9239006288454\%).
\]

This is the proportion of fitted **among-species** variance allocated to the
phylogenetically structured component, conditional on the constructed
Grafen-scaled tree and BM covariance specification. It is not variance
"explained by phylogeny."

## Reviewer Comment 12 consistency after the tutorial revision

The revised phylogenetic tutorial now uses the same equations 11--13 and the
same common `v_tilde` for BM and OU Model A. It reports total I2 as typical
marginal variation not attributed to assumed sampling error under the fitted
model, and labels each component-specific I2 as a fitted marginal-variance
allocation. The BM phylogenetic heritability wording is likewise limited to
the fitted among-species variance allocation. Therefore, this tutorial section
now supports the R12 response without treating I2 or Hphylo2 as variance
explained by phylogeny, pairwise correlation, or a biological mechanism.

## Reproducibility files and checks

- Script: `revision_checks/phylo_generalized_i2_audit.R`
- Numerical outputs: `revision_checks/phylo_generalized_i2_outputs/`
- Validation: `Rscript --vanilla revision_checks/phylo_generalized_i2_audit.R --check`
  returned `PHYLO_GENERALIZED_I2_AUDIT_CHECKS_PASSED`.
