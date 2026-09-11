# OU correctness audit: Moura example

## Scope and reproducibility

This audit concerns only the phylogenetic exponential/OU section of `tutorial_v2.qmd`. It does not alter the spatial analysis, R13 prediction target, Figure 3 architecture, Bayesian example, manuscript, or response letter.

The analysis used the primary `metadat::dat.moura2021` data (1,828 effect sizes, 457 studies, 341 species) and the tree supplied with that data. The tree was transformed with `ape::compute.brlen()` at its default Grafen-type setting. It is ultrametric after this operation, has height 1 in constructed Grafen branch-length units, and does not contain estimated divergence times.

The primary audit script is `revision_checks/ou_correctness_audit.R`. The 101-point fixed-rho profile was fitted independently on totoro using `revision_checks/ou_rho_profile_parallel.R` with at most 100 workers, then copied back as machine-readable output. `revision_checks/summarise_ou_profile.R` derives the likelihood-ratio interval and sensitivity subset without refitting.

## Tree, taxa, and distance checks

The BM matrix is `A = ape::vcv(tree, corr = TRUE)`. The exponential model uses the direct matrix `D_phylo = ape::cophenetic.phylo(tree)`, reordered to the same 341-tip order as the factor supplied to `rma.mv()`.

| Check | Result |
|---|---:|
| Tree ultrametric | yes |
| Tree height | 1.0000000000 constructed Grafen units |
| Missing data taxa / missing tree taxa | 0 / 0 |
| Duplicated tips | 0 |
| Grouping-factor order matches matrix order | yes |
| Direct distance symmetry error | 0 |
| Direct distance diagonal maximum | 0 |
| Missing direct-distance entries | 0 |

The old executable expression was `J - A`, where `J` was an all-ones matrix despite being named `I`; it was not literal identity minus `A`. Literal `I - A` has negative off-diagonal entries (minimum -0.9970588235) and is not a distance matrix. For this particular ultrametric tree,

$$J-A = d/(2h) = d/2,$$

where `d` is direct patristic distance and `h = 1`. The maximum numerical error in this relation was `6.66133814775094e-16`. Thus the old implementation happened to be equivalent after a factor-of-two rescaling of rho, but the tutorial now uses direct patristic distance. The jointly estimated raw-distance rho divided by twice the old normalized-distance rho was `1.00000004488074`; the absolute log-likelihood difference was `1.37690676638158e-08`.

## Exponential parameterisation and external validation

For `metafor::rma.mv(..., struct = "SPEXP")`, the verified correlation parameterisation is

$$\operatorname{cor}(d)=\exp(-d/\rho).$$

With the same direct patristic distance and branch-length scale, the equivalent form is `exp(-alpha * d)` with `alpha = 1/rho`. Here, rho is in constructed Grafen branch-length units and alpha is per constructed Grafen branch-length unit.

At the joint estimate, the direct matrix `exp(-D_phylo/rho)` and `ape::corMartins(value = alpha, ...)` had identical row and column names and ordering. Their maximum absolute difference was `5.55111512312578e-17` and mean absolute difference was `1.77439143310408e-19`; they agree to numerical precision.

## Matched BM and exponential-model comparison

All values below are REML results from otherwise matched models with `test = "t"`. Model A jointly estimates rho and is the valid model for the BM AIC comparison. Model B conditions on Model A's estimated matrix and is not used for that comparison because its nominal AIC excludes the rho parameter.

| Model | Mean [95% CI] | REML log-likelihood | Parameters | AIC | Delta AIC vs BM |
|---|---|---:|---:|---:|---:|
| BM | 0.368165851 [0.113113658, 0.623218044] | -167.672633555 | 5 | 345.345267110 | 0 |
| Exponential Model A, direct distance | 0.351407472 [0.281854510, 0.420960434] | -160.378331969 | 6 | 332.756663937 | -12.588603173 |
| Exponential, old normalized distance, diagnostic only | 0.351407472 [0.281854514, 0.420960430] | -160.378331955 | 6 | 332.756663910 | -12.588603201 |
| Fixed-matrix Model B, not an AIC competitor | 0.351407472 [0.281854510, 0.420960435] | -160.378331952 | 5 | 330.756663905 | -14.588603206 |

The Model A estimate was `rho = 0.03640768450` and `alpha = 27.46672890` per constructed Grafen branch-length unit. Its implied correlations at patristic distances 0, 0.1764705882, 0.7882352941, and 2 were 1, 0.007851350, `3.957473e-10`, and `1.389001e-24`, respectively. The corresponding BM correlations were 1, 0.911764706, 0.605882353, and approximately 0.

| Variance component | BM | Exponential Model A |
|---|---:|---:|
| Study | 0.01915845087 | 0.01573329919 |
| Effect size | 0.01445013805 | 0.01441935300 |
| Non-phylogenetic species | 0.05566175428 | 0.000000006425 |
| Phylogenetic species | 0.05122365948 | 0.10293874680 |

## Rho profile and identification

The fixed-rho profile re-estimated all remaining variance components at 101 log-spaced rho values from 0.02 to 50 times the joint estimate. The profile had an interior maximum at `rho = 0.03640768034` (`alpha = 27.46673204`) and declined on both sides. The likelihood-ratio cutoff was 1.920729410.

The grid-supported 95% range was rho 0.02276756624--0.06295769365. Log-rho interpolation at the likelihood-ratio crossings gave an approximate 95% interval of rho 0.02112830523--0.06782727397, corresponding to alpha 14.74333172--47.32987285 per constructed Grafen unit. This is a grid-interpolated likelihood interval, not a claim of exact continuous-profile endpoints.

The profile is comparatively concentrated rather than flat or boundary-limited on the evaluated grid. It supports estimating a covariance-decay parameter for this model and scale, but does not support a biological interpretation of alpha as the strength of stabilising selection, particularly because the branch lengths are Grafen-type constructed values rather than divergence times.

## Why the confidence interval changes

The pooled mean changed modestly (0.3682 under BM versus 0.3514 under Model A), whereas its standard error changed from 0.130044852 to 0.035463348. This follows from the different fitted species covariance assumptions.

| Quantity | BM | Exponential Model A |
|---|---:|---:|
| Mean off-diagonal correlation | 0.332202103 | 0.009318163 |
| Median off-diagonal correlation | 0 | approximately 0 |
| Mean off-diagonal covariance | 0.017016607 | 0.000959200 |
| Median off-diagonal covariance | 0 | approximately 0 |
| Pooled-mean SE | 0.130044852 | 0.035463348 |

The exponential fit places more variance in the structured species component but, because its fitted correlations decline very quickly on the constructed branch-length scale, treats most species pairs as nearly uncorrelated. BM retains appreciable covariance for more pairs. Consequently, the effective information for the overall mean differs sharply even though the central mean is similar. This is a statistical consequence of covariance specification, not evidence for a particular biological process.

Within the likelihood-supported rho values, the pooled mean ranged from 0.3446016 to 0.3540187 and the pooled-mean SE from 0.0293415 to 0.0460537. Thus the increased precision is itself sensitive to plausible rho values, although all supported profile points gave substantially smaller SEs than the BM model.

## Reviewer-response support

- **R5:** The tutorial now constructs direct patristic distances, verifies taxa/order/matrix invariants, documents the scale conversion that made the old code accidentally equivalent for this tree, states `SPEXP` as `exp(-d/rho)`, and validates the direct exponential matrix against `ape::corMartins()`.
- **R6:** AIC is described as relative fit among candidate covariance models. The fixed-matrix AIC is explicitly excluded from the BM comparison, and variance allocation is described as covariance-specification-sensitive rather than biological evidence.
- **R17:** The tutorial and this audit quantify why similar pooled means have very different confidence intervals, compare fitted covariance structures, and show the profile-supported sensitivity of the pooled-mean SE.
