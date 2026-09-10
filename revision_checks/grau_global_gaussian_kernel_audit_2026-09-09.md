# Grau-Andres global Gaussian-kernel audit

## Scope and geometry

This audit retains the verified 2,361-effect, 393-study, 383-recorded-location
Grau-Andres dataset and its WGS84 ellipsoidal geodesic distance matrix in
kilometres. The Gaussian `metafor` kernel is exactly
`Cor(d) = exp(-d^2 / rho^2)`, so `rho` is its e-folding distance. No tutorial,
manuscript, or response-letter source was changed.

## Primary fits

| Model | Mean [95% CI] | Effect variance | Study variance | Spatial variance | rho (km) | REML logLik | AIC |
|---|---:|---:|---:|---:|---:|---:|---:|
| SPGAU spatial-only | -0.336 [-0.463, -0.209] | 0.79384 | -- | 1.23720 | 0.362 | -4034.341 | 8076.683 |
| SPGAU combined, original optimizer solution | -0.335 [-0.520, -0.150] | 0.77521 | 1.21322 | 0.03658 | 3090.80 | -4017.811 | 8045.622 |

The spatial-only fit has a large spatial variance but an e-folding distance of
only 0.362 km. As with the exponential spatial-only model, this deliberately
restricted model asks the spatial component to absorb heterogeneity that may not
be spatial; it is not evidence for broad-scale autocorrelation.

## Profile and targeted free-refit result

The original one-dimensional profiles were conditional on the original
approximately 3,091-km optimizer branch. On that branch, the profiled range was
very flat (all evaluated values from about 1,545 to 6,182 km were less than
0.375 log-likelihood units below the branch maximum). A separately optimized
fixed `rho = 200 km` point had logLik -4017.701, which exceeded the original
free solution and triggered the requested limited targeted check.

Only three free starting points were used:

| Start | Final rho (km) | Effect variance | Study variance | Spatial variance | REML logLik | AIC | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| Fixed-profile solution at 200 km | 312.357 | 0.77554 | 1.13952 | 0.08302 | -4017.427 | 8044.854 | completed |
| Intermediate 800 km | 312.355 | 0.77554 | 1.13951 | 0.08302 | -4017.427 | 8044.854 | completed |
| Original 3,091 km solution | 3090.796 | 0.77521 | 1.21322 | 0.03658 | -4017.811 | 8045.622 | completed; matrix-class warning only |

The 200- and 800-km starts converged to the same better solution, whereas the
3,091-km start remained at a distinct, inferior stationary solution. Their
log-likelihood difference is only 0.384 despite the roughly tenfold difference
in `rho`, and the variance allocation also changes. Thus there is evidence of
optimizer-dependent local solutions along a shallow likelihood surface. The
Gaussian `rho` is weakly/practically unidentified; neither 312 km nor 3,091 km
should receive substantive interpretation.

Fixing Gaussian spatial variance to zero gives logLik -4018.489. Relative to
the best observed targeted free solution, the loss is only 1.063
log-likelihood units. The additional Gaussian spatial component is therefore
weakly identified after study heterogeneity is included. This is an
identifiability statement, not merely a claim that its point estimate is small.

No broad multi-start grid was run because the limited targeted refits answered
the specific discrepancy and already established optimizer dependence.

## Generalized I2

The actual sampling-variance matrix and intercept-only design reproduce
`v_tilde = (k-p)/tr(P) = 0.111123796928`.

| Model/solution | Total I2 | Effect I2 | Study I2 | Spatial I2 |
|---|---:|---:|---:|---:|
| SPGAU spatial-only | 94.8125% | 37.0579% | -- | 57.7546% |
| SPGAU combined, best observed 312-km solution | 94.7315% | 36.7695% | 54.0258% | 3.9363% |
| SPGAU combined, 3,091-km stationary solution | 94.7979% | 36.2902% | 56.7952% | 1.7125% |

Spatial I2 is only the share of typical marginal variance allocated to the
fitted spatial component. It is not variance explained by geographic distance,
a range estimate, or pairwise correlation. The solution-dependent component
I2 values reinforce the weak-identification conclusion.

## Kernel-level conclusion

The pooled mean remains negative and similar under exponential and Gaussian
decay. Kernel choice does not alter that biological conclusion. In both
combined models, study-level heterogeneity absorbs most structured variance and
the separately added spatial component is weakly identified. The Gaussian audit
additionally demonstrates optimizer-dependent range solutions, so a single
Gaussian range estimate must not be reported as well determined.
