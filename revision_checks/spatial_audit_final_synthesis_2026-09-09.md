# Final numerical synthesis: spatial audit

## What is now settled

All results below come from corrected saved model objects, not historical QMD
output. Grau-Andres uses 2,361 effects, 393 studies, 383 recorded-coordinate
locations, and WGS84 geodesic distances. Scholer uses 949 effects, 205 studies,
454 recorded-coordinate locations, and the same distance convention. In both,
effect, study, and recorded-location identifiers are distinct grouping
structures. The pooled mean is stable within each dataset across the principal
covariance specifications, whereas variance allocation and inference about
spatial structure are model-sensitive.

## Grau-Andres: exponential versus Gaussian

The exponential unstructured-only model has the lowest AIC (8042.979), followed
by exponential combined (8045.604) and spatial-only (8076.949). The combined
exponential fit has a small estimated spatial variance (0.0619; spatial I²
2.93%), but the profile shows the component is also weakly identified: setting
it to zero costs only 0.687 log-likelihood units, and approximately 200-3,000 km
has nearly indistinguishable likelihood.

The Gaussian results lead to the same substantive conclusion. Its spatial-only
model assigns substantial variance to space but collapses to a 0.362-km range,
which does not indicate broad-scale autocorrelation. In the combined model,
targeted refits identify a better approximately 312-km solution and a distinct
approximately 3,091-km stationary solution only 0.384 log-likelihood units
lower. At the best observed solution, spatial variance is 0.0830 and spatial I²
is 3.94%; fixing it to zero costs only 1.063 log-likelihood units. Therefore the
additional component is weakly identified and its rho is optimizer-dependent.
The better observed Gaussian combined AIC is 8044.854, still above the common
unstructured-only model.

Thus the negative pooled mean is robust to exponential versus Gaussian decay,
but neither combined kernel yields a well-determined additional spatial range
once study heterogeneity is included.

## Scholer: exponential hierarchy

The corrected unstructured-only fit gives mean 0.669 [0.558, 0.779] and AIC
1656.479. The combined fit gives 0.657 [0.536, 0.778] and AIC 1657.090. The
restricted spatial-only fit is much worse (AIC 1864.607) and shifts the mean to
0.458 because it lacks study heterogeneity.

Within the spatial-only restriction, spatial variance (0.3227) and rho
(168.71 km) have clear profile maxima. This does not show that the combined
hierarchy needs a strong spatial component. In the combined fit, spatial
variance is small (0.0204; spatial I² 2.94%) and weakly identified: setting it
to zero costs 1.694 log-likelihood units. Its rho profile is extremely broad;
every checked value from 10 to 12,000 km lies within 1.92 log-likelihood units
of the maximum. The primary solution agrees with the profile-grid maximum, so
there is no evidence requiring a targeted Scholer multi-start search.

Scholer's generalized representative sampling variance is
`v_tilde = 0.000725865052242`. Total I² is approximately 99.90% in all three
models. In the combined model, effect, study, and spatial I² are 32.43%, 64.53%,
and 2.94%, respectively. These describe marginal variance allocation, not
variance explained by distance or pairwise spatial correlation.

## Interpretation boundary

A small estimated spatial variance and a weakly identified spatial component
are different statements. A point estimate can be numerically small while a
profile still permits zero and a broad range of alternatives; that is what
happens in both combined analyses. Conversely, a restricted spatial-only model
can estimate a large spatial variance because it is forced to represent
otherwise unmodelled study heterogeneity. Range interpretation must therefore
be tied to the fitted variance, profile likelihood, and model structure—not to
the rho point estimate alone.
