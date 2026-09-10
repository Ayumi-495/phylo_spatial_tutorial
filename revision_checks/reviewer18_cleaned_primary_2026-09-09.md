# Cleaned primary Grau-Andres spatial audit

This companion audit is numerical support only. No QMD, manuscript, or response-letter file was edited.

## Evidence

- published_cleaned retains 2,355 effect sizes, 390 studies, and 380 recorded-coordinate locations after the seven source-publication influential-effect exclusions. The full sensitivity dataset has 2,361, 393, and 383, respectively.
- The cleaned generalized-I2 calculation used the cleaned-data diagonal sampling-variance matrix V and intercept-only X. v_tilde = 0.110937235978 (tr(P) = 21,219.2054296).

| model | component variances | generalized I2 (%) |
|---|---:|---:|
| unstructured-only | effect 0.7516; study 1.1383; total 1.8899 | effect 37.5630; study 56.8924; total 94.4554 |
| spatial-only | effect 0.7630; spatial 1.1495; total 1.9125 | effect 37.7092; spatial 56.8082; total 94.5174 |
| combined | effect 0.7514; study 1.0943; spatial 0.04785; total 1.8935 | effect 37.4865; study 54.5919; spatial 2.3870; total 94.4654 |

- Cleaned primary fits: unstructured-only mean -0.3639 (95% CI -0.4845 to -0.2434), AIC 7944.375; spatial-only mean -0.3558 (-0.4785 to -0.2332), AIC 7968.495; combined mean -0.3624 (-0.5446 to -0.1802), AIC 7946.509. The combined fit retained the matrix/S4 deprecation warning seen in the full audit; no fit reported an optimizer error.
- The targeted profiles contain 43 fixed-parameter fits (9 spatial-only tau2, 11 spatial-only rho, 10 combined tau2, 13 combined rho), saved incrementally on Totoro under revision_checks/reviewer18_cleaned_primary_outputs/profiles/.

## Interpretation

- Spatial-only spatial variance is clearly separated from zero: fixing tau2 = 0 loses 381.047 REML log-likelihood units. Its profile maximum is near rho = 0.05 km, whereas the free fit is rho = 0.1369 km; the difference is 0.0033 log-likelihood units. The range is therefore very short and poorly resolved (the targeted 95%-likelihood region spans 0.005 to 0.5 km).
- In the combined model, tau2 = 0 is only 0.933 REML log-likelihood units below the maximum, so the additional spatial variance is weakly identified. The tau2 profile is broad and ridge-like (0.01 to 0.1 within 0.5 log-likelihood units; 0 to 0.4 within the targeted 95%-likelihood threshold).
- The combined rho profile has its grid maximum near 385 km, but the free estimate is 2,058.9 km and differs by less than 0.001 log-likelihood units at the fitted point. Values from roughly 200 to 12,000 km are within the targeted 95%-likelihood region; rho identifiability is weak and rho is not a well-defined correlation range.
- The limited multi-start check at fixed spatial tau2 = 0.8, 1.0, and 1.3 found no meaningful competing optimum. At tau2 = 1.3, the best alternate rho was 0.151 km and improved the single-start log-likelihood by only 0.0196 units.
- Relative to the full 2,361-effect sensitivity analysis, the AIC ranking is unchanged: unstructured-only (rank 1), combined (rank 2), spatial-only (rank 3). Cleaned means shift by -0.0135 to -0.0218, while total generalized I2 remains about 94.5%. Study I2 decreases modestly and the combined spatial I2 decreases from about 2.93% to 2.39%; the qualitative identifiability conclusions are unchanged.

## Uncertainty

- These are targeted fixed-parameter likelihood profiles, not formal confidence intervals. The cleaned spatial-only rho is effectively a sub-kilometre boundary-scale estimate, and the combined variance/range surface has a broad ridge. The point estimate rho about 2,059 km should not be interpreted as evidence for a resolved range.
- Profile fits carry the same non-fatal Matrix/S4 deprecation warning in some combined evaluations. No warning indicated a failed fit.

## Recommendation

Use published_cleaned as the candidate primary worked example and retain the full 2,361-effect analysis as sensitivity analysis. The pooled mean is robust across covariance specifications, whereas variance allocation and spatial interpretation are model-sensitive. The existing Gaussian results can remain labelled explicitly as a full-data sensitivity; a cleaned-data Gaussian refit would be preferable for a direct kernel comparison beside the cleaned primary, but was not run in this audit.

Outputs: reviewer18_cleaned_primary_outputs/cleaned_generalized_i2.csv, cleaned_profile_results_compiled.csv, cleaned_profile_summary.csv, targeted_multistart_compiled.csv, cleaned_vs_full_generalized_i2.csv, and cleaned_vs_full_aic_ranks.csv.
