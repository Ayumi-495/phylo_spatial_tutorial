# Online tutorial equalto alignment audit, 2026-10-01

The current `tutorial_v2.qmd` was corrected and its five `equalto()` model definitions rerun with CRAN `glmmTMB` 1.1.15.2. All five fits have optimiser convergence code 0, positive-definite Hessians, and no fitting warnings. This verifies implementation and numerical convergence, not model adequacy or parameter identification.

## Scope and source protection

Started from clean `main`, commit `b667a594f23b4d8360d7dfd18192616e7e4b761d`. No manuscript or response letter was edited. At the initial audit closure, no commit, push or publication had been requested or performed; Ayumi subsequently authorised committing and pushing this verified checkpoint. Preserved pre-edit source is `tutorial_before.qmd`; the pre-edit rendered HTML is locally saved in `/private/tmp/ayumi-equalto-index-before.html`.

The supplied assessment refers to an older tutorial. Current global Grau-Andrés exponential/Gaussian and Scholer sections use `metafor`, with no `equalto()` blocks. This task does not restore or rerun their superseded glmmTMB models and does not establish whether their historical sampling-variance assignments were wrong. The current Grau glmmTMB example is the Spain regional exponential fit; it was rerun here. Existing global numerical-audit choices and the OU HOLD were preserved.

## Changes

Each of the five current model preparations converts the effect identifier to a factor, names its sampling variances by observation ID, selects values in factor-level order, and applies matching row and column names. Assertions check unique/nonmissing IDs, finite positive variances, both name orders, and recovery of original observation-level variances from the matrix diagonal. `diag(..., nrow = length(effect_levels))` also handles a single observation correctly.

Both Moura and Lim MA/MR blocks explicitly align `A` to the corresponding phylogenetic factor levels. Lim MR rebuilds its sampling matrix instead of relying on the MA block. The introductory claim that all mismatches cause an error was replaced with explicit alignment guidance.

The current old Moura code sorts A alphabetically after the phylogenetic factor has been aligned to tree-tip order. Both old Moura blocks fail the propto order check under the locally installed glmmTMB 1.1.15. This additional issue is fixed. It is separate from sampling-variance assignment.

## Reruns and comparisons

| Model | New pooled intercept | Baseline | Largest absolute compared difference |
|---|---:|---|---:|
| Moura MA | 0.3681658334 | Pre-existing saved glmmTMB fit | 1.70e-11 |
| Moura MR | 0.3561815242 | Preserved source's printed, rounded fixed effects and components | 4.20e-9 |
| Lim MA | -0.1307540727 | Preserved source rerun with glmmTMB 1.1.15 | 0 |
| Lim MR | -0.1395075251 | Preserved source rerun with glmmTMB 1.1.15 | 0 |
| Spain exponential | -0.0972717246 | Preserved source rerun with glmmTMB 1.1.15 | 0 |

The maximum covers the scalar parameters present in `comparison.csv`, including logLik/AIC where available; it is not a relative error. Moura MR has no saved prior glmmTMB RDS, so its comparison is limited to rounded displayed estimates, with no old logLik/AIC comparison. Spain also reproduces its pre-existing regional CSV: iid variance 0.2137041829, spatial variance 0.6256591500, range 28.96280422 km, REML logLik -241.32069984, AIC 490.64139968.

Known sampling matrices are fixed covariance terms. `known_sampling_variance_first` in the result CSV is their first diagonal element, not an estimated heterogeneity component. The iid variance is `sigma(fit)^2`. `pdHess` does not establish identification of every component or range.

Five actual source constructions also passed permutation tests using character IDs `2,10,1` and explicitly ordered factors. They recovered each observation's sampling variance. Duplicate IDs failed as intended.

## Environment and evidence

CRAN lists 1.1.15.2, published 2026-09-29: <https://cran.r-project.org/package=glmmTMB>. The package and dependencies were installed into `/private/tmp/ayumi-equalto-library`; the user's normal library was not changed. Reproduction requires installing that exact version into an appropriate local library or updating the audit's temporary path.

The macOS binary reports it was built under R 4.6.1; the executing R version is in `new_sessionInfo.txt`. Old saved glmmTMB models cannot safely be queried with the new DLL schema: a VarCorr attempt failed on `combinom_disp_Link`. The old Moura baseline was therefore extracted in a separate session with installed 1.1.15; new RDS objects were inspected with 1.1.15.2.

- `audit.R`: executes actual pre-edit/current QMD model blocks, separately by package environment; RDS outputs are git-ignored.
- `new_results.csv`, `old_results.csv`, `saved_moura_baseline.csv`, `comparison.csv`: full scalar evidence and provenance.
- `old_errors.txt`: both old Moura propto order-check failures.
- `verify_alignment.R`, `compare.R`, `extract_saved_baseline.R`: reproducible checks.
- `new_run.log`: all five fits completed; the first run then hit an empty-error-list logging bug. That logging bug was fixed without refitting, and `verify_run.log` verifies all saved new fits.
- `new_sessionInfo.txt`, `old_sessionInfo.txt`: runtime evidence.

## Review and rendering

Turing `/root/turing_equalto` independently confirmed scope, ordering safeguards, and Spain baseline extraction. Wald `/root/wald_equalto` independently re-extracted new RDS with 1.1.15.2 and old RDS with 1.1.15. All five numerical checks passed within the stated baseline limits. The raw Moura MA RDS maximum difference was 1.64e-11; the CSV difference is 1.70e-11 because of decimal serialization. Nightingale `/root/nightingale_equalto` checked source, scalar evidence, final HTML and package version; no code/numerical corrections were required. Its requested correction of pending report statuses was applied.

Initial sandboxed executing render stopped at an existing OpenTree network chunk before the edited examples. Two network-enabled executing renders with the temporary 1.1.15.2 library completed all 179 render steps. The final run is in `render_final.log`; it includes the final source validation note and `sessionInfo()` reports glmmTMB_1.1.15.2. HTML text, all five displayed matrix constructions, and all local image references passed `verify_html.py`. The older `validate_spatial_tutorial_revision.R` did not pass: it expects a heading absent in both pre-edit and current source, then errors because its label argument is missing. This is a pre-existing validator/source mismatch; no success is claimed for that unrelated legacy check. Browser file-URL preview was denied by its URL policy; no browser screenshot or visual-layout check is claimed. The edit changes code and prose, not figure layouts.

## Continuation

Ayumi subsequently authorised the commit/push of this verified checkpoint. A historical full-data Grau equalto audit, if wanted, must use the superseded code/data and retain its separate provenance. It is not needed to restore removed models to this current tutorial.

Final closure: AYUMI project HANDOFF, worklog and hub updated and read back on 2026-10-01. Four gates met, zero unmet or abandoned. All changes were local and uncommitted at the initial audit closure. The subsequent publication checkpoint is recorded in AYUMI continuation notes.
