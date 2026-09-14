# Figure 3 profile-likelihood audit

`reviewer13_figure3_profile_audit.R` profiles one constrained parameter at a
time from the saved final `rma.mv` fit, allowing the remaining model parameters
to be re-optimized by `metafor::confint.rma.mv()`. It does not overwrite or
otherwise alter the saved final fits. Each target has a CSV audit row and an RDS
containing the raw `confint` profile result. `collate` rebuilds the summary CSV
from completed target rows without rerunning profiles.

The four Moura variance intervals are newly computed, full-precision profile-
likelihood intervals from the final BM model. They are not copied from the
rounded `confint()` display in `tutorial_v2.qmd`.

The global `rho` searches are preserved even when a likelihood-ratio crossing
reaches a search boundary. Such results are marked `finite_two_sided = FALSE`
and are not plotted as finite two-sided confidence intervals.

`reviewer13_global_profile_boundary_audit.R` recovers the boundary signs from
the raw global `confint.rma` objects and writes
`global_profile_boundary_audit.csv`. The supplementary figure draws bars only
for profiles with finite two-sided intervals; boundary-limited components and
spatial ranges are shown as point estimates with their identification status.
