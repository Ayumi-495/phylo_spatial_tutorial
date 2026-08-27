# REM revision plan

Target: revised submission by 22 November 2026.

This file converts the decision letter into auditable changes. Work is developed on `codex/rem-revision-audit`; `main` remains unchanged until the revised analysis and text have been reviewed.

## Order of work

1. Validate the statistical parameterisations.
2. Refit the minimum model comparisons requested by the reviewers.
3. Revise conceptual explanations and interpretations.
4. Add diagnostics and reproducibility safeguards.
5. Prepare the reporting checklist and point-by-point response.
6. Render the complete tutorial and verify all reported values against fresh output.

## Revision matrix

| Work package | Decision-letter concern | Required change | Evidence of completion | Status |
| --- | --- | --- | --- | --- |
| Sampling covariance | Sampling covariance was conflated with effect-level and study-level random effects | Define the sampling-error matrix `V`, explain when off-diagonal terms are required, and state that random effects do not reconstruct sampling covariance | Conceptual section in `tutorial_v2.qmd`; worked non-diagonal `V` example or justified sensitivity analysis | In progress |
| Confounding | Correlated random effects were described as controlling confounding | Replace causal/confounding language with a statement about covariance and uncertainty; state the need for measured moderators or a causal design | Revised interpretation throughout | In progress |
| Phylogenetic OU | Verify the OU construction and avoid overinterpreting AIC | Check the distance transformation against patristic distance and `ape::corMartins`; use the direct model AIC that counts the fitted range parameter; distinguish covariance fit from evolutionary-process evidence | `analysis/revision_audit_ou.R`, CI log, revised OU section | In progress |
| Branch lengths | Grafen branch lengths can be inappropriate for evolutionary-rate interpretation | State that the example uses constructed branch lengths and that decay parameters are not in evolutionary-time units; recommend dated or otherwise justified branch lengths for biological interpretation | Revised tree and OU text | In progress |
| Phylogenetic uncertainty | One tree was treated as known | Add guidance and a compact workflow for repeating analyses across a posterior or bootstrap sample of trees and pooling summaries | New tutorial subsection | Not started |
| Spatial distance | Global Web Mercator distances and range interpretation were problematic | Use geodesic distances in kilometres for the global example; quantify projection distortion; define `SPEXP` as `exp(-d/rho)` | `analysis/revision_audit_spatial.R`, CI log, revised spatial section | In progress |
| Spatial identifiability | The spatial-only model imposed a strong assumption | Compare unstructured-only, spatial-only, and spatial-plus-unstructured models with study structure retained; report range uncertainty and boundary behaviour | Three-model comparison for full and sensitivity datasets | In progress |
| Extreme effect sizes | Spatial variance appeared implausibly large | Audit the Hedges' d distribution, trace extreme values to source calculations, and report a transparent sensitivity analysis | Data-audit table plus source-level verification notes | In progress |
| Cross-package agreement | Similarity across packages was asserted without aligned models | Align sampling error, unstructured heterogeneity, spatial/phylogenetic covariance, residual variance, and distance units before comparing estimates | Parameter-mapping table and matched-model results | Not started |
| Bayesian practice | Two chains, weak prior reporting, and incomplete diagnostics | Use four chains, report priors and seed, examine R-hat and effective sample sizes, and add posterior predictive checks | Updated code and diagnostic summary | Not started |
| Model diagnostics | Diagnostics were too limited | Add convergence checks, profile or interval checks for range parameters, residual/influence checks where supported, and warnings for weak identification | Diagnostic subsection and saved outputs | Not started |
| Heterogeneity | I-squared formulation and interpretation needed clarification | Define denominator and each variance component; avoid interpreting I-squared as a biological variance fraction without qualification | Revised equations and reporting text | Not started |
| Prediction | Prediction intervals were missing | Report prediction intervals for pooled examples and explain their target population | Updated outputs and text | Not started |
| Decision framework | Readers need to know when phylogenetic or spatial models are warranted | Add an evidence-based decision table covering design, replication, available coordinates/tree, expected dependence, and sensitivity analyses | New decision table | Not started |
| Reproducibility | Seeds, package versions, matrix matching, and logical constants were inconsistent | Add seeds, four-chain defaults, explicit `REML = TRUE`, matrix-name assertions, package versions, and a clean render workflow | CI render and reproducibility section | Not started |
| Reporting and workflow | The tutorial lacked a compact synthesis and reporting checklist | Add workflow and reporting checklists covering effect-size construction, `V`, random effects, covariance matrices, diagnostics, uncertainty, and sensitivity analyses | Final tutorial checklists | Not started |
| Terminology and presentation | Terminology, CI/CrI labels, likelihood assumptions, moderators, bias, and captions need tightening | Apply the decision-letter line edits after the statistical revisions are fixed | Full-text pass and rendered HTML review | Not started |

## Analysis acceptance criteria

### OU audit

- `1 - A` is shown to equal patristic distance divided by twice the ultrametric tree height for this tree.
- The fitted exponential matrix is numerically equal to `ape::corMartins` after converting distance units.
- The model comparison uses the AIC of the model that jointly estimates the range parameter.
- No sentence treats an AIC difference as proof of a biological OU process.

### Spatial audit

- Distances are geodesic and reported in kilometres.
- Repeated observations at one coordinate share a location identifier.
- The three requested covariance structures are fitted with the same fixed effects, sampling variances, and study structure.
- `rho` is reported as a distance scale with an interval and a clear weak-identification warning when applicable.
- Extreme Hedges' d values are enumerated, checked at source level, and handled only through a labelled sensitivity analysis.

### Sampling covariance

- The tutorial distinguishes `V` from random-effect covariance.
- A diagonal `V` is explicitly described as an assumption of independent sampling errors.
- Shared controls, repeated measurements, and multiple outcomes are identified as possible sources of off-diagonal sampling covariance.
- Any assumed correlation is justified and varied in sensitivity analyses.

## Final deliverables

- Revised `tutorial_v2.qmd` and rendered `index.html`.
- Reproducible audit scripts and machine-readable model summaries.
- A concise reporting checklist in the tutorial.
- A point-by-point response letter in which every response identifies the changed analysis or exact section.
- A clean pull request for review; no automatic merge into `main`.
