# Gates: cleaned Gaussian kernel sensitivity

OWNS: revision_checks/cleaned_gaussian_kernel_audit.R, revision_checks/validate_cleaned_gaussian_kernel.R, revision_checks/cleaned_gaussian_kernel_audit_2026-09-10.md, revision_checks/cleaned_gaussian_kernel_outputs/**, revision_checks/GATES_cleaned_gaussian_kernel.md, revision_checks/SPATIAL_AUDIT_HANDOFF_2026-09-09.md, revision_checks/validate_spatial_audit_handoff.R

Scope: Fit only the cleaned-data Gaussian spatial-only and targeted combined models, validate their inputs and outputs, and add their source-of-truth record to the existing spatial handoff without editing tutorial, manuscript, or response sources.

- [x] G1: The saved `published_cleaned` object is used directly and its 2,355 effect sizes, 390 studies, 380 sites, and spatial matrix order are verified before fitting.
  CHECK: Rscript revision_checks/validate_cleaned_gaussian_kernel.R inputs
  EXPECT: CLEANED_GAUSSIAN_INPUTS_VALIDATED
  CWD: ..
  EVIDENCE: `Rscript revision_checks/validate_cleaned_gaussian_kernel.R inputs` exited 0 and printed `CLEANED_GAUSSIAN_INPUTS_VALIDATED` on 2026-09-10.

- [x] G2: The Gaussian spatial-only and the three specified combined starts are saved with the requested common model hierarchy and finite fitted summaries.
  CHECK: Rscript revision_checks/validate_cleaned_gaussian_kernel.R fits
  EXPECT: CLEANED_GAUSSIAN_FITS_VALIDATED
  CWD: ..
  EVIDENCE: `Rscript revision_checks/validate_cleaned_gaussian_kernel.R fits` exited 0 and printed `CLEANED_GAUSSIAN_FITS_VALIDATED` on 2026-09-10.

- [x] G3: The limited combined multi-start and zero-spatial-variance restriction quantify Gaussian spatial-component identification without a broad profile/grid search.
  CHECK: Rscript revision_checks/validate_cleaned_gaussian_kernel.R identification
  EXPECT: CLEANED_GAUSSIAN_IDENTIFICATION_VALIDATED
  CWD: ..
  EVIDENCE: `Rscript revision_checks/validate_cleaned_gaussian_kernel.R identification` exited 0 and printed `CLEANED_GAUSSIAN_IDENTIFICATION_VALIDATED` on 2026-09-10.

- [x] G4: Any reported Gaussian generalized I2 uses the validated cleaned-data v_tilde and has additive component partitions.
  CHECK: Rscript revision_checks/validate_cleaned_gaussian_kernel.R i2
  EXPECT: CLEANED_GAUSSIAN_I2_VALIDATED
  CWD: ..
  EVIDENCE: `Rscript revision_checks/validate_cleaned_gaussian_kernel.R i2` exited 0 and printed `CLEANED_GAUSSIAN_I2_VALIDATED` on 2026-09-10.

- [x] G5: The updated handoff cites all new source artifacts and all newly modified paths remain audit-only.
  CHECK: Rscript revision_checks/validate_cleaned_gaussian_kernel.R handoff
  EXPECT: CLEANED_GAUSSIAN_HANDOFF_VALIDATED
  CWD: ..
  EVIDENCE: `Rscript revision_checks/validate_spatial_audit_handoff.R` and `Rscript revision_checks/validate_cleaned_gaussian_kernel.R handoff` each exited 0 and printed their success tokens on 2026-09-10.
