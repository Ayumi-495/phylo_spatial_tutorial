# Gates: Reviewer 18 influential-effect audit

OWNS: revision_checks/reviewer18_influential_effects.R, revision_checks/reviewer18_compile_comparison.R, revision_checks/reviewer18_validate_outputs.R, revision_checks/reviewer18_influential_effects_2026-09-09.md, revision_checks/reviewer18_influential_effects_outputs/**

Scope: Recover the source-publication influential-effect exclusions, fit the three corrected spatial models to the published-cleaned spatial data, and compare them with the verified full-data fits without editing tutorial or manuscript sources.

- [x] G1: the archived code selections and Dryad rows are reconciled without inventing a replacement for a missing record
  CHECK: Rscript revision_checks/reviewer18_validate_outputs.R provenance
  EXPECT: PROVENANCE_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision; path=8b3a4328aec5/19 entries; output=PROVENANCE_VALIDATED

- [x] G2: the all-spatially-usable and published-cleaned datasets have verified row, study, site, and exclusion counts
  CHECK: Rscript revision_checks/reviewer18_validate_outputs.R datasets
  EXPECT: DATASETS_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision; path=8b3a4328aec5/19 entries; output=DATASETS_VALIDATED

- [x] G3: all three published-cleaned metafor fits completed with the required common data, design, sampling variance, and geodesic site matrix
  CHECK: Rscript revision_checks/reviewer18_validate_outputs.R fits
  EXPECT: FITS_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision; path=8b3a4328aec5/19 entries; output=FITS_VALIDATED

- [x] G4: the saved comparison and report answer the pooled-effect, model-ranking, variance-allocation, spatial-interpretation, and Spain-impact questions
  CHECK: Rscript revision_checks/reviewer18_validate_outputs.R report
  EXPECT: REPORT_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision; path=8b3a4328aec5/19 entries; output=REPORT_VALIDATED

- [x] G5: tutorial_v2.qmd has no tracked modifications in the revision worktree
  CHECK: Rscript revision_checks/reviewer18_validate_outputs.R source
  EXPECT: SOURCE_UNCHANGED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision; path=8b3a4328aec5/19 entries; output=SOURCE_UNCHANGED
