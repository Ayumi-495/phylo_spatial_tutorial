# Gates: finalize spatial numerical audits before document revision

OWNS: revision_checks/spatial_audit_GATES.md, revision_checks/validate_spatial_audit_final.R, revision_checks/grau_global_gaussian_*, revision_checks/gaussian_global_outputs/**, revision_checks/scholer_spatial_*, revision_checks/scholer_spatial_audit_outputs/**, revision_checks/spatial_audit_final_synthesis_*.md

Scope: settle the Grau-Andres Gaussian and Scholer exponential model results, profiles, identifiability conclusions, and requested I2 values without modifying tutorial or manuscript sources.

- [x] G0: this ledger states executable outcome checks that can fail
  CHECK: node /Users/ayumi/.codex/skills/unlazy/scripts/gate-lint.mjs spatial_audit_GATES.md
  EXPECT: LINT OK
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=LINT OK

- [x] G1: limited Grau Gaussian free refits distinguish a common optimum, local optima, or a flat optimizer-dependent ridge and preserve the zero-spatial-variance comparison
  CHECK: Rscript validate_spatial_audit_final.R grau_gaussian
  EXPECT: GRAU_GAUSSIAN_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=GRAU_GAUSSIAN_VALIDATED

- [x] G2: corrected Scholer unstructured-only, spatial-only, and combined primary fits are saved with the audited data and grouping structure
  CHECK: Rscript validate_spatial_audit_final.R scholer_primary
  EXPECT: SCHOLER_PRIMARY_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=SCHOLER_PRIMARY_VALIDATED

- [x] G3: Scholer spatial-only and combined variance and rho profiles are saved and support an evidence-qualified identifiability assessment
  CHECK: Rscript validate_spatial_audit_final.R scholer_profiles
  EXPECT: SCHOLER_PROFILES_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=SCHOLER_PROFILES_VALIDATED

- [x] G4: generalized Scholer I2 uses the actual V and X and its component values add to total I2 for every retained model
  CHECK: Rscript validate_spatial_audit_final.R scholer_i2
  EXPECT: SCHOLER_I2_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=SCHOLER_I2_VALIDATED

- [x] G5: retained Grau Gaussian I2 uses the verified common v_tilde and matches the finalized variance components
  CHECK: Rscript validate_spatial_audit_final.R grau_gaussian_i2
  EXPECT: GRAU_GAUSSIAN_I2_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=GRAU_GAUSSIAN_I2_VALIDATED

- [x] G6: the final synthesis compares Grau exponential, Grau Gaussian, and Scholer exponential without overstating spatial variance or range identification
  CHECK: Rscript validate_spatial_audit_final.R synthesis
  EXPECT: SPATIAL_SYNTHESIS_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=SPATIAL_SYNTHESIS_VALIDATED

- [x] G7: tutorial_v2.qmd remains identical to spatial-completion commit 0e57293 and no manuscript or response-letter file is changed
  CHECK: Rscript validate_spatial_audit_final.R protected_sources
  EXPECT: PROTECTED_SOURCES_UNCHANGED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=PROTECTED_SOURCES_UNCHANGED
