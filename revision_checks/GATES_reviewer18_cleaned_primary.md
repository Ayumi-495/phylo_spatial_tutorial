# Gates: Reviewer 18 cleaned-primary numerical support

OWNS: revision_checks/reviewer18_cleaned_i2.R, revision_checks/reviewer18_cleaned_profile_point.R, revision_checks/reviewer18_cleaned_profile_launcher.sh, revision_checks/reviewer18_cleaned_targeted_multistart.R, revision_checks/reviewer18_cleaned_multistart_launcher.sh, revision_checks/reviewer18_cleaned_collect_profiles.R, revision_checks/reviewer18_cleaned_validate.R, revision_checks/reviewer18_cleaned_primary_2026-09-09.md, revision_checks/reviewer18_cleaned_primary_outputs/**

Scope: Add generalized I2 and targeted cleaned-data SPEXP profiles without refitting the completed primary models or editing tutorial/manuscript sources.

- [x] G1: cleaned generalized v_tilde and component I2 use the saved models' actual vi and X and component I2 adds to total I2
  CHECK: cd .. && Rscript revision_checks/reviewer18_cleaned_validate.R i2
  EXPECT: CLEANED_I2_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=CLEANED_I2_VALIDATED

- [x] G2: all declared targeted profile points completed and the fixed primary parameter points reproduce the free-fit maxima within tolerance
  CHECK: cd .. && Rscript revision_checks/reviewer18_cleaned_validate.R profiles
  EXPECT: CLEANED_PROFILES_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=CLEANED_PROFILES_VALIDATED

- [x] G3: cleaned versus full conclusions and the Gaussian recommendation are supported by saved machine-readable comparisons
  CHECK: cd .. && Rscript revision_checks/reviewer18_cleaned_validate.R comparison
  EXPECT: CLEANED_COMPARISON_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=CLEANED_COMPARISON_VALIDATED

- [x] G4: no tutorial, manuscript, response-letter, or cleaned Gaussian analysis was modified or run
  CHECK: cd .. && Rscript revision_checks/reviewer18_cleaned_validate.R scope
  EXPECT: CLEANED_SCOPE_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=CLEANED_SCOPE_VALIDATED

- [x] G5: the companion report contains every requested result and labels evidence, interpretation, uncertainty, and recommendation
  CHECK: cd .. && Rscript revision_checks/reviewer18_cleaned_validate.R report
  EXPECT: CLEANED_REPORT_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=CLEANED_REPORT_VALIDATED

- [x] G6: the authoritative spatial handoff contains every required section and every named local source-of-truth path exists
  CHECK: cd .. && Rscript revision_checks/validate_spatial_audit_handoff.R
  EXPECT: SPATIAL_AUDIT_HANDOFF_PATHS_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=8b3a4328aec5/19 entries; output=SPATIAL_AUDIT_HANDOFF_PATHS_VALIDATED
