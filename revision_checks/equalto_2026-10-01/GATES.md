# Gates: current tutorial equalto alignment

Scope: Align all five current tutorial equalto models, verify using glmmTMB 1.1.15.2, preserve prior evidence, and render the revised tutorial. Removed global Grau/Scholer models and the manuscript are outside this edit.

- [x] G1: All five model blocks construct value-aligned named sampling matrices and aligned phylogenetic matrices.
  CHECK: Rscript revision_checks/equalto_2026-10-01/verify_alignment.R
  EXPECT: FIVE_BLOCKS_ALIGNMENT_PERMUTATION_AND_DUPLICATE_CHECKS_PASSED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial; path=b9a6a80b7e09/23 entries; output=FIVE_BLOCKS_ALIGNMENT_PERMUTATION_AND_DUPLICATE_CHECKS_PASSED
- [x] G2: All five models run under glmmTMB 1.1.15.2, with old/new comparisons and convergence evidence recorded.
  CHECK: Rscript revision_checks/equalto_2026-10-01/compare.R
  EXPECT: FIVE_MODELS_BASELINE_COMPARISON_PASSED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial; path=b9a6a80b7e09/23 entries; output=Warning message: | package ‘glmmTMB’ was built under R version 4.6.1
- [x] G3: Revised index.html renders, source and displayed results agree, and independent review passes.
  CHECK: python3 revision_checks/equalto_2026-10-01/verify_html.py
  EXPECT: FINAL_HTML_VERSION_ALIGNMENT_AND_LOCAL_IMAGE_REFERENCES_PASSED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial; path=b9a6a80b7e09/23 entries; output=FINAL_HTML_VERSION_ALIGNMENT_AND_LOCAL_IMAGE_REFERENCES_PASSED
- [x] G4: Current state and specialist execution are saved in a project report and AYUMI continuation notes.
  EVIDENCE: Read-back of AYUMI HANDOFF.md, worklog.md and project hub.md verified 2026-10-01 and equalto task record; REPORT.md retains actual Turing/Wald/Nightingale handles and completed final render.
