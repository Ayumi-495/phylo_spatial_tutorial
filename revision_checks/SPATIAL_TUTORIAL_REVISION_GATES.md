# Gates: spatial tutorial revision

OWNS: tutorial_v2.qmd, revision_checks/validate_spatial_tutorial_revision.R, revision_checks/SPATIAL_TUTORIAL_REVISION_GATES.md

Scope: revise only the spatial tutorial section from finalized audit outputs, without refitting models or changing manuscript and response-letter sources.

- [x] G0: this ledger has executable outcome checks
  CHECK: node /Users/ayumi/.codex/skills/unlazy/scripts/gate-lint.mjs SPATIAL_TUTORIAL_REVISION_GATES.md
  EXPECT: LINT OK
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=c73383804948/19 entries; output=LINT OK

- [x] G1: the revised spatial QMD matches audited primary, sensitivity, Gaussian, Spain, and Scholer evidence
  CHECK: Rscript validate_spatial_tutorial_revision.R spatial_tutorial
  EXPECT: SPATIAL_TUTORIAL_REVISION_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=c73383804948/19 entries; output=SPATIAL_TUTORIAL_REVISION_VALIDATED

- [x] G2: tutorial and audit artefacts render without a broad hidden spatial wrapper or global projected spatial example
  CHECK: Rscript validate_spatial_tutorial_revision.R spatial_tutorial
  EXPECT: SPATIAL_TUTORIAL_REVISION_VALIDATED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=c73383804948/19 entries; output=SPATIAL_TUTORIAL_REVISION_VALIDATED

- [x] G3: no manuscript or response-letter source changed from checkpoint 3684290
  CHECK: git diff --name-only 3684290 -- | Rscript -e 'p <- readLines(file("stdin")); ext <- tools::file_ext(p); bad <- grepl("manuscript|response", basename(p), ignore.case = TRUE) & ext %in% c("md", "qmd", "tex", "docx", "pdf"); stopifnot(!any(bad)); cat("PROTECTED_WRITING_SOURCES_UNCHANGED\\n")'
  EXPECT: PROTECTED_WRITING_SOURCES_UNCHANGED
  EVIDENCE: exit=0; shell=/bin/sh; cwd=/Users/ayumi/Library/CloudStorage/GoogleDrive-ayumi.mizuno5@gmail.com/My Drive/research/phylo_spatial_tutorial_revision/revision_checks; path=c73383804948/19 entries; output=PROTECTED_WRITING_SOURCES_UNCHANGED
