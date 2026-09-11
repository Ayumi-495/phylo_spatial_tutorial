# Gates: OU correctness revision 2026-09-11

OWNS: tutorial_v2.qmd, revision_checks/OU_CORRECTNESS_AUDIT_2026-09-11.md, revision_checks/ou_correctness_audit.R, revision_checks/ou_rho_profile_parallel.R, revision_checks/summarise_ou_profile.R, revision_checks/ou_correctness_outputs/**, revision_checks/GATES_ou_correctness_2026-09-11.md

Scope: directly construct and validate the Moura phylogenetic exponential correlation, document its identification limits, and revise only the corresponding OU tutorial section.

- [x] G1: numerical OU audit reproduces the Moura distance, fit, and profile checks from source data
  CHECK: Rscript revision_checks/ou_correctness_audit.R --check
  EXPECT: OU_CORRECTNESS_AUDIT_CHECKS_PASSED
  EVIDENCE: `OU_CORRECTNESS_AUDIT_CHECKS_PASSED` on 2026-09-11.

- [x] G2: tutorial OU text uses direct patristic distance and does not make prohibited biological claims
  CHECK: Rscript revision_checks/ou_correctness_audit.R --check-tutorial
  EXPECT: OU_TUTORIAL_CHECKS_PASSED
  EVIDENCE: `OU_TUTORIAL_CHECKS_PASSED` on 2026-09-11.

- [x] G3: the tutorial renders to HTML without executing chunks
  CHECK: quarto render tutorial_v2.qmd --to html --no-execute
  EXPECT: Output created: index.html
  EVIDENCE: `quarto render tutorial_v2.qmd --to html --no-execute` completed on 2026-09-11.

- [x] G4: an executing render reaches only the pre-existing missing Moura BM meta-regression RDS boundary
  CHECK: Rscript revision_checks/ou_correctness_audit.R --check-executing-render
  EXPECT: EXPECTED_EXECUTING_RENDER_BLOCK_CONFIRMED
  EVIDENCE: `EXPECTED_EXECUTING_RENDER_BLOCK_CONFIRMED` on 2026-09-11; chunk 49 failed only at `moura2021_BM_meta_reg.rds`.

- [x] G5: diff scope is limited to the OU audit artifacts and focused tutorial OU section
  EVIDENCE: `git diff --check` passed; generated HTML/PNG validation artifacts were removed before review.
