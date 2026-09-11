# Gates: Moura BM meta-regression RDS restoration 2026-09-11

OWNS: Rdata/tutorial_v2/moura2021_BM_meta_reg.rds, revision_checks/MOURA_BM_META_REG_RDS_RESTORATION_2026-09-11.md, revision_checks/moura_bm_meta_reg_restore.R, revision_checks/moura_bm_meta_reg_restore_outputs/**, revision_checks/GATES_moura_bm_meta_reg_restore_2026-09-11.md

Scope: recover an existing object if present; otherwise regenerate only the BM meta-regression object shown in `tutorial_v2.qmd`, verify it against the displayed results, and identify the next executing-render boundary. Do not revise R9 or alter scientific model specification.

- [x] G1: all local and Git recovery locations for the named RDS have been searched before regeneration
  EVIDENCE: main, revision, Codex worktree, ignored files, CloudStorage repository locations, all Git refs, and the unreachable commit found by `git fsck` had no target copy.

- [x] G2: the saved RDS is either a recovered copy or a clean-session fit that matches the tutorial's BM meta-regression specification and reported diagnostics
  CHECK: Rscript revision_checks/moura_bm_meta_reg_restore.R --check
  EXPECT: MOURA_BM_META_REG_RDS_RESTORE_CHECKS_PASSED
  EVIDENCE: `MOURA_BM_META_REG_RDS_RESTORE_CHECKS_PASSED` after a `Rscript --vanilla` fit and independent saved-object check.

- [x] G3: a full executing render either succeeds or records the next missing precomputed artifact without creating a placeholder
  CHECK: Rscript revision_checks/moura_bm_meta_reg_restore.R --render-check
  EXPECT: MOURA_BM_META_REG_RENDER_CHECK_COMPLETED
  EVIDENCE: `MOURA_BM_META_REG_RENDER_CHECK_COMPLETED: FULL_RENDER_SUCCEEDED` on 2026-09-11.

- [x] G4: all `readRDS()` dependencies named in the current tutorial are inventoried with their local availability
  CHECK: Rscript revision_checks/moura_bm_meta_reg_restore.R --inventory-check
  EXPECT: MOURA_BM_META_REG_RDS_INVENTORY_COMPLETED
  EVIDENCE: `MOURA_BM_META_REG_RDS_INVENTORY_COMPLETED`; one executable BM object is present and four inactive RDS references remain absent.

- [x] G5: the validated ignored RDS has identical copies in the original main and revision worktrees
  CHECK: Rscript --vanilla revision_checks/moura_bm_meta_reg_restore.R --main-copy-check
  EXPECT: MOURA_BM_META_REG_MAIN_COPY_MATCHES
  EVIDENCE: `MOURA_BM_META_REG_MAIN_COPY_MATCHES`; SHA-256 `aacff59d7ede332ceedcb170055fc9e7dfd113fef420fde49fec46aff871eac4` in both locations.
