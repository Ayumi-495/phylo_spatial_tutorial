# Tutorial RDS inventory

This inventory covers every `.rds` filename referred to by the current `tutorial_v2.qmd`, including one path-only reference. RDS objects remain local and Git-ignored. Generation scripts and this inventory are repository source artifacts to be committed later, but no commit has been made in this pass.

| Filename | Intended contents | Status | Current location | Generating script | Required by current executing render? | Expected revised replacement? |
|---|---|---|---|---|---|---|
| `phylo_eg1_meta_BM.rds` | Moura example BM `metafor::rma.mv` meta-analysis without moderators; an optional precomputed version of the displayed BM fit | Missing | None found in local worktrees or Git history | No standalone script found; model is displayed in `tutorial_v2.qmd` | No, the loader is commented | No planned replacement; decide later whether the optional loader should be retained |
| `moura2021_BM_meta_reg.rds` | Moura BM `metafor::rma.mv` meta-regression with `temporally.pooled`, study/effect/species/phylogenetic random effects, and BM `A` matrix | Available | `Rdata/tutorial_v2/` in both original main and revision worktrees | `revision_checks/moura_bm_meta_reg_restore.R` (permanent precomputation and validation recipe) | Yes | No |
| `moura2021_BM_brms.rds` | Moura BM `brms` meta-analysis used for posterior visualisation | Missing | None found locally | Existing fitting code is in `tutorial_v2.qmd`; no separate script yet | No, chunk is `eval: false` | Yes, scheduled for replacement during R9; do not restore the old fit |
| `phylo_eg2_metafor_BM_mr.rds` | Lim et al. (2014) BM `metafor::rma.mv` meta-regression with `environment` moderator | Missing | None found locally | Existing fitting code is in `tutorial_v2.qmd` | No, its current reference is an `eval: false` path-only assignment | Undecided; needed only if this metafor visualisation remains in the final tutorial |
| `phylo_eg2_brms_mr.rds` | Lim et al. (2014) BM `brms` meta-regression with `environment`, phylogenetic, species, and sampling-variance covariance terms | Missing | None found locally | Existing fitting code is in `tutorial_v2.qmd` | No, chunk is `eval: false` | Undecided; needed only if this brms meta-regression visualisation remains in the final tutorial |

`phylo_eg2_metafor_BM_mr.rds` merits later code review: the current `eval: false` chunk assigns its path with `here()` rather than loading it with `readRDS()`. This pass makes no change to that dormant code.
