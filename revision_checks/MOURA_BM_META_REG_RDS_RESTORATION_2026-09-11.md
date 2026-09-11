# Moura BM meta-regression RDS restoration

## Recovery search before regeneration

No copy of `moura2021_BM_meta_reg.rds` was found in either the original main worktree, the `revision-online-tutorial` worktree, or the additional local Codex worktree. Searches included tracked, untracked, and ignored filesystem locations. `Rdata/` was absent from both main and revision worktrees before restoration.

Git history was searched across all refs with path history and object-name queries. No reachable Git object contained the target path. The one unreachable commit found by `git fsck` changed only `tutorial_v2.qmd` and did not contain the target or any `Rdata/tutorial_v2/*.rds` path. No recoverable existing copy was found.

## Object restored

`Rdata/tutorial_v2/moura2021_BM_meta_reg.rds` now contains a single `metafor::rma.mv` object. It is the exact model shown earlier in the tutorial under the in-memory name `phylo_eg1.1_meta_BM`:

```r
rma.mv(
  yi, vi,
  mods = ~ temporally.pooled,
  random = list(~ 1 | study.id, ~ 1 | effect.size.id,
                ~ 1 | species.id, ~ 1 | species.id.phy),
  R = list(species.id.phy = vcv(compute.brlen(dat.moura2021$tree), corr = TRUE)),
  data = dat_moura2021,
  sparse = TRUE,
  method = "REML"
)
```

The RDS loader assigns this same fitted object to `moura2021_BM_meta_reg` only so that the following plotting code can use a descriptive name. The tutorial now states this relationship beside the loader.

The model was fit through `Rscript --vanilla` in a clean R session. The saved object was then checked independently without refitting. All displayed tutorial values matched within `0.0005`:

| Quantity | Observed |
|---|---:|
| Intercept | 0.3561815250 |
| `temporally.pooledyes` coefficient | 0.0395335885 |
| Study variance | 0.0194046609 |
| Effect-size variance | 0.0144882390 |
| Non-phylogenetic species variance | 0.0539505408 |
| Phylogenetic species variance | 0.0519847583 |
| REML log-likelihood | -165.9738169596 |
| AIC | 343.9476339192 |

The machine-readable comparison and object manifest are in `revision_checks/moura_bm_meta_reg_restore_outputs/`.

## Render and RDS inventory

After restoration, the full executing command `quarto render tutorial_v2.qmd --to html` completed successfully through all 171 chunks. No subsequent missing artifact was encountered.

The current tutorial refers to five RDS filenames, four through `readRDS()` and one as a path-only assignment. Availability after this restoration is:

| RDS | Available locally | Executed in current full render |
|---|---|---|
| `phylo_eg1_meta_BM.rds` | no | no, commented example |
| `moura2021_BM_meta_reg.rds` | yes | yes |
| `moura2021_BM_brms.rds` | no | no, `eval: false` |
| `phylo_eg2_metafor_BM_mr.rds` | no | no, `eval: false` path-only assignment |
| `phylo_eg2_brms_mr.rds` | no | no, `eval: false` |

The missing Bayesian objects were not regenerated or modified, in accordance with the instruction not to begin the R9 revision.
