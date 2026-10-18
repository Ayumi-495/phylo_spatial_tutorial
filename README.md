# A unified framework for phylogenetic and spatial meta-analysis

Online tutorial accompanying *A unified framework for phylogenetic and spatial meta-analysis: concepts, implementation, and practical guidance*, by Ayumi Mizuno, Coralie Williams, Malgorzata Lagisz, Alistair M. Senior and Shinichi Nakagawa.

[Read the tutorial](https://ayumi-495.github.io/phylo_spatial_tutorial/).

## Contents

| Path | Contents |
| --- | --- |
| `tutorial_v2.qmd` | Tutorial text and analysis code |
| `data/` | Example datasets, trees and stored figure inputs |
| `scripts/` | Four figure-generation scripts and their [input requirements](scripts/README.md) |
| `figs/tutorial/` | Tutorial figures |
| `index.html`, `tutorial_v2_files/` | Website and supporting assets |
| `phylo_spatial.Rproj` | RStudio project |

## Run the tutorial

Read the rendered tutorial without installing R. To try its examples, install the packages in the **Setup** section; `phytools` and `geosphere` are also used later. The CSV and Nexus inputs are read directly from commit-fixed GitHub URLs, and the Moura and Lim examples use the `metadat` package. These data-loading examples do not require cloning the repository.

For local work, download and extract the repository ZIP or clone it, open `phylo_spatial.Rproj`, and run the setup and relevant preparation chunks before their model chunks. Bayesian refitting requires a working Stan toolchain; examples using `cmdstanr` require CmdStan.

Several examples load precomputed model objects from `Rdata/tutorial_v2/`. These large RDS files are excluded from Git and are supplied separately during revision. A GitHub download alone is therefore insufficient for a full executing render. The required filenames are:

```text
moura2021_BM_brms.rds
moura2021_BM_meta_reg.rds
moura2021_BM_meta_reg_brms.rds
lim2014_BM_brms_vcv.rds
lim2014_BM_brms_se.rds
phylo_eg2_brms_mr.rds
spain_spatial_only_brms.rds
```

With the required objects and software installed, render from the repository root:

```sh
quarto render tutorial_v2.qmd
```

Rendering executes enabled chunks that load saved model objects and run an Open Tree of Life example that requires internet access. Chunks marked `eval: false` do not run automatically.

## Figures and reproducibility

The four [figure scripts](scripts/README.md) use retained tables or fitted objects and do not fit or profile models. Their CSV inputs are included under `data/figure_inputs/`. Model-based figures additionally require separately supplied RDS files under `Rdata/figure_inputs/`. All four figure-generation scripts have been run successfully with the stored inputs and reproduce the supplied PNG files byte-for-byte. A full tutorial render has also been verified. The final archive check remains pending.

R and package versions are recorded by `sessionInfo()` in the tutorial. Stochastic fitting code includes explicit seeds. The setup code does not pin package versions; the tutorial's `equalto()` examples were checked with `glmmTMB` 1.1.15.2. Original study references and data sources are given in the tutorial; please cite those sources when reusing their data.

The publication version will be archived through Zenodo after the reproducibility materials are finalised. Its permanent DOI will be added here and in the tutorial once available. The live website and repository may change after publication.
