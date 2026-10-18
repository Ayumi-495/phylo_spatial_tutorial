# Figure scripts

Run these independent scripts from the repository root with `Rscript scripts/<filename>.R`. They recreate figures from stored inputs and do not fit or profile models. They can overwrite PNGs under `figs/tutorial/`.

| Script | Required inputs |
| --- | --- |
| [replot_brms_parameter_figures.R](replot_brms_parameter_figures.R) | Six posterior-draw CSVs under `data/figure_inputs/`; writes Bayesian parameter figures. |
| [create_visualisation_interval_figures.R](create_visualisation_interval_figures.R) | Four interval CSVs and the Spain Bayesian summary CSV; writes package result figures and the Spain comparison. No fitted models are loaded. |
| [create_phylogenetic_orchard_figures.R](create_phylogenetic_orchard_figures.R) | Moura/Lim interval CSVs plus three fitted-model RDS files; writes phylogenetic orchard figures. |
| [create_spatial_result_figures.R](create_spatial_result_figures.R) | Six Grau-Andrés/Scholer fits and fourteen profile CSV/RDS pairs; writes their spatial result figures. |

For example:

```sh
Rscript scripts/replot_brms_parameter_figures.R
```

## Stored tables

The following paths are relative to `data/figure_inputs/`. Only the tables consumed by these scripts are included.

Posterior draws:

- `moura_bm/moura_bm_brms_parameter_draws.csv`
- `moura_mr/posterior_parameter_draws.csv`
- `lim_vcv/posterior_parameter_draws.csv`
- `lim_se/posterior_parameter_draws.csv`
- `lim_mr/posterior_parameter_draws.csv`
- `spain/brms_parameter_draws.csv`

Intervals and posterior summaries:

- `intervals/moura_frequentist_intervals.csv`
- `intervals/lim_frequentist_intervals.csv`
- `intervals/spain_metafor_intervals.csv`
- `intervals/spain_glmmTMB_wald_intervals.csv`
- `spain/brms_spatial_only_result.csv`

For each of `grau_andres/profiles/` and `scholer/profiles/`, the spatial script reads seven CSVs:

```text
unstructured_only_sigma2_1.csv
unstructured_only_sigma2_2.csv
spatial_only_sigma2_1.csv
spatial_only_tau2_1.csv
combined_sigma2_1.csv
combined_sigma2_2.csv
combined_tau2_1.csv
```

## Saved model and profile objects

RDS inputs belong under `Rdata/figure_inputs/`, separately from the tutorial's own `Rdata/tutorial_v2/` objects. RDS files are excluded from Git and are available as `reproducibility-models-a1f2ee1.zip` in the [Zenodo archive](https://doi.org/10.5281/zenodo.23454135). Extract the ZIP into the repository root. GitHub alone therefore does not provide every input needed by the two model-based figure scripts.

The phylogenetic orchard script requires these three files under `Rdata/figure_inputs/intervals/`:

```text
moura_metafor_bm.rds
moura_glmmTMB_bm.rds
lim_metafor_bm_meta_regression.rds
```

The spatial script requires `unstructured_only.rds`, `spatial_only.rds` and `combined.rds` under each of `Rdata/figure_inputs/grau_andres/` and `Rdata/figure_inputs/scholer/`. It also needs the fourteen profile objects under the corresponding `profiles/` directories, with the same basenames as the CSVs and extension `.rds`. These objects contain interval-boundary signs used in the figures.

The spatial model/profile RDS files are included in that archive bundle. Their CSV summaries do not preserve the interval-boundary signs. All four scripts have been run successfully with the stored inputs, producing 21 PNG files identical to the supplied versions. A full tutorial render and the archived source/model bundles have also been verified.
