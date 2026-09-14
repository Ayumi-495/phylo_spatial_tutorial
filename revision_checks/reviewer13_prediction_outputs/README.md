# Reviewer 13 prediction audit

Target: the latent true effect for a new effect-size observation from a new study and a marginal new species-level realization. It includes fitted study, effect-size, non-phylogenetic species, and phylogenetic species heterogeneity, plus uncertainty in the pooled mean. It excludes future sampling error.

The frequentist interval is a plug-in t prediction interval conditional on estimated variance components. It does not separately propagate uncertainty in those component estimates.

The legacy Figure 3 thin line came from orchaRd::orchard_plot() through pred_interval_esmeans(), not a direct metafor call. For this simple BM model it agrees with metafor::predict.rma() with default newvi = 0.

The Bayesian distribution was constructed explicitly from joint posterior draws and four new marginal random-effect draws. brms::posterior_predict() was not used because it simulates observed outcomes with sampling error.
