#!/usr/bin/env Rscript

# Fit and validate the retained phylogenetic brms examples used in tutorial_v2.qmd.
# Each expensive fit is saved locally as a Git-ignored RDS only after HMC checks.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("moura_mr", "lim_vcv_ma", "lim_se_ma", "lim_mr")) {
  stop("Usage: Rscript brms_consistency_precompute.R {moura_mr|lim_vcv_ma|lim_se_ma|lim_mr}", call. = FALSE)
}
model_id <- args[[1L]]
root <- normalizePath(".")

suppressPackageStartupMessages({
  library(ape)
  library(brms)
  library(ggplot2)
  library(metadat)
  library(metafor)
  library(posterior)
  library(rstan)
  library(tidybayes)
})

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

configurations <- list(
  moura_mr = list(
    rds = "moura2021_BM_meta_reg_brms.rds", seed = 20260912L,
    label = "Moura et al. BM meta-regression",
    ppc_group = "temporally.pooled"
  ),
  lim_vcv_ma = list(
    rds = "lim2014_BM_brms_vcv.rds", seed = 20260913L,
    label = "Lim et al. BM meta-analysis (VCV formulation)",
    ppc_group = NULL
  ),
  lim_se_ma = list(
    rds = "lim2014_BM_brms_se.rds", seed = 20260914L,
    label = "Lim et al. BM meta-analysis (se formulation)",
    ppc_group = NULL
  ),
  lim_mr = list(
    rds = "phylo_eg2_brms_mr.rds", seed = 20260915L,
    label = "Lim et al. BM meta-regression",
    ppc_group = "environment"
  )
)
config <- configurations[[model_id]]
# The Lim VCV, corrected se(), and meta-regression formulations are rerun
# unchanged after divergent transitions; only their integrator targets are more
# conservative.
adapt_delta <- if (model_id %in% c("lim_vcv_ma", "lim_se_ma", "lim_mr")) 0.99 else 0.95
fit_config <- list(chains = 4L, cores = 4L, iter = 6000L, warmup = 2000L,
                   adapt_delta = adapt_delta, max_treedepth = 15L)
# Preserve the failed 0.95 candidates and their diagnostics alongside the
# sampler-only 0.99 reruns, rather than overwriting either audit record.
run_id <- if (identical(model_id, "lim_vcv_ma")) {
  "lim_vcv_ma_adapt_delta_0.99"
} else if (identical(model_id, "lim_se_ma")) {
  "lim_se_ma_adapt_delta_0.99"
} else if (identical(model_id, "lim_mr")) {
  "lim_mr_adapt_delta_0.99"
} else {
  model_id
}
out_dir <- file.path(root, "revision_checks", "brms_consistency_outputs", run_id)
rds_path <- file.path(root, "Rdata", "tutorial_v2", config$rds)
candidate_path <- file.path(out_dir, sub("\\.rds$", "_candidate.rds", config$rds))

prepare_moura <- function() {
  dat <- dat.moura2021$dat
  dat$effect.size.id <- factor(seq_len(nrow(dat)))
  dat$species.id.phy <- dat$species.id
  dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  tree <- compute.brlen(dat.moura2021$tree)
  A <- vcv(tree, corr = TRUE)
  tips <- rownames(A)
  dat$species.id.phy <- factor(as.character(dat$species.id.phy), levels = tips)
  dat$species.id <- factor(dat$species.id)
  dat$study.id <- factor(dat$study.id)
  dat$temporally.pooled <- factor(dat$temporally.pooled)
  assert(!anyNA(dat$species.id.phy), "Moura phylogenetic species are not aligned to A.")
  list(dat = dat, A = A)
}

prepare_lim <- function() {
  dat <- dat.lim2014$o_o_unadj
  tree <- compute.brlen(dat.lim2014$o_o_unadj_tree)
  dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  A <- vcv(tree, corr = TRUE)
  tips <- rownames(A)
  dat$phy <- factor(as.character(dat$species), levels = tips)
  dat$species <- factor(dat$species)
  dat$id <- factor(seq_len(nrow(dat)))
  dat$environment <- factor(dat$environment)
  assert(!anyNA(dat$phy), "Lim phylogenetic species are not aligned to A.")
  V <- diag(dat$vi)
  rownames(V) <- colnames(V) <- levels(dat$id)
  list(dat = dat, A = A, V = V)
}

priors_for <- function(groups, has_slope, fixed_vcv) {
  priors <- c(
    set_prior("normal(0, 1)", class = "Intercept"),
    set_prior("exponential(1)", class = "sigma")
  )
  for (group in groups) priors <- c(priors, set_prior("exponential(1)", class = "sd", group = group))
  if (has_slope) priors <- c(priors, set_prior("normal(0, 1)", class = "b"))
  if (fixed_vcv) priors <- c(priors, set_prior("constant(1)", class = "sd", group = "id"))
  priors
}

make_specification <- function(id) {
  if (identical(id, "moura_mr")) {
    x <- prepare_moura()
    formula <- bf(yi | se(sqrt(vi), sigma = TRUE) ~ 1 + temporally.pooled +
                    (1 | study.id) + (1 | species.id) +
                    (1 | gr(species.id.phy, cov = A)))
    return(list(
      dat = x$dat, data2 = list(A = x$A), formula = formula,
      priors = priors_for(c("study.id", "species.id", "species.id.phy"), TRUE, FALSE),
      parameters = list(
        c("b_Intercept", "Pooled mean (Fisher's Z)", "Fixed effect"),
        c("b_temporally.pooledyes", "Temporally pooled contrast", "Fixed effect"),
        c("sigma", "IID effect-size variance", "Variance"),
        c("sd_study.id__Intercept", "Study variance", "Variance"),
        c("sd_species.id__Intercept", "Non-phylogenetic species variance", "Variance"),
        c("sd_species.id.phy__Intercept", "Phylogenetic species variance", "Variance")
      ),
      metafor = function() rma.mv(yi, vi, mods = ~ temporally.pooled,
        random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
        R = list(species.id.phy = x$A), data = x$dat, method = "REML", sparse = TRUE),
      metafor_values = function(fit) c(
        "Pooled mean (Fisher's Z)" = fit$b[[1L]],
        "Temporally pooled contrast" = fit$b[[2L]],
        "IID effect-size variance" = fit$sigma2[[2L]],
        "Study variance" = fit$sigma2[[1L]],
        "Non-phylogenetic species variance" = fit$sigma2[[3L]],
        "Phylogenetic species variance" = fit$sigma2[[4L]]
      ),
      group = "temporally.pooled"
    ))
  }

  x <- prepare_lim()
  if (identical(id, "lim_vcv_ma")) {
    formula <- bf(yi ~ 1 + (1 | species) + (1 | gr(phy, cov = A)) + (1 | gr(id, cov = V)))
    return(list(
      dat = x$dat, data2 = list(A = x$A, V = x$V), formula = formula,
      priors = priors_for(c("species", "phy"), FALSE, TRUE),
      parameters = list(c("b_Intercept", "Pooled mean (Fisher's Z)", "Fixed effect"),
                        c("sigma", "IID effect-size variance", "Variance"),
                        c("sd_species__Intercept", "Non-phylogenetic species variance", "Variance"),
                        c("sd_phy__Intercept", "Phylogenetic species variance", "Variance")),
      metafor = function() rma.mv(yi, vi, random = list(~ 1 | id, ~ 1 | phy, ~ 1 | species),
        R = list(phy = x$A), data = x$dat, method = "REML", sparse = TRUE), group = NULL
      , metafor_values = function(fit) c(
        "Pooled mean (Fisher's Z)" = fit$b[[1L]], "IID effect-size variance" = fit$sigma2[[1L]],
        "Non-phylogenetic species variance" = fit$sigma2[[3L]], "Phylogenetic species variance" = fit$sigma2[[2L]]
      )
    ))
  }
  if (identical(id, "lim_se_ma")) {
    # se(sqrt(vi), sigma = TRUE) supplies V and represents iid heterogeneity once as sigma.
    formula <- bf(yi | se(sqrt(vi), sigma = TRUE) ~ 1 + (1 | species) + (1 | gr(phy, cov = A)))
    return(list(
      dat = x$dat, data2 = list(A = x$A), formula = formula,
      priors = priors_for(c("species", "phy"), FALSE, FALSE),
      parameters = list(c("b_Intercept", "Pooled mean (Fisher's Z)", "Fixed effect"),
                        c("sigma", "IID effect-size variance", "Variance"),
                        c("sd_species__Intercept", "Non-phylogenetic species variance", "Variance"),
                        c("sd_phy__Intercept", "Phylogenetic species variance", "Variance")),
      metafor = function() rma.mv(yi, vi, random = list(~ 1 | id, ~ 1 | phy, ~ 1 | species),
        R = list(phy = x$A), data = x$dat, method = "REML", sparse = TRUE), group = NULL
      , metafor_values = function(fit) c(
        "Pooled mean (Fisher's Z)" = fit$b[[1L]], "IID effect-size variance" = fit$sigma2[[1L]],
        "Non-phylogenetic species variance" = fit$sigma2[[3L]], "Phylogenetic species variance" = fit$sigma2[[2L]]
      )
    ))
  }
  formula <- bf(yi ~ 1 + environment + (1 | species) + (1 | gr(phy, cov = A)) + (1 | gr(id, cov = V)))
  list(
    dat = x$dat, data2 = list(A = x$A, V = x$V), formula = formula,
    priors = priors_for(c("species", "phy"), TRUE, TRUE),
    parameters = list(c("b_Intercept", "Intercept (Fisher's Z)", "Fixed effect"),
                      c("b_environmentwild", "Wild-environment contrast", "Fixed effect"),
                      c("sigma", "IID effect-size variance", "Variance"),
                      c("sd_species__Intercept", "Non-phylogenetic species variance", "Variance"),
                      c("sd_phy__Intercept", "Phylogenetic species variance", "Variance")),
    metafor = function() rma.mv(yi, vi, mods = ~ environment,
      random = list(~ 1 | id, ~ 1 | phy, ~ 1 | species),
      R = list(phy = x$A), data = x$dat, method = "REML", sparse = TRUE),
    metafor_values = function(fit) c(
      "Intercept (Fisher's Z)" = fit$b[[1L]], "Wild-environment contrast" = fit$b[[2L]],
      "IID effect-size variance" = fit$sigma2[[1L]],
      "Non-phylogenetic species variance" = fit$sigma2[[3L]], "Phylogenetic species variance" = fit$sigma2[[2L]]
    ), group = "environment"
  )
}

summarise_draws <- function(fit, specs) {
  draws_array <- posterior::as_draws_array(rstan::extract(fit$fit, permuted = FALSE, inc_warmup = FALSE))
  all <- as.data.frame(posterior::summarise_draws(draws_array, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(all)[names(all) == "variable"] <- "parameter"
  names(all)[names(all) == "posterior::rhat"] <- "rhat"
  names(all)[names(all) == "posterior::ess_bulk"] <- "ess_bulk"
  names(all)[names(all) == "posterior::ess_tail"] <- "ess_tail"
  key_names <- vapply(specs, `[[`, character(1), 1L)
  key <- all[match(key_names, all$parameter), ]
  assert(!anyNA(key$parameter), "Expected posterior parameter is absent from the fitted model.")

  sampler <- as.data.frame(brms::nuts_params(fit))
  energies <- split(sampler$Value[sampler$Parameter == "energy__"], sampler$Chain[sampler$Parameter == "energy__"])
  bfmi <- vapply(energies, function(x) mean(diff(x)^2) / stats::var(x), numeric(1))
  depth <- sampler$Value[sampler$Parameter == "treedepth__"]
  hmc <- data.frame(
    max_rhat = max(all$rhat[is.finite(all$rhat)]),
    min_bulk_ess = min(all$ess_bulk[is.finite(all$ess_bulk)]),
    min_tail_ess = min(all$ess_tail[is.finite(all$ess_tail)]),
    divergences = sum(sampler$Parameter == "divergent__" & sampler$Value > 0),
    max_treedepth_seen = max(depth), max_treedepth = fit_config$max_treedepth,
    min_bfmi = min(bfmi)
  )
  list(draws_array = draws_array, all = all, key = key, hmc = hmc,
       bfmi = data.frame(chain = as.integer(names(bfmi)), bfmi = unname(bfmi)))
}

write_parameter_plot <- function(draws_array, specs) {
  draw_df <- posterior::as_draws_df(draws_array)
  pieces <- lapply(specs, function(s) {
    values <- draw_df[[s[[1L]]]]
    if (identical(s[[3L]], "Variance")) values <- values^2
    data.frame(value = values, parameter = s[[2L]], panel = s[[3L]])
  })
  plot_data <- do.call(rbind, pieces)
  plot_data$parameter <- factor(plot_data$parameter, levels = rev(unique(plot_data$parameter)))
  p <- ggplot(plot_data, aes(x = value, y = parameter)) +
    tidybayes::stat_halfeye(.width = c(0.5, 0.95), point_interval = "median_qi",
                             fill = "#73A9AD", color = "#24535A") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
    facet_grid(panel ~ ., scales = "free_x", space = "free_y") +
    labs(x = "Posterior distribution (median, 50% and 95% credible intervals)", y = NULL) +
    theme_classic(base_size = 12)
  ggsave(file.path(out_dir, "posterior_parameter_distributions.png"), p, width = 10, height = 7, dpi = 180)
  write.csv(plot_data, file.path(out_dir, "posterior_parameter_draws.csv"), row.names = FALSE)
}

write_ppc <- function(fit, dat, group, seed) {
  set.seed(seed + 100L)
  yrep <- brms::posterior_predict(fit, ndraws = 500L)
  stats <- rbind(
    data.frame(statistic = "mean", observed = mean(dat$yi), replicated_median = median(rowMeans(yrep)),
      lower_95 = unname(quantile(rowMeans(yrep), .025)), upper_95 = unname(quantile(rowMeans(yrep), .975)),
      p_replicated_ge_observed = mean(rowMeans(yrep) >= mean(dat$yi))),
    data.frame(statistic = "standard_deviation", observed = sd(dat$yi), replicated_median = median(apply(yrep, 1L, sd)),
      lower_95 = unname(quantile(apply(yrep, 1L, sd), .025)), upper_95 = unname(quantile(apply(yrep, 1L, sd), .975)),
      p_replicated_ge_observed = mean(apply(yrep, 1L, sd) >= sd(dat$yi)))
  )
  if (!is.null(group)) {
    g <- dat[[group]]
    lev <- levels(factor(g))
    assert(length(lev) == 2L, "Moderator PPC requires two factor levels.")
    index1 <- which(g == lev[[1L]]); index2 <- which(g == lev[[2L]])
    contrasts <- rowMeans(yrep[, index2, drop = FALSE]) - rowMeans(yrep[, index1, drop = FALSE])
    observed <- mean(dat$yi[index2]) - mean(dat$yi[index1])
    stats <- rbind(stats, data.frame(statistic = paste0("mean_contrast_", lev[[2L]], "_minus_", lev[[1L]]),
      observed = observed, replicated_median = median(contrasts), lower_95 = unname(quantile(contrasts, .025)),
      upper_95 = unname(quantile(contrasts, .975)), p_replicated_ge_observed = mean(contrasts >= observed)))
  }
  write.csv(stats, file.path(out_dir, "posterior_predictive_summary.csv"), row.names = FALSE)
  png(file.path(out_dir, "posterior_predictive_density_overlay.png"), width = 1800, height = 1200, res = 180)
  print(brms::pp_check(fit, type = "dens_overlay", ndraws = 100L))
  dev.off()
}

write_comparison <- function(fit, spec) {
  frequentist <- spec$metafor()
  draws <- posterior::as_draws_df(posterior::as_draws_array(rstan::extract(fit$fit, permuted = FALSE, inc_warmup = FALSE)))
  values <- vapply(spec$parameters, function(s) {
    x <- draws[[s[[1L]]]]
    if (identical(s[[3L]], "Variance")) x <- x^2
    median(x)
  }, numeric(1))
  rows <- data.frame(parameter = vapply(spec$parameters, `[[`, character(1), 2L),
                     brms_posterior_median = values)
  rows$metafor_estimate <- unname(spec$metafor_values(frequentist)[rows$parameter])
  rows$difference <- rows$brms_posterior_median - rows$metafor_estimate
  write.csv(rows, file.path(out_dir, "metafor_comparison.csv"), row.names = FALSE)
}

diagnostics_pass <- function(x) {
  h <- x$hmc
  all(x$key$rhat <= 1.01) && all(x$key$ess_bulk >= 400) && all(x$key$ess_tail >= 400) &&
    h$divergences == 0L && h$max_treedepth_seen < h$max_treedepth && h$min_bfmi >= .3
}

assert(!file.exists(rds_path), paste("Refusing to overwrite existing RDS:", rds_path))
assert(!file.exists(candidate_path), paste("Candidate fit already exists:", candidate_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
spec <- make_specification(model_id)
brms::validate_prior(spec$priors, spec$formula, data = spec$dat, data2 = spec$data2, family = gaussian())

fit <- brm(formula = spec$formula, family = gaussian(), data = spec$dat, data2 = spec$data2,
  prior = spec$priors, backend = "cmdstanr", chains = fit_config$chains, cores = fit_config$cores,
  iter = fit_config$iter, warmup = fit_config$warmup, seed = config$seed,
  control = list(adapt_delta = fit_config$adapt_delta, max_treedepth = fit_config$max_treedepth),
  save_pars = save_pars(all = TRUE), refresh = 100L)
saveRDS(fit, candidate_path)

diagnostics <- summarise_draws(fit, spec$parameters)
write.csv(diagnostics$all, file.path(out_dir, "all_parameter_diagnostics.csv"), row.names = FALSE)
write.csv(diagnostics$key, file.path(out_dir, "key_parameter_diagnostics.csv"), row.names = FALSE)
write.csv(diagnostics$hmc, file.path(out_dir, "hmc_diagnostics.csv"), row.names = FALSE)
write.csv(diagnostics$bfmi, file.path(out_dir, "bfmi_by_chain.csv"), row.names = FALSE)
write_parameter_plot(diagnostics$draws_array, spec$parameters)
write_ppc(fit, spec$dat, spec$group, config$seed)
write_comparison(fit, spec)
writeLines(c(
  config$label, paste("Formula:", deparse(spec$formula)),
  "Priors: normal(0, 1) for intercept and moderators; exponential(1) for estimated SDs; fixed VCV sampling term constant(1) where used.",
  paste("chains", fit_config$chains), paste("iter", fit_config$iter), paste("warmup", fit_config$warmup),
  paste("seed", config$seed), paste("adapt_delta", fit_config$adapt_delta), paste("max_treedepth", fit_config$max_treedepth)
), file.path(out_dir, "manifest.txt"))

assert(diagnostics_pass(diagnostics), "Diagnostics failed: outputs saved but candidate RDS was not promoted.")
dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
assert(file.rename(candidate_path, rds_path), "Could not promote validated RDS.")
cat("BRMS_CONSISTENCY_FIT_VALIDATED:", model_id, "\n")
