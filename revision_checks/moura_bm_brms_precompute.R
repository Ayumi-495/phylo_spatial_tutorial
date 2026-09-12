#!/usr/bin/env Rscript

# Precompute the Moura et al. BM brms meta-analysis shown in tutorial_v2.qmd.
# The fitted object is intentionally local and Git-ignored. This script saves it
# only after the configured HMC diagnostics pass and writes all rendered inputs.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "--fit"
root <- normalizePath(".")
rds_path <- file.path(root, "Rdata", "tutorial_v2", "moura2021_BM_brms.rds")
out_dir <- file.path(root, "revision_checks", "moura_bm_brms_precompute_outputs")
candidate_path <- file.path(out_dir, "moura2021_BM_brms_candidate.rds")

assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

config <- list(
  chains = 4L,
  cores = 4L,
  # A 3,000-iteration candidate still left study-SD bulk ESS below 400.
  # Retain 4,000 draws per chain before deciding whether more sampling helps.
  iter = 6000L,
  warmup = 2000L,
  seed = 20260911L,
  # A first run at brms defaults had 6 divergences and 3,893/4,000 post-warmup
  # transitions at treedepth 10. At adapt_delta = 0.95, the longer candidate
  # had adequate mixing but one remaining divergence. These sampler-only
  # settings address that diagnosed geometry without changing the statistical
  # model or its priors.
  adapt_delta = 0.99,
  max_treedepth = 15L
)

key_parameters <- c(
  "b_Intercept",
  "sd_study.id__Intercept",
  "sd_effect.size.id__Intercept",
  "sd_species.id__Intercept",
  "sd_species.id.phy__Intercept"
)

parameter_labels <- c(
  b_Intercept = "pooled_mean",
  sd_study.id__Intercept = "study_sd",
  sd_effect.size.id__Intercept = "effect_size_sd",
  sd_species.id__Intercept = "species_nonphylogenetic_sd",
  sd_species.id.phy__Intercept = "species_phylogenetic_sd"
)

prepare_data <- function() {
  suppressPackageStartupMessages({
    library(ape)
    library(metadat)
    library(metafor)
  })
  dat <- dat.moura2021$dat
  dat$species.id.phy <- dat$species.id
  dat$effect.size.id <- factor(seq_len(nrow(dat)))
  dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  tree <- compute.brlen(dat.moura2021$tree)
  A <- vcv(tree, corr = TRUE)
  tip_order <- rownames(A)
  dat$species.id.phy <- factor(as.character(dat$species.id.phy), levels = tip_order)
  assert(!anyNA(dat$species.id.phy) && identical(levels(dat$species.id.phy), tip_order),
         "Species factor does not match the Brownian-motion correlation matrix.")
  list(dat = dat, A = A)
}

model_formula <- function() {
  brms::bf(
    yi | se(sqrt(vi), sigma = FALSE) ~ 1 +
      (1 | study.id) +
      (1 | effect.size.id) +
      (1 | species.id) +
      (1 | gr(species.id.phy, cov = A))
  )
}

model_priors <- function() {
  c(
    brms::set_prior("normal(0, 1)", class = "Intercept"),
    brms::set_prior("exponential(1)", class = "sd", group = "study.id"),
    brms::set_prior("exponential(1)", class = "sd", group = "effect.size.id"),
    brms::set_prior("exponential(1)", class = "sd", group = "species.id"),
    brms::set_prior("exponential(1)", class = "sd", group = "species.id.phy")
  )
}

summarise_vector <- function(x) {
  c(median = stats::median(x),
    q2.5 = unname(stats::quantile(x, 0.025)),
    q97.5 = unname(stats::quantile(x, 0.975)))
}

extract_diagnostics <- function(fit) {
  # brms 2.23.0 stores a stanfit object even with backend = "cmdstanr" here.
  # Extracting the unpermuted array retains the four chain identities needed
  # for rank-normalised R-hat and bulk/tail ESS.
  draws_all <- posterior::as_draws_array(
    rstan::extract(fit$fit, permuted = FALSE, inc_warmup = FALSE)
  )
  all_summary <- posterior::summarise_draws(
    draws_all, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  )
  names(all_summary)[names(all_summary) == "variable"] <- "parameter"
  all_summary <- as.data.frame(all_summary)
  names(all_summary)[names(all_summary) == "posterior::rhat"] <- "rhat"
  names(all_summary)[names(all_summary) == "posterior::ess_bulk"] <- "ess_bulk"
  names(all_summary)[names(all_summary) == "posterior::ess_tail"] <- "ess_tail"

  draws_key <- posterior::subset_draws(draws_all, variable = key_parameters)
  key_summary <- posterior::summarise_draws(
    draws_key, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  )
  names(key_summary)[names(key_summary) == "variable"] <- "parameter"
  key_summary <- as.data.frame(key_summary)
  names(key_summary)[names(key_summary) == "posterior::rhat"] <- "rhat"
  names(key_summary)[names(key_summary) == "posterior::ess_bulk"] <- "ess_bulk"
  names(key_summary)[names(key_summary) == "posterior::ess_tail"] <- "ess_tail"

  key_array <- posterior::as_draws_array(draws_key)
  for (chain in seq_len(dim(key_array)[2L])) {
    key_summary[[paste0("chain_", chain, "_mean")]] <- vapply(
      seq_along(key_parameters),
      function(i) mean(key_array[, chain, key_parameters[[i]]]), numeric(1)
    )
    key_summary[[paste0("chain_", chain, "_sd")]] <- vapply(
      seq_along(key_parameters),
      function(i) stats::sd(key_array[, chain, key_parameters[[i]]]), numeric(1)
    )
  }

  # nuts_params() supports both rstan and cmdstanr brmsfit backends.
  nuts <- as.data.frame(brms::nuts_params(fit))
  divergent <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  treedepth <- nuts$Value[nuts$Parameter == "treedepth__"]
  energy <- nuts[nuts$Parameter == "energy__", c("Chain", "Value")]
  energy_by_chain <- split(energy$Value, energy$Chain)
  bfmi <- vapply(energy_by_chain, function(e) {
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))

  hmc <- data.frame(
    max_rhat_all_parameters = max(all_summary$rhat[is.finite(all_summary$rhat)]),
    min_bulk_ess_all_parameters = min(all_summary$ess_bulk[is.finite(all_summary$ess_bulk)]),
    min_tail_ess_all_parameters = min(all_summary$ess_tail[is.finite(all_summary$ess_tail)]),
    divergent_transitions = divergent,
    max_observed_treedepth = max(treedepth),
    configured_max_treedepth = config$max_treedepth,
    min_bfmi = min(bfmi),
    stringsAsFactors = FALSE
  )
  bfmi_table <- data.frame(chain = as.integer(names(bfmi)), bfmi = unname(bfmi))
  list(all = all_summary, key = key_summary, hmc = hmc, bfmi = bfmi_table)
}

posterior_outputs <- function(fit, dat) {
  draws <- posterior::as_draws_df(posterior::as_draws_array(
    rstan::extract(fit$fit, pars = key_parameters,
                   permuted = FALSE, inc_warmup = FALSE)
  ))
  values <- list(
    pooled_mean = draws$b_Intercept,
    study_variance = draws$sd_study.id__Intercept^2,
    effect_size_variance = draws$sd_effect.size.id__Intercept^2,
    species_nonphylogenetic_variance = draws$sd_species.id__Intercept^2,
    species_phylogenetic_variance = draws$sd_species.id.phy__Intercept^2
  )
  estimates <- do.call(rbind, lapply(names(values), function(component) {
    stats <- summarise_vector(values[[component]])
    data.frame(component = component, estimate = stats[["median"]],
               lower_95 = stats[["q2.5"]], upper_95 = stats[["q97.5"]],
               interval = "95% credible interval", row.names = NULL)
  }))

  parameter_draws <- rbind(
    data.frame(value = values$pooled_mean, parameter = "Pooled mean (Fisher's Z)", panel = "Fixed effect"),
    data.frame(value = values$study_variance, parameter = "Study variance", panel = "Variance"),
    data.frame(value = values$effect_size_variance, parameter = "Effect-size variance", panel = "Variance"),
    data.frame(value = values$species_nonphylogenetic_variance,
               parameter = "Non-phylogenetic species variance", panel = "Variance"),
    data.frame(value = values$species_phylogenetic_variance,
               parameter = "Phylogenetic species variance", panel = "Variance")
  )
  parameter_draws$parameter <- factor(parameter_draws$parameter,
                                      levels = rev(unique(parameter_draws$parameter)))
  parameter_plot <- ggplot2::ggplot(parameter_draws, ggplot2::aes(x = value, y = parameter)) +
    tidybayes::stat_halfeye(.width = c(0.5, 0.95), point_interval = "median_qi",
                             fill = "#73A9AD", color = "#24535A") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
    ggplot2::facet_grid(panel ~ ., scales = "free_x", space = "free_y") +
    ggplot2::labs(x = "Posterior distribution (median, 50% and 95% credible intervals)", y = NULL) +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(file.path(out_dir, "moura_bm_brms_parameter_distributions.png"),
                  parameter_plot, width = 10, height = 7, dpi = 180)
  utils::write.csv(parameter_draws,
                   file.path(out_dir, "moura_bm_brms_parameter_draws.csv"), row.names = FALSE)

  set.seed(config$seed)
  yrep <- brms::posterior_predict(fit, ndraws = 500L)
  ppc <- rbind(
    data.frame(
      statistic = "mean", observed = mean(dat$yi),
      replicated_median = stats::median(rowMeans(yrep)),
      replicated_lower_95 = unname(stats::quantile(rowMeans(yrep), 0.025)),
      replicated_upper_95 = unname(stats::quantile(rowMeans(yrep), 0.975)),
      p_replicated_ge_observed = mean(rowMeans(yrep) >= mean(dat$yi))
    ),
    data.frame(
      statistic = "standard_deviation", observed = stats::sd(dat$yi),
      replicated_median = stats::median(apply(yrep, 1L, stats::sd)),
      replicated_lower_95 = unname(stats::quantile(apply(yrep, 1L, stats::sd), 0.025)),
      replicated_upper_95 = unname(stats::quantile(apply(yrep, 1L, stats::sd), 0.975)),
      p_replicated_ge_observed = mean(apply(yrep, 1L, stats::sd) >= stats::sd(dat$yi))
    )
  )

  ppc_path <- file.path(out_dir, "moura_bm_brms_ppc_ecdf_overlay.png")
  grDevices::png(ppc_path, width = 1800, height = 1200, res = 180)
  print(brms::pp_check(fit, type = "dens_overlay", ndraws = 100L))
  grDevices::dev.off()
  list(estimates = estimates, ppc = ppc, ppc_path = ppc_path)
}

frequentist_comparison <- function(dat, A, brms_estimates) {
  meta_fit <- metafor::rma.mv(
    yi, vi,
    random = list(~ 1 | study.id, ~ 1 | effect.size.id,
                  ~ 1 | species.id, ~ 1 | species.id.phy),
    R = list(species.id.phy = A), data = dat, sparse = TRUE,
    method = "REML", test = "t"
  )

  dat$g <- factor(1)
  VCV <- diag(dat$vi)
  rownames(VCV) <- colnames(VCV) <- dat$effect.size.id
  tmb_fit <- glmmTMB::glmmTMB(
    yi ~ 1 + equalto(0 + effect.size.id | g, VCV) +
      (1 | study.id) + (1 | species.id) +
      propto(0 + species.id.phy | g, A),
    data = dat, REML = TRUE
  )
  tmb_sd <- exp(tmb_fit$fit$par[names(tmb_fit$fit$par) == "theta"])
  assert(length(tmb_sd) == 3L,
         "Unexpected glmmTMB random-effect parameter count; inspect the current model mapping.")

  component <- c("pooled_mean", "study_variance", "effect_size_variance",
                 "species_nonphylogenetic_variance", "species_phylogenetic_variance")
  metafor_estimate <- c(meta_fit$b[[1L]], meta_fit$sigma2)
  glmmTMB_estimate <- c(
    stats::coef(summary(tmb_fit))$cond["(Intercept)", "Estimate"],
    tmb_sd[[1L]]^2, stats::sigma(tmb_fit)^2, tmb_sd[[2L]]^2, tmb_sd[[3L]]
  )
  brms_estimate <- brms_estimates$estimate[match(component, brms_estimates$component)]
  data.frame(
    component = component,
    metafor_REML_estimate = metafor_estimate,
    glmmTMB_REML_estimate = glmmTMB_estimate,
    brms_posterior_median = brms_estimate,
    brms_minus_metafor = brms_estimate - metafor_estimate,
    brms_minus_glmmTMB = brms_estimate - glmmTMB_estimate,
    stringsAsFactors = FALSE
  )
}

old_displayed_comparison <- function(brms_estimates) {
  # The old fit object is absent. These are the rounded means/SDs printed in
  # tutorial_v2.qmd before R9; variance values below are their squared SDs.
  old <- data.frame(
    component = c("pooled_mean", "study_variance", "effect_size_variance",
                  "species_nonphylogenetic_variance", "species_phylogenetic_variance"),
    old_displayed_value = c(0.37, 0.14^2, 0.12^2, 0.23^2, 0.27^2),
    stringsAsFactors = FALSE
  )
  new <- brms_estimates[, c("component", "estimate")]
  names(new)[2L] <- "new_posterior_median"
  out <- merge(old, new, by = "component", sort = FALSE)
  out$difference_new_minus_old_display <- out$new_posterior_median - out$old_displayed_value
  out$comparison_note <- "Old value is rounded display; variance is squared from displayed SD."
  out
}

write_outputs <- function(fit, dat, A) {
  diagnostics <- extract_diagnostics(fit)
  posterior <- posterior_outputs(fit, dat)
  comparison <- frequentist_comparison(dat, A, posterior$estimates)
  old_comparison <- old_displayed_comparison(posterior$estimates)

  write.csv(diagnostics$all, file.path(out_dir, "all_parameter_diagnostics.csv"), row.names = FALSE)
  write.csv(diagnostics$key, file.path(out_dir, "key_parameter_diagnostics.csv"), row.names = FALSE)
  write.csv(diagnostics$hmc, file.path(out_dir, "hmc_diagnostics.csv"), row.names = FALSE)
  write.csv(diagnostics$bfmi, file.path(out_dir, "bfmi_by_chain.csv"), row.names = FALSE)
  write.csv(posterior$estimates, file.path(out_dir, "posterior_estimates.csv"), row.names = FALSE)
  write.csv(posterior$ppc, file.path(out_dir, "posterior_predictive_check_summary.csv"), row.names = FALSE)
  write.csv(comparison, file.path(out_dir, "cross_package_comparison.csv"), row.names = FALSE)
  write.csv(old_comparison, file.path(out_dir, "old_displayed_brms_comparison.csv"), row.names = FALSE)

  list(diagnostics = diagnostics, posterior = posterior,
       comparison = comparison, old_comparison = old_comparison)
}

diagnostics_pass <- function(diagnostics) {
  key <- diagnostics$key
  hmc <- diagnostics$hmc
  all(key$rhat <= 1.01) && all(key$ess_bulk >= 400) && all(key$ess_tail >= 400) &&
    hmc$divergent_transitions == 0L &&
    hmc$max_observed_treedepth < hmc$configured_max_treedepth &&
    hmc$min_bfmi >= 0.3
}

write_manifest <- function() {
  writeLines(c(
    "Moura BM brms precomputation",
    "Formula: yi | se(sqrt(vi), sigma = FALSE) ~ 1 + (1 | study.id) + (1 | effect.size.id) + (1 | species.id) + (1 | gr(species.id.phy, cov = A))",
    "Known sampling variances: vi enters as fixed observation-level variance through se(sqrt(vi)); sigma is fixed to zero.",
    "Priors: Intercept normal(0, 1); each group-level SD exponential(1).",
    paste("Chains:", config$chains),
    paste("Cores:", config$cores),
    paste("Iterations per chain:", config$iter),
    paste("Warmup per chain:", config$warmup),
    paste("Seed:", config$seed),
    paste("adapt_delta:", config$adapt_delta),
    paste("max_treedepth:", config$max_treedepth),
    paste("R:", R.version.string),
    paste("brms:", as.character(packageVersion("brms"))),
    paste("cmdstanr:", as.character(packageVersion("cmdstanr"))),
    paste("CmdStan:", as.character(cmdstanr::cmdstan_version())),
    paste("posterior:", as.character(packageVersion("posterior"))),
    paste("metafor:", as.character(packageVersion("metafor"))),
    paste("glmmTMB:", as.character(packageVersion("glmmTMB")))
  ), file.path(out_dir, "manifest.txt"))
}

if (identical(mode, "--fit")) {
  assert(!file.exists(rds_path), paste("Refusing to overwrite existing RDS:", rds_path))
  assert(!file.exists(candidate_path),
         paste("A completed candidate fit is present; validate it before fitting again:", candidate_path))
  suppressPackageStartupMessages({
    library(brms)
    library(posterior)
    library(bayesplot)
    library(ggplot2)
    library(tidybayes)
    library(glmmTMB)
    library(cmdstanr)
    library(rstan)
  })
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prepared <- prepare_data()
  formula <- model_formula()
  priors <- model_priors()
  brms::validate_prior(priors, formula, data = prepared$dat,
                       data2 = list(A = prepared$A), family = gaussian())
  fit <- brms::brm(
    formula = formula, family = gaussian(), data = prepared$dat,
    data2 = list(A = prepared$A), prior = priors,
    backend = "cmdstanr", chains = config$chains, cores = config$cores,
    iter = config$iter, warmup = config$warmup, seed = config$seed,
    control = list(adapt_delta = config$adapt_delta,
                   max_treedepth = config$max_treedepth),
    save_pars = brms::save_pars(all = TRUE), refresh = 100L
  )
  # Preserve completed sampling before reporting steps, which may fail for API
  # reasons unrelated to the posterior simulation itself.
  saveRDS(fit, candidate_path)
  results <- write_outputs(fit, prepared$dat, prepared$A)
  write_manifest()
  assert(diagnostics_pass(results$diagnostics),
         "Sampler diagnostics did not meet validation criteria; diagnostics were saved but no RDS was written.")
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  assert(file.rename(candidate_path, rds_path), "Could not promote the validated candidate RDS.")
  cat("MOURA_BM_BRMS_PRECOMPUTATION_PASSED\n")
  quit(status = 0L)
}

if (identical(mode, "--candidate-check")) {
  assert(file.exists(candidate_path), paste("Missing candidate RDS:", candidate_path))
  assert(!file.exists(rds_path), paste("Final RDS already exists:", rds_path))
  suppressPackageStartupMessages({
    library(brms)
    library(posterior)
    library(bayesplot)
    library(ggplot2)
    library(tidybayes)
    library(glmmTMB)
    library(cmdstanr)
    library(rstan)
  })
  prepared <- prepare_data()
  fit <- readRDS(candidate_path)
  assert(inherits(fit, "brmsfit"), "Candidate object is not a brmsfit.")
  results <- write_outputs(fit, prepared$dat, prepared$A)
  write_manifest()
  assert(diagnostics_pass(results$diagnostics), "Candidate RDS fails diagnostic validation.")
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  assert(file.rename(candidate_path, rds_path), "Could not promote the validated candidate RDS.")
  cat("MOURA_BM_BRMS_CANDIDATE_VALIDATED\n")
  quit(status = 0L)
}

if (identical(mode, "--check")) {
  assert(file.exists(rds_path), paste("Missing RDS:", rds_path))
  suppressPackageStartupMessages({
    library(brms)
    library(posterior)
    library(bayesplot)
    library(ggplot2)
    library(tidybayes)
    library(glmmTMB)
    library(cmdstanr)
    library(rstan)
  })
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prepared <- prepare_data()
  fit <- readRDS(rds_path)
  assert(inherits(fit, "brmsfit"), "Saved object is not a brmsfit.")
  results <- write_outputs(fit, prepared$dat, prepared$A)
  write_manifest()
  assert(diagnostics_pass(results$diagnostics), "Saved RDS fails diagnostic validation.")
  cat("MOURA_BM_BRMS_EXISTING_RDS_VALIDATED\n")
  quit(status = 0L)
}

stop("Unknown mode: ", mode, call. = FALSE)
