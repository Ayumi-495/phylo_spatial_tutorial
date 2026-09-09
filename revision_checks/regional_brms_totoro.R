# Regional Spain brms audit for the cross-package comparison.
# Run on Totoro only. This script is separate from tutorial_v2.qmd.
# Model: known diagonal sampling variances + iid effect-size heterogeneity
# (brms sigma) + one exponential GP over recorded coordinate locations.

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript regional_brms_totoro.R <spain_prepared.rds> <output_dir>")
}
prepared_path <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260908L
chains <- 4L
cores <- 4L
threads_per_chain <- 10L
iter <- 3000L
warmup <- 1500L
adapt_delta <- 0.95
max_treedepth <- 12L

p <- readRDS(prepared_path)
dat <- p$dat
stopifnot(nrow(dat) == 186L, nlevels(dat$study_id) == 30L,
          nlevels(dat$site_id) == 32L,
          all(is.finite(dat$x_km)), all(is.finite(dat$y_km)),
          all(is.finite(dat$d_Hedges)), all(dat$var_Hedges > 0))

# Same observation-level model as the frequentist audit. With sigma = TRUE,
# brms uses sqrt(se_i^2 + sigma^2), so sigma is the iid effect-size SD and
# the known sampling SE remains a separate, fixed component.
form <- brms::bf(
  d_Hedges | se(sqrt(var_Hedges), sigma = TRUE) ~
    1 + gp(x_km, y_km, cov = "exponential", scale = FALSE)
)

# Explicit weakly informative priors on the common km scale. The exponential
# brms GP uses covariance sdgp^2 * exp(-distance / lscale), so lscale is the
# e-folding range in kilometres when scale = FALSE.
priors <- c(
  brms::prior(normal(0, 1), class = "Intercept"),
  brms::prior(student_t(3, 0, 1), class = "sigma"),
  brms::prior(student_t(3, 0, 1), class = "sdgp", coef = "gpx_kmy_km"),
  brms::prior(lognormal(log(50), 1), class = "lscale", coef = "gpx_kmy_km")
)

writeLines(c(
  paste("R:", R.version.string),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("cmdstanr:", as.character(packageVersion("cmdstanr"))),
  paste("CmdStan:", as.character(cmdstanr::cmdstan_version())),
  paste("seed:", seed), paste("chains:", chains), paste("cores:", cores),
  paste("threads_per_chain:", threads_per_chain), paste("iter:", iter),
  paste("warmup:", warmup), paste("adapt_delta:", adapt_delta),
  paste("max_treedepth:", max_treedepth),
  "data: Spain-labelled subset; 186 effect sizes, 30 studies, 32 recorded-coordinate locations.",
  "formula: d_Hedges | se(sqrt(var_Hedges), sigma=TRUE) ~ 1 + gp(x_km, y_km, cov=exponential, scale=FALSE).",
  "sampling model: Normal(mu, sqrt(se_i^2 + sigma^2)); sigma is iid effect-size SD.",
  "spatial model: exponential GP covariance sdgp^2 * exp(-distance/lscale) over unique recorded coordinate pairs.",
  "coordinates: WGS84 Lambert Conformal Conic projected x_km/y_km from the shared prepared regional audit file.",
  "No study-level intercept was included in this regional spatial-only comparison.",
  paste("priors:", paste(capture.output(print(priors)), collapse = " "))
), file.path(out_dir, "brms_metadata.txt"))

# Save generated Stan code before sampling so the exact parameterisation and
# observation-level variance decomposition remain inspectable if sampling fails.
stancode <- brms::make_stancode(form, data = dat, family = gaussian(), prior = priors,
                                backend = "cmdstanr")
writeLines(stancode, file.path(out_dir, "brms_stancode.stan"))

fit <- brms::brm(
  formula = form,
  data = dat,
  family = gaussian(),
  prior = priors,
  backend = "cmdstanr",
  chains = chains,
  cores = cores,
  threads = brms::threading(threads = threads_per_chain, static = TRUE),
  iter = iter,
  warmup = warmup,
  seed = seed,
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  refresh = 100
)

# Save the completed model immediately after sampling, before diagnostics or
# plotting, so a completed fit survives any downstream reporting failure.
saveRDS(fit, file.path(out_dir, "brms_spatial_only.rds"))

post <- posterior::as_draws_df(fit)
summ <- brms::posterior_summary(fit)
write.csv(as.data.frame(summ), file.path(out_dir, "brms_posterior_summary.csv"))

pick_summary <- function(pattern, label) {
  nm <- rownames(summ)
  hit <- grep(pattern, nm, value = TRUE)
  if (length(hit) != 1L) stop("Expected one summary row for ", label, "; got: ", paste(hit, collapse = ", "))
  x <- summ[hit, , drop = FALSE]
  ci_lo_col <- if ("l-95% CI" %in% colnames(summ)) "l-95% CI" else "Q2.5"
  ci_hi_col <- if ("u-95% CI" %in% colnames(summ)) "u-95% CI" else "Q97.5"
  rhat_col <- if ("Rhat" %in% colnames(summ)) "Rhat" else NA_character_
  bulk_col <- if ("Bulk_ESS" %in% colnames(summ)) "Bulk_ESS" else NA_character_
  tail_col <- if ("Tail_ESS" %in% colnames(summ)) "Tail_ESS" else NA_character_
  data.frame(parameter = label, estimate = x[1, "Estimate"],
             est_error = x[1, "Est.Error"],
             ci_lb = x[1, ci_lo_col], ci_ub = x[1, ci_hi_col],
             rhat = if (is.na(rhat_col)) NA_real_ else x[1, rhat_col],
             bulk_ess = if (is.na(bulk_col)) NA_real_ else x[1, bulk_col],
             tail_ess = if (is.na(tail_col)) NA_real_ else x[1, tail_col],
             stringsAsFactors = FALSE)
}

summary_rows <- rbind(
  pick_summary("^b_Intercept$", "pooled_mean"),
  pick_summary("^sigma$", "iid_effect_sd"),
  pick_summary("^sdgp_gpx_kmy_km$|^sdgp\\(gpx_kmy_km\\)$", "spatial_sd"),
  pick_summary("^lscale_gpx_kmy_km$|^lscale\\(gpx_kmy_km\\)$", "rho_km")
)
summary_rows$spatial_variance <- NA_real_
summary_rows$iid_effect_variance <- NA_real_
summary_rows$spatial_variance[summary_rows$parameter == "spatial_sd"] <-
  summary_rows$estimate[summary_rows$parameter == "spatial_sd"]^2
summary_rows$iid_effect_variance[summary_rows$parameter == "iid_effect_sd"] <-
  summary_rows$estimate[summary_rows$parameter == "iid_effect_sd"]^2
write.csv(summary_rows, file.path(out_dir, "brms_spatial_only_result.csv"), row.names = FALSE)

# HMC diagnostics: R-hat/ESS from posterior_summary, plus divergences and
# maximum treedepth from NUTS sampler parameters.
diag <- posterior::summarise_draws(post,
  posterior::rhat, posterior::ess_bulk, posterior::ess_tail)
diag_rhat <- diag[[grep("rhat$", names(diag), value = TRUE)[1L]]]
diag_bulk <- diag[[grep("ess_bulk$", names(diag), value = TRUE)[1L]]]
diag_tail <- diag[[grep("ess_tail$", names(diag), value = TRUE)[1L]]]
sampler <- posterior::as_draws_df(brms::nuts_params(fit))
divergences <- sum(sampler$Parameter == "divergent__" & sampler$Value > 0)
treedepth_values <- sampler$Value[sampler$Parameter == "treedepth__"]
max_treedepth_seen <- if (length(treedepth_values)) max(treedepth_values) else NA_real_
max_rhat <- max(diag_rhat, na.rm = TRUE)
min_bulk_ess <- min(diag_bulk, na.rm = TRUE)
min_tail_ess <- min(diag_tail, na.rm = TRUE)
diagnostics <- data.frame(
  chains = chains, cores = cores, threads_per_chain = threads_per_chain,
  max_rhat = max_rhat, min_bulk_ess = min_bulk_ess,
  min_tail_ess = min_tail_ess, divergences = divergences,
  max_treedepth_seen = max_treedepth_seen,
  max_treedepth = max_treedepth,
  rhat_ok = is.finite(max_rhat) && max_rhat < 1.01,
  ess_ok = is.finite(min_bulk_ess) && is.finite(min_tail_ess) &&
    min_bulk_ess >= 400 && min_tail_ess >= 400,
  divergences_ok = divergences == 0,
  treedepth_ok = is.na(max_treedepth_seen) || max_treedepth_seen < max_treedepth,
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(out_dir, "brms_diagnostics.csv"), row.names = FALSE)
write.csv(as.data.frame(diag), file.path(out_dir, "brms_draw_diagnostics.csv"), row.names = FALSE)

# Basic posterior predictive check saved as a PNG and as replicated draws.
pp <- brms::posterior_predict(fit, ndraws = min(200L, nrow(post)), seed = seed + 1L)
write.csv(pp, file.path(out_dir, "brms_posterior_predictive_draws.csv"), row.names = FALSE)
png(file.path(out_dir, "brms_pp_check_dens_overlay.png"), width = 1800, height = 1200, res = 180)
print(brms::pp_check(fit, type = "dens_overlay", ndraws = min(100L, nrow(post))))
dev.off()

writeLines(c(
  paste("max R-hat:", max_rhat), paste("minimum bulk ESS:", min_bulk_ess),
  paste("minimum tail ESS:", min_tail_ess), paste("divergences:", divergences),
  paste("maximum treedepth seen:", max_treedepth_seen),
  "The known sampling standard errors were combined with the additional iid heterogeneity as sqrt(se_i^2 + sigma^2), not multiplied by sigma.",
  "The exponential GP range is lscale in km because scale=FALSE and x_km/y_km were supplied."
), file.path(out_dir, "brms_diagnostics.txt"))

message("brms Spain regional fit and diagnostics saved.")
