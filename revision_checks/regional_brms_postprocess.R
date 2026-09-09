# Post-process a completed regional brms fit without refitting it.
suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript regional_brms_postprocess.R <brms_output_dir>")
}
out_dir <- args[[1L]]
fit_path <- file.path(out_dir, "brms_spatial_only.rds")
if (!file.exists(fit_path)) stop("Missing completed fit: ", fit_path)
fit <- readRDS(fit_path)
seed <- 20260908L

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

diag <- posterior::summarise_draws(post,
  posterior::rhat, posterior::ess_bulk, posterior::ess_tail)
diag_rhat <- diag[[grep("rhat$", names(diag), value = TRUE)[1L]]]
diag_bulk <- diag[[grep("ess_bulk$", names(diag), value = TRUE)[1L]]]
diag_tail <- diag[[grep("ess_tail$", names(diag), value = TRUE)[1L]]]
sampler <- as.data.frame(brms::nuts_params(fit))
divergences <- sum(sampler$Parameter == "divergent__" & sampler$Value > 0)
treedepth_values <- sampler$Value[sampler$Parameter == "treedepth__"]
max_treedepth_seen <- if (length(treedepth_values)) max(treedepth_values) else NA_real_
max_rhat <- max(diag_rhat, na.rm = TRUE)
min_bulk_ess <- min(diag_bulk, na.rm = TRUE)
min_tail_ess <- min(diag_tail, na.rm = TRUE)
diagnostics <- data.frame(
  max_rhat = max_rhat, min_bulk_ess = min_bulk_ess,
  min_tail_ess = min_tail_ess, divergences = divergences,
  max_treedepth_seen = max_treedepth_seen, max_treedepth = 12L,
  rhat_ok = is.finite(max_rhat) && max_rhat < 1.01,
  ess_ok = is.finite(min_bulk_ess) && is.finite(min_tail_ess) &&
    min_bulk_ess >= 400 && min_tail_ess >= 400,
  divergences_ok = divergences == 0,
  treedepth_ok = is.na(max_treedepth_seen) || max_treedepth_seen < 12L,
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(out_dir, "brms_diagnostics.csv"), row.names = FALSE)
write.csv(as.data.frame(diag), file.path(out_dir, "brms_draw_diagnostics.csv"), row.names = FALSE)

pp <- brms::posterior_predict(fit, ndraws = min(200L, nrow(post)), seed = seed + 1L)
write.csv(pp, file.path(out_dir, "brms_posterior_predictive_draws.csv"), row.names = FALSE)
png(file.path(out_dir, "brms_pp_check_dens_overlay.png"), width = 1800, height = 1200, res = 180)
print(brms::pp_check(fit, type = "dens_overlay", ndraws = min(100L, nrow(post))))
dev.off()

writeLines(c(
  paste("max R-hat:", max_rhat), paste("minimum bulk ESS:", min_bulk_ess),
  paste("minimum tail ESS:", min_tail_ess), paste("divergences:", divergences),
  paste("maximum treedepth seen:", max_treedepth_seen),
  "The known sampling standard errors were combined with the additional iid heterogeneity as sqrt(se_i^2 + sigma^2).",
  "The exponential GP range is lscale in km because scale=FALSE and x_km/y_km were supplied."
), file.path(out_dir, "brms_diagnostics.txt"))
message("brms post-processing saved.")
