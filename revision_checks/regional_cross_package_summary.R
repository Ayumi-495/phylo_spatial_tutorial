# Create a compact, inspectable regional cross-package comparison table.
suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript regional_cross_package_summary.R <audit_output_dir>")
out_dir <- args[[1L]]
mf <- read.csv(file.path(out_dir, "metafor_spatial_only_result.csv"), stringsAsFactors = FALSE)
tmb <- read.csv(file.path(out_dir, "glmmTMB_spatial_only_result.csv"), stringsAsFactors = FALSE)
fit <- readRDS(file.path(out_dir, "brms_output", "brms_spatial_only.rds"))
prep <- readRDS(file.path(out_dir, "spain_prepared.rds"))
draws <- as.data.frame(posterior::as_draws_df(fit))

q <- function(x) as.numeric(quantile(x, c(0.025, 0.5, 0.975), names = FALSE))
b_mean <- q(draws$b_Intercept)
b_iid <- q(draws$sigma^2)
b_sp <- q(draws$sdgp_gpx_kmy_km^2)
b_rho <- q(draws$lscale_gpx_kmy_km)

comparison <- data.frame(
  package = c("metafor", "glmmTMB", "brms"),
  n = c(mf$n, tmb$n, nrow(fit$data)),
  studies = c(mf$studies, tmb$studies, nlevels(prep$dat$study_id)),
  sites = c(mf$sites, tmb$sites, nrow(unique(fit$data[c("x_km", "y_km")]))),
  pooled_mean = c(mf$mean, tmb$mean, b_mean[2]),
  mean_lb = c(mf$ci_lb, tmb$ci_lb, b_mean[1]),
  mean_ub = c(mf$ci_ub, tmb$ci_ub, b_mean[3]),
  iid_effect_variance = c(mf$iid_effect_variance, tmb$iid_effect_variance, b_iid[2]),
  iid_effect_variance_lb = c(NA, NA, b_iid[1]),
  iid_effect_variance_ub = c(NA, NA, b_iid[3]),
  spatial_variance = c(mf$spatial_variance, tmb$spatial_variance, b_sp[2]),
  spatial_variance_lb = c(NA, NA, b_sp[1]),
  spatial_variance_ub = c(NA, NA, b_sp[3]),
  spatial_sd = c(tmb$spatial_sd, tmb$spatial_sd, b_sp[2]^0.5),
  rho_km = c(mf$rho_km, tmb$rho_km, b_rho[2]),
  rho_lb_km = c(NA, NA, b_rho[1]),
  rho_ub_km = c(NA, NA, b_rho[3]),
  logLik_REML = c(mf$logLik_REML, tmb$logLik_REML, NA),
  AIC_REML = c(mf$AIC_REML, tmb$AIC_REML, NA),
  diagnostic_status = c(
    "completed; Matrix class-coercion warning",
    "optimizer code 0; pdHess TRUE; diagnose() Wald caution",
    "max R-hat 1.0023; min bulk ESS 1121; min tail ESS 1708; divergences 0; max treedepth 7/12; PPC saved"
  ),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(out_dir, "regional_cross_package_comparison.csv"), row.names = FALSE)

writeLines(c(
  "The metafor and glmmTMB rows use REML and the same known diagonal sampling VCV; metafor used REMLf=FALSE to match the equalto likelihood convention.",
  "The brms row is posterior inference, not a REML/AIC fit; its intervals are 95% credible intervals.",
  "brms sigma^2 and sdgp^2/rho summaries are posterior quantiles transformed draw-by-draw, not squares of interval endpoints.",
  "The brms GP has 32 unique coordinate inputs because gp() groups repeated x_km/y_km values; no study-level intercept was included."
), file.path(out_dir, "regional_cross_package_comparison_notes.txt"))
