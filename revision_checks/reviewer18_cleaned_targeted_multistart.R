# Limited fixed-tau2 multi-start diagnostic for suspicious spatial-only profile
# points. This is not a broad free-model starting-value search.

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript reviewer18_cleaned_targeted_multistart.R <prepared_rds> <primary_rds> <output_dir> <fixed_tau2> <rho_start>")
}
p <- readRDS(args[[1L]])
primary <- readRDS(args[[2L]])
out_dir <- args[[3L]]
fixed_tau2 <- as.numeric(args[[4L]])
rho_start <- as.numeric(args[[5L]])
stopifnot(inherits(primary, "rma.mv"), primary$k == 2355L,
          fixed_tau2 > 0, rho_start > 0,
          identical(rownames(p$distance_km), levels(p$dat$site_id)),
          identical(colnames(p$distance_km), levels(p$dat$site_id)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

warnings <- character()
started <- proc.time()[["elapsed"]]
fit <- withCallingHandlers(
  tryCatch(rma.mv(
    yi = d_Hedges, V = var_Hedges,
    random = list(~1 | effect_id, ~site_id | const),
    struct = "SPEXP", dist = list(site_id = p$distance_km),
    tau2 = fixed_tau2,
    data = p$dat, method = "REML", test = "t", sparse = TRUE,
    control = list(sigma2.init = primary$sigma2, rho.init = rho_start)
  ), error = function(e) e),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
elapsed <- proc.time()[["elapsed"]] - started

if (inherits(fit, "error")) {
  row <- data.frame(fixed_tau2 = fixed_tau2, rho_start = rho_start,
                    final_rho_km = NA_real_, logLik_REML = NA_real_,
                    status = "error", elapsed_seconds = elapsed,
                    warnings = paste(c(conditionMessage(fit), unique(warnings)), collapse = " | "))
} else {
  row <- data.frame(fixed_tau2 = fixed_tau2, rho_start = rho_start,
                    final_rho_km = fit$rho[1], logLik_REML = fit$fit.stats["ll", "REML"],
                    status = "completed_no_explicit_optimizer_status",
                    elapsed_seconds = elapsed,
                    warnings = paste(unique(warnings), collapse = " | "))
  saveRDS(fit, file.path(out_dir, sprintf("spatial_only_tau2_%g_rhostart_%g.rds",
                                         fixed_tau2, rho_start)))
}
write.csv(row, file.path(out_dir, sprintf("spatial_only_tau2_%g_rhostart_%g.csv",
                                          fixed_tau2, rho_start)), row.names = FALSE)
print(row)
