# One fixed-parameter profile-likelihood point for a saved published_cleaned
# SPEXP rma.mv model. Independent invocations are safe to parallelize.

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript reviewer18_cleaned_profile_point.R <prepared_rds> <model_rds> <output_dir> <model> <tau2|rho> <value>")
}
prepared_rds <- args[[1L]]
model_rds <- args[[2L]]
out_dir <- args[[3L]]
model_name <- args[[4L]]
component <- args[[5L]]
value <- as.numeric(args[[6L]])
stopifnot(model_name %in% c("spatial_only", "combined"),
          component %in% c("tau2", "rho"), is.finite(value), value >= 0,
          file.exists(prepared_rds), file.exists(model_rds))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The saved model call refers to p$dat and p$distance_km. Recreate p exactly
# from the saved preparation object before update() evaluates that call.
p <- readRDS(prepared_rds)
stopifnot(nrow(p$dat) == 2355L, nlevels(p$dat$study_id) == 390L,
          nlevels(p$dat$site_id) == 380L,
          identical(rownames(p$distance_km), levels(p$dat$site_id)),
          identical(colnames(p$distance_km), levels(p$dat$site_id)))
fit <- readRDS(model_rds)
stopifnot(inherits(fit, "rma.mv"), fit$k == 2355L,
          isTRUE(all.equal(as.numeric(fit$vi), p$dat$var_Hedges)),
          identical(dim(fit$X), c(2355L, 1L)),
          all(as.numeric(fit$X) == 1))

warnings <- character()
started <- proc.time()[["elapsed"]]
result <- withCallingHandlers(
  tryCatch({
    if (component == "tau2") update(fit, tau2 = value) else update(fit, rho = value)
  }, error = function(e) e),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
elapsed <- proc.time()[["elapsed"]] - started

if (inherits(result, "error")) {
  row <- data.frame(
    model = model_name, component = component, value = value,
    mean = NA_real_, effect_variance = NA_real_, study_variance = NA_real_,
    spatial_variance = NA_real_, rho_km = NA_real_,
    logLik_REML = NA_real_, AIC_REML = NA_real_, status = "error",
    elapsed_seconds = elapsed,
    warnings = paste(c(conditionMessage(result), unique(warnings)), collapse = " | "),
    stringsAsFactors = FALSE
  )
} else {
  row <- data.frame(
    model = model_name, component = component, value = value,
    mean = as.numeric(result$b[1]), effect_variance = result$sigma2[1],
    study_variance = if (model_name == "combined") result$sigma2[2] else NA_real_,
    spatial_variance = result$tau2[1], rho_km = result$rho[1],
    logLik_REML = as.numeric(result$fit.stats["ll", "REML"]),
    AIC_REML = as.numeric(result$fit.stats["AIC", "REML"]),
    status = "completed_no_explicit_optimizer_status",
    elapsed_seconds = elapsed,
    warnings = paste(unique(warnings), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

safe_value <- gsub("[^0-9A-Za-z_.-]", "_", format(value, scientific = FALSE,
                                                       trim = TRUE, digits = 15))
path <- file.path(out_dir, paste0(model_name, "_", component, "_", safe_value, ".csv"))
write.csv(row, path, row.names = FALSE)
cat(path, "\n")
