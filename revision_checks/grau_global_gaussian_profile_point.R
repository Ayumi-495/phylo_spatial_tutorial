# Run and save one fixed-parameter profile-likelihood point for a saved rma.mv model.
# Independent points can be distributed across a bounded number of Totoro workers.

suppressPackageStartupMessages({
  library(metafor)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript totoro_profile_point.R <model_rds> <data_csv> <output_dir> <model_name> <tau2|rho> <value>")
}
model_rds <- args[[1]]
data_csv <- args[[2]]
out_dir <- args[[3]]
model_name <- args[[4]]
component <- args[[5]]
value <- as.numeric(args[[6]])
stopifnot(component %in% c("tau2", "rho"), is.finite(value), value >= 0)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The saved rma.mv call refers to dat and distance_km. Reconstruct these objects
# identically so each profile point is self-contained and safely distributable.
dat <- read.csv(data_csv, stringsAsFactors = FALSE)
dat <- dat[complete.cases(dat$latitude, dat$longitude), ]
stopifnot(nrow(dat) == 2361L, all(dat$var_Hedges > 0))
dat$effect_id <- factor(seq_len(nrow(dat)))
dat$study_id <- factor(dat$study_id)
dat$site_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
site_levels <- sort(unique(dat$site_key))
dat$site_id <- factor(dat$site_key, levels = site_levels)
dat$const <- factor("all_sites")
stopifnot(nlevels(dat$study_id) == 393L, nlevels(dat$site_id) == 383L)

site_lookup <- unique(dat[c("site_id", "site_key", "latitude", "longitude")])
site_lookup <- site_lookup[match(site_levels, site_lookup$site_key), ]
coordinates_lonlat <- as.matrix(site_lookup[c("longitude", "latitude")])
distance_km <- geosphere::distm(coordinates_lonlat, fun = geosphere::distGeo) / 1000
rownames(distance_km) <- colnames(distance_km) <- site_levels
stopifnot(identical(rownames(distance_km), levels(dat$site_id)))
stopifnot(identical(colnames(distance_km), levels(dat$site_id)))

# The saved audited fit was created from a preparation object named `p`.
# Recreate its distance member so update() can evaluate the original call.
p <- list(dat = dat, distance_km = distance_km)

fit <- readRDS(model_rds)
started <- proc.time()[["elapsed"]]
warning_messages <- character()
result <- withCallingHandlers(
  tryCatch({
    if (component == "tau2") {
      update(fit, tau2 = value)
    } else {
      update(fit, rho = value)
    }
  }, error = function(e) e),
  warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
elapsed <- proc.time()[["elapsed"]] - started

safe_value <- gsub("[^0-9A-Za-z_.-]", "_", format(value, scientific = FALSE, trim = TRUE))
path <- file.path(out_dir, paste0(model_name, "_", component, "_", safe_value, ".csv"))
if (inherits(result, "error")) {
  row <- data.frame(model = model_name, component = component, value = value,
                    logLik_REML = NA_real_, AIC_REML = NA_real_,
                    status = "error",
                    elapsed_seconds = elapsed,
                    warnings = paste(c(conditionMessage(result), unique(warning_messages)), collapse = " | "))
} else {
  row <- data.frame(model = model_name, component = component, value = value,
                    logLik_REML = as.numeric(result$fit.stats["ll", "REML"]),
                    AIC_REML = as.numeric(result$fit.stats["AIC", "REML"]),
                    status = "completed_no_explicit_optimizer_status",
                    elapsed_seconds = elapsed,
                    warnings = paste(unique(warning_messages), collapse = " | "))
}
write.csv(row, path, row.names = FALSE)
message(path)
