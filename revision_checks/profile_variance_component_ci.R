#!/usr/bin/env Rscript
# Obtain one 95% profile-likelihood interval for a saved finalized rma.mv fit.
# Called independently for each variance component so completed intervals are
# saved immediately. This script never overwrites the saved model object.

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: profile_variance_component_ci.R <audit_dir> <model> <type> <index>")
}

audit_dir <- args[[1]]
model_name <- args[[2]]
type <- args[[3]]
index <- as.integer(args[[4]])
if (!type %in% c("sigma2", "tau2")) stop("type must be sigma2 or tau2")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
input <- file.path(root, "revision_checks", audit_dir, paste0(model_name, ".rds"))
out_dir <- file.path(root, "revision_checks", audit_dir, "variance_profile_ci")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(input)
started <- Sys.time()
ci <- do.call(confint, c(list(object = fit, level = 0.95, time = TRUE), setNames(list(index), type)))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
row <- as.data.frame(ci$random[1, , drop = FALSE])
names(row) <- c("estimate", "ci_lb", "ci_ub")
row$model <- model_name
row$parameter_type <- type
row$parameter_index <- index
row$elapsed_seconds <- elapsed
row$source_model <- input
row <- row[, c("model", "parameter_type", "parameter_index", "estimate", "ci_lb", "ci_ub", "elapsed_seconds", "source_model")]

stem <- sprintf("%s_%s_%d", model_name, type, index)
saveRDS(ci, file.path(out_dir, paste0(stem, ".rds")))
write.csv(row, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)
message("Saved ", stem)
