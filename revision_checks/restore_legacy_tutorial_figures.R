#!/usr/bin/env Rscript

# Historical compatibility entry point. The former implementation recreated
# component panels from point estimates alone, so re-running it could silently
# replace interval-bearing public figures with incomplete ones. Public tutorial
# result figures are now generated only from the saved interval artifacts.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
source(file.path(root, "revision_checks", "create_visualisation_interval_figures.R"))
