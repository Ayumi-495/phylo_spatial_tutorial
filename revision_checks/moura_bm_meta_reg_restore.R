#!/usr/bin/env Rscript

# Precompute and validate the exact metafor BM meta-regression displayed in tutorial_v2.qmd.
# This script was first used for recovery, and is retained as the permanent
# precomputation recipe for moura2021_BM_meta_reg.rds.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "--fit"
root <- normalizePath(".")
target <- file.path(root, "Rdata", "tutorial_v2", "moura2021_BM_meta_reg.rds")
out_dir <- file.path(root, "revision_checks", "moura_bm_meta_reg_restore_outputs")
validation_path <- file.path(out_dir, "moura2021_BM_meta_reg_validation.csv")
manifest_path <- file.path(out_dir, "moura2021_BM_meta_reg_manifest.txt")
inventory_path <- file.path(out_dir, "tutorial_readRDS_inventory.csv")
render_log_path <- file.path(out_dir, "executing_render_after_moura_bm_restore.log")
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)

inventory <- function() {
  qmd <- paste(readLines(file.path(root, "tutorial_v2.qmd"), warn = FALSE), collapse = "\n")
  read_pattern <- 'readRDS\\(here\\("Rdata", "tutorial_v2", "([^"]+\\.rds)"\\)\\)'
  path_pattern <- 'here\\("Rdata", "tutorial_v2", "([^"]+\\.rds)"\\)'
  read_paths <- regmatches(qmd, gregexpr(read_pattern, qmd, perl = TRUE))[[1L]]
  all_paths <- regmatches(qmd, gregexpr(path_pattern, qmd, perl = TRUE))[[1L]]
  read_files <- sub(read_pattern, "\\1", read_paths, perl = TRUE)
  files <- unique(sub(path_pattern, "\\1", all_paths, perl = TRUE))
  files <- files[nzchar(files)]
  data.frame(file = files,
    path = file.path("Rdata", "tutorial_v2", files),
    reference_type = ifelse(files %in% read_files, "readRDS", "path_only"),
    exists = file.exists(file.path(root, "Rdata", "tutorial_v2", files)),
    stringsAsFactors = FALSE)
}

if (identical(mode, "--inventory-check")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table <- inventory()
  write.csv(table, inventory_path, row.names = FALSE)
  cat("MOURA_BM_META_REG_RDS_INVENTORY_COMPLETED\n")
  quit(status = 0L)
}

if (identical(mode, "--main-copy-check")) {
  main_root <- sub("_revision$", "", root)
  main_target <- file.path(main_root, "Rdata", "tutorial_v2", "moura2021_BM_meta_reg.rds")
  assert(file.exists(target) && file.exists(main_target), "Revision or original-main RDS copy is missing.")
  assert(identical(unname(tools::md5sum(target)), unname(tools::md5sum(main_target))),
         "Revision and original-main RDS copies differ.")
  cat("MOURA_BM_META_REG_MAIN_COPY_MATCHES\n")
  quit(status = 0L)
}

if (identical(mode, "--render-check")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result <- system2("quarto", c("render", "tutorial_v2.qmd", "--to", "html"), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  writeLines(c(paste("exit_status:", status), result), render_log_path)
  if (identical(status, 0L)) {
    cat("MOURA_BM_META_REG_RENDER_CHECK_COMPLETED: FULL_RENDER_SUCCEEDED\n")
  } else {
    rds_lines <- result[grepl("\\.rds", result)]
    assert(length(rds_lines) > 0L, "Executing render failed without identifying an RDS dependency.")
    cat("MOURA_BM_META_REG_RENDER_CHECK_COMPLETED: NEXT_RDS_BOUNDARY\n")
    cat(paste(rds_lines, collapse = "\n"), "\n")
  }
  quit(status = 0L)
}

expected <- c(
  intercept = 0.3562,
  temporally_pooled_yes = 0.0395,
  study_variance = 0.0194,
  effect_size_variance = 0.0145,
  species_nonphylogenetic_variance = 0.0540,
  species_phylogenetic_variance = 0.0520,
  REML_logLik = -165.9738,
  AIC = 343.9476
)

fit_model <- function() {
  suppressPackageStartupMessages({ library(ape); library(metadat); library(metafor) })
  dat <- dat.moura2021$dat
  dat$species.id.phy <- dat$species.id
  dat$effect.size.id <- factor(seq_len(nrow(dat)))
  dat <- escalc(measure = "ZCOR", ri = ri, ni = ni, data = dat)
  tree <- compute.brlen(dat.moura2021$tree)
  A <- vcv(tree, corr = TRUE)
  tip_order <- rownames(A)
  dat$species.id.phy <- factor(as.character(dat$species.id.phy), levels = tip_order)
  assert(!anyNA(dat$species.id.phy) && identical(levels(dat$species.id.phy), tip_order),
         "Species factor does not match BM correlation-matrix order.")
  rma.mv(yi, vi,
    mods = ~ temporally.pooled,
    random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
    R = list(species.id.phy = A),
    data = dat,
    sparse = TRUE,
    method = "REML")
}

metrics <- function(fit) {
  data.frame(
    metric = names(expected),
    expected_tutorial_display = unname(expected),
    observed = c(fit$beta[[1L]], fit$beta[[2L]], fit$sigma2[[1L]], fit$sigma2[[2L]],
                 fit$sigma2[[3L]], fit$sigma2[[4L]], fit$fit.stats["ll", "REML"],
                 fit$fit.stats["AIC", "REML"]),
    stringsAsFactors = FALSE)
}

if (identical(mode, "--fit")) {
  assert(!file.exists(target), paste("Refusing to overwrite existing target:", target))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- fit_model()
  table <- metrics(fit)
  table$absolute_difference <- abs(table$observed - table$expected_tutorial_display)
  table$within_display_rounding <- table$absolute_difference < 5e-4
  assert(all(table$within_display_rounding), "Refit does not reproduce the tutorial values to displayed precision.")
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, target)
  write.csv(table, validation_path, row.names = FALSE)
  writeLines(c(
    "Saved object: moura2021_BM_meta_reg",
    "Object class: rma.mv",
    "Source fit name in tutorial code: phylo_eg1.1_meta_BM",
    "Relationship: the saved object is that same fitted model, saved under the loader name moura2021_BM_meta_reg.",
    "Specification: rma.mv(yi, vi, mods = ~ temporally.pooled; random study.id, effect.size.id, species.id, species.id.phy; R species.id.phy = vcv(compute.brlen(dat.moura2021$tree), corr = TRUE); sparse = TRUE; method = REML).",
    paste("R version:", R.version.string),
    paste("metafor version:", as.character(packageVersion("metafor"))),
    paste("ape version:", as.character(packageVersion("ape"))),
    paste("metadat version:", as.character(packageVersion("metadat")))
  ), manifest_path)
  cat("MOURA_BM_META_REG_RDS_FIT_COMPLETED\n")
  quit(status = 0L)
}

if (identical(mode, "--validate-existing")) {
  assert(file.exists(target), "Missing restored RDS.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- readRDS(target)
  assert(inherits(fit, "rma.mv"), "Restored object is not an rma.mv fit.")
  table <- metrics(fit)
  table$absolute_difference <- abs(table$observed - table$expected_tutorial_display)
  table$within_display_rounding <- table$absolute_difference < 5e-4
  assert(all(table$within_display_rounding), "Saved fit does not reproduce tutorial values to displayed precision.")
  write.csv(table, validation_path, row.names = FALSE)
  cat("MOURA_BM_META_REG_EXISTING_RDS_VALIDATED\n")
  quit(status = 0L)
}

if (identical(mode, "--check")) {
  assert(file.exists(target), "Missing restored RDS.")
  assert(file.exists(validation_path) && file.exists(manifest_path), "Missing restoration validation artifacts.")
  fit <- readRDS(target)
  assert(inherits(fit, "rma.mv"), "Restored object is not an rma.mv fit.")
  table <- metrics(fit)
  assert(all(abs(table$observed - table$expected_tutorial_display) < 5e-4),
         "Saved fit does not reproduce tutorial values to displayed precision.")
  cat("MOURA_BM_META_REG_RDS_RESTORE_CHECKS_PASSED\n")
  quit(status = 0L)
}

stop("Unknown mode: ", mode, call. = FALSE)
