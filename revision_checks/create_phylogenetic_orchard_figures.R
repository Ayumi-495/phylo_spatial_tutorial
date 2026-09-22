#!/usr/bin/env Rscript

# Render the phylogenetic tutorial's frequentist orchard plots from saved fits.
# The saved models are produced by precompute_visualisation_intervals.R; this
# script does not fit or profile a model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(glmmTMB)
  library(metafor)
  library(orchaRd)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript revision_checks/create_phylogenetic_orchard_figures.R", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
inputs <- file.path(root, "revision_checks", "visualisation_interval_outputs")
fig_dir <- file.path(root, "figs", "tutorial")
stopifnot(dir.exists(fig_dir))

read_fit <- function(filename, class_name, n) {
  path <- file.path(inputs, filename)
  if (!file.exists(path)) stop("Missing precomputed fit: ", path, call. = FALSE)
  fit <- readRDS(path)
  if (!inherits(fit, class_name)) stop("Unexpected saved-model class: ", path, call. = FALSE)
  rows <- if (inherits(fit, "rma.mv")) fit$k else nrow(fit$frame)
  if (!identical(as.integer(rows), as.integer(n))) stop("Unexpected model row count: ", path, call. = FALSE)
  fit
}

check_interval <- function(observed, expected, tolerance, label) {
  if (!isTRUE(all.equal(as.numeric(observed), as.numeric(expected), tolerance = tolerance))) {
    stop("Orchard estimate or confidence interval disagrees with the saved interval artifact: ",
         label, call. = FALSE)
  }
}

saved_intervals <- function(filename) {
  read.csv(file.path(inputs, filename), check.names = FALSE, stringsAsFactors = FALSE)
}

save_orchard <- function(summary, filename, width, height, title, group = NULL) {
  plot <- orchaRd::orchard_plot(
    summary,
    mod = if (is.null(group)) "1" else group,
    group = if (is.null(group)) "study.id" else "id",
    xlab = "Effect size (Fisher's z)",
    k = TRUE,
    g = FALSE,
    k.pos = "right",
    legend.pos = "bottom.out",
    twig.size = 0.5,
    branch.size = 1.1,
    trunk.size = 0.5
  ) +
    labs(title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
  if (is.null(group)) plot <- plot + scale_x_discrete(labels = "Overall effect")
  ggsave(file.path(fig_dir, filename), plot, width = width, height = height,
         units = "in", dpi = 220, bg = "white")
}

moura_meta <- read_fit("moura_metafor_bm.rds", "rma.mv", 1828L)
moura_tmb <- read_fit("moura_glmmTMB_bm.rds", "glmmTMB", 1828L)
lim_meta <- read_fit("lim_metafor_bm_meta_regression.rds", "rma.mv", 170L)
moura_intervals <- saved_intervals("moura_frequentist_intervals.csv")
lim_intervals <- saved_intervals("lim_frequentist_intervals.csv")

# Saved rma.mv fits omit the formula attribute needed by orchaRd's plotting
# method. Annotate only in-memory copies; the saved fits are unchanged.
moura_meta$formula <- ~ 1
lim_meta$formula <- ~ environment

moura_meta_orchard <- orchaRd::mod_results(moura_meta, mod = "1", group = "study.id")
moura_meta_row <- moura_intervals[moura_intervals$package == "metafor" &
                                    moura_intervals$parameter == "intrcpt", ]
stopifnot(nrow(moura_meta_row) == 1L, nrow(moura_meta_orchard$mod_table) == 1L)
check_interval(moura_meta_orchard$mod_table[1L, c("estimate", "lowerCL", "upperCL")],
               moura_meta_row[1L, c("estimate", "ci_lb", "ci_ub")], 1e-5, "Moura metafor")
save_orchard(moura_meta_orchard, "moura_metafor_orchard.png", 9, 5.6,
             "Moura et al.: BM meta-analysis (metafor)")

# Use the same effect-size rows and sampling variances as the saved metafor fit.
# Conversion constructs an rma.mv-compatible plotting object, not a new fit.
for (field in c("yi", "study.id", "effect.size.id", "species.id", "species.id.phy")) {
  if (!identical(as.character(moura_tmb$frame[[field]]),
                 as.character(moura_meta$data[[field]]))) {
    stop("Moura glmmTMB and metafor rows disagree on ", field, call. = FALSE)
  }
}
moura_tmb_rma <- orchaRd::glmmTMB_to_rma(
  moura_tmb, yi = "yi", vi = "vi", data = moura_meta$data,
  measure = "GEN", test = "z"
)
moura_tmb_orchard <- orchaRd::mod_results(moura_tmb_rma, mod = "1", group = "study.id")
moura_tmb_row <- moura_intervals[moura_intervals$package == "glmmTMB" &
                                   moura_intervals$parameter == "(Intercept)", ]
stopifnot(nrow(moura_tmb_row) == 1L, nrow(moura_tmb_orchard$mod_table) == 1L)
check_interval(moura_tmb_orchard$mod_table[1L, c("estimate", "lowerCL", "upperCL")],
               moura_tmb_row[1L, c("estimate", "ci_lb", "ci_ub")], 1e-4, "Moura glmmTMB")
if (max(abs(as.numeric(moura_tmb_orchard$mod_table[1L, c("lowerPR", "upperPR")]) -
            as.numeric(moura_meta_orchard$mod_table[1L, c("lowerPR", "upperPR")]))) > 0.001) {
  stop("Converted Moura prediction interval disagrees with the matched metafor model.", call. = FALSE)
}
save_orchard(moura_tmb_orchard, "moura_glmmtmb_orchard.png", 9, 5.6,
             "Moura et al.: BM meta-analysis (glmmTMB)")

lim_meta_orchard <- orchaRd::mod_results(lim_meta, mod = "environment", group = "id")
lim_captive <- lim_meta_orchard$mod_table[lim_meta_orchard$mod_table$name == "Captive", ]
lim_intercept <- lim_intervals[lim_intervals$package == "metafor" &
                                lim_intervals$parameter == "intrcpt", ]
stopifnot(nrow(lim_captive) == 1L, nrow(lim_intercept) == 1L,
          nrow(lim_meta_orchard$mod_table) == 2L)
check_interval(lim_captive[1L, c("estimate", "lowerCL", "upperCL")],
               lim_intercept[1L, c("estimate", "ci_lb", "ci_ub")], 1e-5, "Lim metafor")
save_orchard(lim_meta_orchard, "lim_metafor_orchard.png", 7.5, 5,
             "Lim et al.: BM meta-regression (metafor)", group = "environment")

message("PHYLOGENETIC_ORCHARD_FIGURES_VALIDATED_AND_WRITTEN")
