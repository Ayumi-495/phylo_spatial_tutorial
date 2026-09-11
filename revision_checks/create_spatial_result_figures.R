#!/usr/bin/env Rscript
# Recreate tutorial result figures from saved, finalized spatial-audit objects.
# This script deliberately does not fit or profile any model.

suppressPackageStartupMessages({
  library(metafor)
  library(orchaRd)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "figs", "tutorial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_model <- function(...) readRDS(file.path(root, ...))

# The audit objects were saved without the formula slot used by orchaRd.
# Restoring the intercept-only formula in memory enables plotting only; it does
# not modify the fitted object on disk or trigger a refit.
orchard_panel <- function(model, x_limits, title) {
  model$formula <- ~ 1
  # Some historical Scholer study labels are not valid UTF-8. The plot only
  # uses this variable to group effect sizes, so replace display labels with a
  # one-to-one ASCII factor in memory while preserving group membership.
  study_labels <- as.character(model$data$study_id)
  if (anyNA(iconv(study_labels, from = "", to = "UTF-8", sub = NA))) {
    model$data$study_id <- factor(sprintf("study_%03d", as.integer(factor(study_labels))))
  }
  orchaRd::orchard_plot(model, group = "study_id", xlab = "Effect size") +
    scale_color_manual(values = "#CDAD00") +
    scale_fill_manual(values = "#FFD700") +
    scale_y_continuous(limits = x_limits) +
    coord_flip() +
    labs(title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

profile_component <- function(audit_dir, model, type, index, label) {
  path <- file.path(root, "revision_checks", audit_dir, "variance_profile_ci",
                    sprintf("%s_%s_%d.csv", model, type, index))
  if (!file.exists(path)) stop("Missing saved profile interval: ", path)
  ci_object <- readRDS(sub("\\.csv$", ".rds", path))
  read_csv(path, show_col_types = FALSE) |>
    transmute(component = label, variance = estimate, ci_lb, ci_ub,
              ci_lb_sign = ci_object$lb.sign, ci_ub_sign = ci_object$ub.sign)
}

variance_panel <- function(dat, levels, y_limit) {
  dat <- dat |>
    mutate(component = factor(component, levels = levels)) |>
    filter(!is.na(variance)) |>
    mutate(ci_ub_plot = pmin(ci_ub, y_limit),
           upper_truncated = ci_ub_sign == ">" | ci_ub > y_limit)
  ggplot(dat, aes(x = component, y = variance)) +
    geom_point(size = 3, colour = "#528B8B") +
    geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub_plot), width = 0.22, colour = "#528B8B") +
    geom_segment(data = filter(dat, upper_truncated),
                 aes(x = component, xend = component, y = y_limit * 0.90, yend = y_limit),
                 inherit.aes = FALSE, colour = "#528B8B",
                 arrow = grid::arrow(length = grid::unit(0.16, "cm"))) +
    geom_text(aes(label = sprintf("%.3f", variance)), hjust = -0.15, size = 3.4) +
    coord_flip(clip = "off") +
    scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = NULL, y = "Estimated variance") +
    theme_classic(base_size = 12) +
    theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
}

save_two_panel <- function(model, components, levels, x_limits, y_limit, title, filename) {
  p1 <- orchard_panel(model, x_limits, title)
  p2 <- variance_panel(components, levels, y_limit)
  ggsave(file.path(out_dir, filename), p1 / p2,
         width = 8.5, height = 7.0, dpi = 220, bg = "white")
}

model_names <- c(
  unstructured_only = "Unstructured-only model",
  spatial_only = "Spatial-only model",
  combined = "Combined model"
)

# Grau-Andres: published_cleaned primary analysis.
grau_dir <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")
grau <- lapply(names(model_names), function(nm) readRDS(file.path(grau_dir, paste0(nm, ".rds"))))
names(grau) <- names(model_names)
grau_components <- list(
  unstructured_only = bind_rows(
    profile_component("reviewer18_influential_effects_outputs", "unstructured_only", "sigma2", 1, "Effect-size"),
    profile_component("reviewer18_influential_effects_outputs", "unstructured_only", "sigma2", 2, "Study")),
  spatial_only = bind_rows(
    profile_component("reviewer18_influential_effects_outputs", "spatial_only", "sigma2", 1, "Effect-size"),
    profile_component("reviewer18_influential_effects_outputs", "spatial_only", "tau2", 1, "Spatial")),
  combined = bind_rows(
    profile_component("reviewer18_influential_effects_outputs", "combined", "sigma2", 1, "Effect-size"),
    profile_component("reviewer18_influential_effects_outputs", "combined", "sigma2", 2, "Study"),
    profile_component("reviewer18_influential_effects_outputs", "combined", "tau2", 1, "Spatial"))
)
for (nm in names(model_names)) {
  save_two_panel(grau[[nm]], grau_components[[nm]], c("Effect-size", "Study", "Spatial"),
                 c(-20, 20), 1.50, paste("Grau-Andrés:", model_names[[nm]]),
                 paste0("grau_andres_", nm, "_results.png"))
}

# Scholer: the same component order and axes are retained across all three
# figures so that the three specifications can be compared visually.
scholer_dir <- file.path(root, "revision_checks", "scholer_spatial_audit_outputs")
scholer <- lapply(names(model_names), function(nm) readRDS(file.path(scholer_dir, paste0(nm, ".rds"))))
names(scholer) <- names(model_names)
scholer_components <- list(
  unstructured_only = bind_rows(
    profile_component("scholer_spatial_audit_outputs", "unstructured_only", "sigma2", 1, "Effect-size"),
    profile_component("scholer_spatial_audit_outputs", "unstructured_only", "sigma2", 2, "Study")),
  spatial_only = bind_rows(
    profile_component("scholer_spatial_audit_outputs", "spatial_only", "sigma2", 1, "Effect-size"),
    profile_component("scholer_spatial_audit_outputs", "spatial_only", "tau2", 1, "Spatial")),
  combined = bind_rows(
    profile_component("scholer_spatial_audit_outputs", "combined", "sigma2", 1, "Effect-size"),
    profile_component("scholer_spatial_audit_outputs", "combined", "sigma2", 2, "Study"),
    profile_component("scholer_spatial_audit_outputs", "combined", "tau2", 1, "Spatial"))
)
for (nm in names(model_names)) {
  save_two_panel(scholer[[nm]], scholer_components[[nm]], c("Effect-size", "Study", "Spatial"),
                 c(-2, 4), 0.52, paste("Scholer et al.:", model_names[[nm]]),
                 paste0("scholer_", nm, "_results.png"))
}

# Spain: a direct, matched comparison uses the same two-panel grammar. The
# upper panel distinguishes frequentist 95% CIs from the brms 95% CrI.
spain_dir <- file.path(root, "revision_checks", "regional_cross_package_audit_outputs")
mf <- read_csv(file.path(spain_dir, "metafor_spatial_only_result.csv"), show_col_types = FALSE)
gt <- read_csv(file.path(spain_dir, "glmmTMB_spatial_only_result.csv"), show_col_types = FALSE)
br <- read_csv(file.path(spain_dir, "brms_output", "brms_spatial_only_result.csv"), show_col_types = FALSE)
br_mean <- filter(br, parameter == "pooled_mean")
br_effect <- filter(br, parameter == "iid_effect_sd")
br_spatial <- filter(br, parameter == "spatial_sd")

spain_mean <- tibble(
  package = factor(c("metafor", "glmmTMB", "brms"), levels = c("metafor", "glmmTMB", "brms")),
  estimate = c(mf$mean, gt$mean, br_mean$estimate),
  lower = c(mf$ci_lb, gt$ci_lb, br_mean$ci_lb),
  upper = c(mf$ci_ub, gt$ci_ub, br_mean$ci_ub),
  interval = c("95% CI", "95% CI", "95% CrI")
)
p1_spain <- ggplot(spain_mean, aes(y = package, x = estimate, colour = interval)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.75) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("95% CI" = "#528B8B", "95% CrI" = "#CDAD00")) +
  scale_x_continuous(limits = c(-0.65, 0.55)) +
  labs(title = "Spain regional subset: spatial-only common target", x = "Pooled effect size", y = NULL, colour = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")

spain_var <- tibble(
  package = factor(rep(c("metafor", "glmmTMB", "brms"), each = 2), levels = c("metafor", "glmmTMB", "brms")),
  component = factor(rep(c("Effect-size", "Spatial"), 3), levels = c("Effect-size", "Spatial")),
  estimate = c(mf$iid_effect_variance, mf$spatial_variance,
               gt$iid_effect_variance, gt$spatial_variance,
               br_effect$iid_effect_variance, br_spatial$spatial_variance),
  lower = c(NA, NA, NA, NA, br_effect$ci_lb^2, br_spatial$ci_lb^2),
  upper = c(NA, NA, NA, NA, br_effect$ci_ub^2, br_spatial$ci_ub^2)
)
p2_spain <- ggplot(spain_var, aes(x = component, y = estimate, colour = package, shape = package)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "Estimated variance", colour = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "spain_cross_package_results.png"), p1_spain / p2_spain,
       width = 8.5, height = 6.4, dpi = 220, bg = "white")

message("Wrote result figures to: ", out_dir)
