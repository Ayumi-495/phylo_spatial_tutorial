#!/usr/bin/env Rscript

# Create public tutorial figures from the interval artifacts produced by
# precompute_visualisation_intervals.R. This script does not fit, profile, or
# modify any statistical model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path.", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
out_dir <- file.path(root, "figs", "tutorial")
interval_dir <- file.path(root, "revision_checks", "visualisation_interval_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
read_interval <- function(name) {
  path <- file.path(interval_dir, name)
  assert(file.exists(path), paste("Missing interval artifact:", path))
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("analysis", "package", "parameter", "estimate", "ci_lb", "ci_ub", "interval_method")
  assert(all(required %in% names(dat)), paste("Malformed interval artifact:", path))
  assert(all(is.finite(as.matrix(dat[c("estimate", "ci_lb", "ci_ub")]))) &&
           all(dat$ci_lb <= dat$estimate) && all(dat$estimate <= dat$ci_ub),
         paste("Invalid interval ordering:", path))
  dat
}

estimate_panel <- function(dat, title, x_label, zero = FALSE, colour = "#528B8B") {
  p <- ggplot(dat, aes(x = estimate, y = parameter)) +
    geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.16,
                  orientation = "y", colour = colour, linewidth = 0.75) +
    geom_point(size = 3, colour = colour) +
    labs(title = title, x = x_label, y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
  if (zero) p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45")
  p
}

save_package_figure <- function(intervals, package, fixed_parameters, component_parameters,
                                title, filename, range_parameter = NULL) {
  fixed <- intervals[intervals$package == package & intervals$parameter %in% fixed_parameters, , drop = FALSE]
  components <- intervals[intervals$package == package & intervals$parameter %in% component_parameters, , drop = FALSE]
  assert(nrow(fixed) == length(fixed_parameters), paste("Missing fixed-effect intervals for", title))
  assert(nrow(components) == length(component_parameters), paste("Missing variance intervals for", title))
  fixed$parameter <- factor(fixed$parameter, levels = rev(fixed_parameters))
  components$parameter <- factor(components$parameter, levels = rev(component_parameters))
  p_fixed <- estimate_panel(fixed, title, "Effect size (95% confidence interval)", zero = TRUE, colour = "#CDAD00")
  p_var <- estimate_panel(components, NULL, "Variance (95% confidence interval)")
  panels <- p_fixed / p_var
  if (!is.null(range_parameter)) {
    range <- intervals[intervals$package == package & intervals$parameter == range_parameter, , drop = FALSE]
    assert(nrow(range) == 1L, paste("Missing range interval for", title))
    range$parameter <- factor(range$parameter, levels = range_parameter)
    panels <- panels / estimate_panel(range, NULL, "Range in km (95% confidence interval)")
  }
  height <- if (is.null(range_parameter)) 6.4 else 8.2
  ggsave(file.path(out_dir, filename), panels, width = 8.5, height = height, dpi = 220, bg = "white")
}

extract_spain_glmmtmb <- function() {
  fit_path <- file.path(root, "revision_checks", "regional_cross_package_audit_outputs", "glmmTMB_spatial_only.rds")
  assert(file.exists(fit_path), paste("Missing saved Spain glmmTMB model:", fit_path))
  fit <- readRDS(fit_path)
  assert(inherits(fit, "glmmTMB") && isTRUE(fit$sdr$pdHess), "Spain glmmTMB fit is not usable for Wald intervals.")
  par <- fit$fit$par
  cov_fixed <- fit$sdr$cov.fixed
  assert(identical(names(par), colnames(cov_fixed)), "Spain glmmTMB parameter order is inconsistent.")
  idx <- c(which(names(par) == "betadisp"), which(names(par) == "theta"))
  assert(length(idx) == 3L, "Expected dispersion plus two spatial theta parameters.")
  se <- sqrt(diag(cov_fixed)[idx])
  component <- c("IID effect-size variance", "Spatial variance", "Exponential range (km)")
  transform <- c("variance", "variance", "range")
  estimate <- exp(par[idx])
  lower <- exp(par[idx] - qnorm(0.975) * se)
  upper <- exp(par[idx] + qnorm(0.975) * se)
  estimate[transform == "variance"] <- estimate[transform == "variance"]^2
  lower[transform == "variance"] <- lower[transform == "variance"]^2
  upper[transform == "variance"] <- upper[transform == "variance"]^2
  beta <- glmmTMB::fixef(fit)$cond
  beta_se <- sqrt(diag(stats::vcov(fit)$cond))[names(beta)]
  fixed <- data.frame(
    analysis = "spain", package = "glmmTMB", parameter = names(beta), estimate = as.numeric(beta),
    ci_lb = as.numeric(beta - qnorm(0.975) * beta_se), ci_ub = as.numeric(beta + qnorm(0.975) * beta_se),
    interval_method = "95% Wald confidence interval", stringsAsFactors = FALSE
  )
  components <- data.frame(
    analysis = "spain", package = "glmmTMB", parameter = component, estimate = estimate,
    ci_lb = lower, ci_ub = upper,
    interval_method = "95% Wald confidence interval on log SD/log range scale", stringsAsFactors = FALSE
  )
  out <- rbind(fixed, components)
  assert(all(out$ci_lb <= out$estimate) && all(out$estimate <= out$ci_ub), "Invalid Spain glmmTMB intervals.")
  write.csv(out, file.path(interval_dir, "spain_glmmTMB_wald_intervals.csv"), row.names = FALSE)
  out
}

read_spain_brms <- function() {
  path <- file.path(root, "revision_checks", "regional_cross_package_audit_outputs", "brms_output", "brms_spatial_only_result.csv")
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  labels <- c(
    pooled_mean = "(Intercept)",
    iid_effect_variance = "IID effect-size variance",
    spatial_variance = "Spatial variance",
    rho_km = "Exponential range (km)"
  )
  out <- data.frame(
    analysis = "spain", package = "brms", parameter = unname(labels[dat$parameter]),
    estimate = dat$estimate, ci_lb = dat$ci_lb, ci_ub = dat$ci_ub,
    interval_method = "95% posterior credible interval", stringsAsFactors = FALSE
  )
  assert(!anyNA(out$parameter) && all(out$ci_lb <= out$estimate) && all(out$estimate <= out$ci_ub),
         "Malformed Spain brms interval artifact.")
  out
}

spain_comparison_panel <- function(dat, parameters, title, x_label, log_scale = FALSE) {
  panel <- dat[dat$parameter %in% parameters, , drop = FALSE]
  assert(nrow(panel) == length(parameters) * 3L, paste("Incomplete Spain comparison panel:", title))
  panel$parameter <- factor(panel$parameter, levels = parameters)
  panel$package <- factor(panel$package, levels = c("metafor", "glmmTMB", "brms"))
  p <- ggplot(panel, aes(x = estimate, y = package, colour = package, shape = package)) +
    geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.14,
                  orientation = "y", linewidth = 0.7) +
    geom_point(size = 2.8) +
    facet_wrap(~ parameter, nrow = 1L, scales = "free_x") +
    labs(title = title, x = x_label, y = NULL, colour = NULL, shape = NULL) +
    scale_colour_manual(
      values = c(metafor = "#E76F6A", glmmTMB = "#619CFF", brms = "#00BA38"),
      breaks = c("metafor", "glmmTMB", "brms"),
      labels = c("metafor (95% CI)", "glmmTMB (95% Wald CI)", "brms (95% CrI)")
    ) +
    scale_shape_manual(
      values = c(metafor = 16, glmmTMB = 15, brms = 17),
      breaks = c("metafor", "glmmTMB", "brms"),
      labels = c("metafor (95% CI)", "glmmTMB (95% Wald CI)", "brms (95% CrI)")
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")
  if (log_scale) p <- p + scale_x_log10()
  p
}

moura <- read_interval("moura_frequentist_intervals.csv")
lim <- read_interval("lim_frequentist_intervals.csv")
spain_metafor <- read_interval("spain_metafor_intervals.csv")
spain_glmmtmb <- extract_spain_glmmtmb()
spain_brms <- read_spain_brms()

save_package_figure(
  moura, "glmmTMB", "(Intercept)",
  c("Study variance", "Effect-size variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic"),
  "Moura et al.: BM meta-analysis (glmmTMB)", "tmb_eg1_1.png"
)
save_package_figure(
  lim, "metafor", c("intrcpt", "environmentwild"),
  c("Effect-size variance", "Species variance, phylogenetic", "Species variance, non-phylogenetic"),
  "Lim et al.: BM meta-regression (metafor)", "metafor_lim.png"
)
save_package_figure(
  lim, "glmmTMB", c("(Intercept)", "environmentwild"),
  c("Effect-size variance", "Species variance, phylogenetic", "Species variance, non-phylogenetic"),
  "Lim et al.: BM meta-regression (glmmTMB)", "lim_glmmtmb_results.png"
)
save_package_figure(
  spain_metafor, "metafor", "intrcpt",
  c("IID effect-size variance", "Spatial variance"),
  "Spain regional subset: spatial-only model (metafor)", "spain_metafor_spatial_only.png",
  range_parameter = "Exponential range (km)"
)
save_package_figure(
  spain_glmmtmb, "glmmTMB", "(Intercept)",
  c("IID effect-size variance", "Spatial variance"),
  "Spain regional subset: spatial-only model (glmmTMB)", "spain_glmmtmb_spatial_only.png",
  range_parameter = "Exponential range (km)"
)

spain_metafor_comparison <- spain_metafor
spain_metafor_comparison$parameter[spain_metafor_comparison$parameter == "intrcpt"] <- "(Intercept)"
spain_all <- rbind(spain_metafor_comparison, spain_glmmtmb, spain_brms)
p_fixed <- spain_comparison_panel(spain_all, "(Intercept)",
                                  "Spain regional subset: spatial-only common target",
                                  "Pooled effect size (95% CI or CrI)") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45")
p_variance <- spain_comparison_panel(spain_all, c("IID effect-size variance", "Spatial variance"),
                                     NULL, "Variance (95% CI or CrI)")
p_range <- spain_comparison_panel(spain_all, "Exponential range (km)", NULL,
                                  "Range in km (95% CI or CrI)", log_scale = TRUE)
spain_comparison <- (p_fixed / p_variance / p_range) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "spain_cross_package_results.png"), spain_comparison,
       width = 10, height = 9.2, dpi = 220, bg = "white")

writeLines(c(
  "Tutorial visualisation figures regenerated from stored interval artifacts.",
  "No statistical model was fitted or profiled by this script.",
  "glmmTMB Spain intervals use the saved model's Wald covariance on transformed parameter scales.",
  "metafor Moura, Lim, and Spain component intervals were supplied by precompute_visualisation_intervals.R.",
  "brms Spain intervals are read from the saved posterior summary."
), file.path(interval_dir, "figure_generation_manifest.txt"))
message("VISUALISATION_INTERVAL_FIGURES_WRITTEN")
