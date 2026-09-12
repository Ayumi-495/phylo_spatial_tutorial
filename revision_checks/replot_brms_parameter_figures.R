#!/usr/bin/env Rscript

# Replot final Bayesian posterior parameter figures from saved draw tables only.
# This presentation pass does not read or modify an RDS, refit a model, or alter
# posterior draws, intervals, PPCs, or diagnostics.
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidybayes)
  library(patchwork)
})

figures <- c(
  "revision_checks/moura_bm_brms_precompute_outputs/moura_bm_brms_parameter_draws.csv" =
    "revision_checks/moura_bm_brms_precompute_outputs/moura_bm_brms_parameter_distributions.png",
  "revision_checks/brms_consistency_outputs/moura_mr/posterior_parameter_draws.csv" =
    "revision_checks/brms_consistency_outputs/moura_mr/posterior_parameter_distributions.png",
  "revision_checks/brms_consistency_outputs/lim_vcv_ma_adapt_delta_0.99/posterior_parameter_draws.csv" =
    "revision_checks/brms_consistency_outputs/lim_vcv_ma_adapt_delta_0.99/posterior_parameter_distributions.png",
  "revision_checks/brms_consistency_outputs/lim_se_ma_adapt_delta_0.99/posterior_parameter_draws.csv" =
    "revision_checks/brms_consistency_outputs/lim_se_ma_adapt_delta_0.99/posterior_parameter_distributions.png",
  "revision_checks/brms_consistency_outputs/lim_mr_adapt_delta_0.99/posterior_parameter_draws.csv" =
    "revision_checks/brms_consistency_outputs/lim_mr_adapt_delta_0.99/posterior_parameter_distributions.png",
  "revision_checks/regional_cross_package_audit_outputs/brms_output/brms_parameter_draws.csv" =
    "revision_checks/regional_cross_package_audit_outputs/brms_output/brms_parameter_distributions.png"
)

make_panel <- function(data, label, caption = NULL) {
  data$parameter <- factor(data$parameter, levels = rev(unique(data$parameter)))
  ggplot(data, aes(x = value, y = parameter)) +
    tidybayes::stat_halfeye(.width = c(0.5, 0.95), point_interval = "median_qi",
                             fill = "#73A9AD", color = "#24535A") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
    labs(title = label,
         x = "Posterior distribution (median, 50% and 95% credible intervals)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

for (input in names(figures)) {
  draws <- read.csv(input, stringsAsFactors = FALSE)
  stopifnot(all(c("value", "parameter", "panel") %in% names(draws)))

  # Fixed effects and spatial range retain their within-group axes.  Each
  # variance component is its own one-row small multiple so a wide component
  # cannot compress the others; all values remain on the original variance
  # scale and each row includes its own full 50% and 95% credible intervals.
  non_variance <- split(draws[draws$panel != "Variance", ],
                        draws$panel[draws$panel != "Variance"])
  non_variance_panels <- Map(
    function(data, label) make_panel(data, label),
    non_variance, names(non_variance)
  )
  variance <- draws[draws$panel == "Variance", ]
  variance_parameters <- unique(variance$parameter)
  variance_panels <- lapply(variance_parameters, function(parameter) {
    data <- variance[variance$parameter == parameter, ]
    make_panel(data, paste("Variance:", unique(data$parameter)))
  })
  panels <- c(non_variance_panels, variance_panels)
  figure <- patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      caption = "Each variance row uses its own horizontal scale; compare numerical values and intervals, not horizontal widths, across variance rows."
    )
  ggsave(figures[[input]], figure, width = 10,
         height = max(4.0, 2.15 * length(panels) + 0.4), dpi = 180)
}
