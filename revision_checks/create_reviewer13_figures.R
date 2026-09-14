#!/usr/bin/env Rscript

# Creates the requested replacement Figure 3 and the primary-cleaned global
# supplementary figure. It reads saved results only and does not refit models.

suppressPackageStartupMessages({
  library(posterior); library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(metafor); library(orchaRd)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript revision_checks/create_reviewer13_figures.R", call. = FALSE)
figure_args <- commandArgs(trailingOnly = TRUE)
if (length(figure_args) > 1L || (length(figure_args) == 1L && !figure_args[[1L]] %in% c("all", "main"))) {
  stop("Usage: Rscript revision_checks/create_reviewer13_figures.R [all|main]", call. = FALSE)
}
write_supplementary <- !length(figure_args) || figure_args[[1L]] == "all"
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."), mustWork = TRUE)
audit_dir <- file.path(root, "revision_checks", "reviewer13_prediction_outputs")
profile_dir <- file.path(root, "revision_checks", "reviewer13_figure3_profile_outputs")
fig_dir <- file.path(root, "figs", "tutorial")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
close_enough <- function(x, y, tolerance=1e-8) isTRUE(all.equal(unname(x), unname(y), tolerance=tolerance))
orchard_ready_fit <- function(fit) {
  # The immutable saved rma.mv objects have no formula attribute. orchaRd needs
  # it for its plot-data method, so annotate an in-memory copy only; no refit
  # or saved-model mutation occurs.
  out <- fit
  out$formula <- ~ 1
  out
}
summaries <- function(x) {
  q <- quantile(x, c(.025, .25, .5, .75, .975), names = FALSE)
  tibble(lo95=q[1], lo50=q[2], mid=q[3], hi50=q[4], hi95=q[5])
}
density_plot <- function(dat, title, xlab) {
  ss <- dat |> group_by(parameter) |> summarise(
    summaries(value),
    glyph_y = 0.12 * max(stats::density(value)$y),
    .groups="drop"
  )
  ggplot(dat, aes(value)) + geom_density(fill="#73A9AD", colour="#24535A", alpha=.82) +
    geom_segment(data=ss, aes(x=lo95,xend=hi95,y=glyph_y,yend=glyph_y), inherit.aes=FALSE, linewidth=.9, colour="#163B40") +
    geom_segment(data=ss, aes(x=lo50,xend=hi50,y=glyph_y,yend=glyph_y), inherit.aes=FALSE, linewidth=3.0, colour="#163B40") +
    geom_point(data=ss, aes(x=mid,y=glyph_y), inherit.aes=FALSE, size=2.8, shape=21, stroke=.55, fill="white", colour="#163B40") +
    facet_wrap(~parameter, scales="free", ncol=1) +
    labs(title=title, x=xlab, y="Posterior density") +
    theme_classic(base_size=10.5) + theme(plot.title=element_text(face="bold"), strip.background=element_blank())
}

# A. Moura phylogenetic BM model.
mf_path <- file.path(audit_dir, "frequentist_latent_prediction_interval.csv")
latent_path <- file.path(audit_dir, "brms_latent_prediction_draws.csv")
assert(file.exists(mf_path) && file.exists(latent_path), "Run reviewer13_prediction_audit.R first.")
mf <- read_csv(mf_path, show_col_types=FALSE); latent <- read_csv(latent_path, show_col_types=FALSE)
assert(nrow(mf)==1L && nrow(latent)>1000L, "Unexpected Reviewer 13 audit output.")
moura_profile <- read_csv(file.path(profile_dir, "figure3_profile_likelihood_summary.csv"), show_col_types=FALSE) |>
  filter(grepl("^moura_", target))
assert(nrow(moura_profile)==4L && all(moura_profile$finite_two_sided) &&
         all(is.finite(moura_profile$ci_lb)) && all(is.finite(moura_profile$ci_ub)),
       "The full-precision Moura profile-likelihood audit must contain four finite two-sided intervals.")
moura_mf_fit <- readRDS(file.path(root, "revision_checks", "ou_correctness_outputs", "baseline_fit_objects.rds"))$fit_bm
assert(inherits(moura_mf_fit, "rma.mv") && moura_mf_fit$k == 1828L, "Unexpected final Moura BM metafor fit.")
moura_orchard <- orchaRd::mod_results(orchard_ready_fit(moura_mf_fit), mod="1", group="study.id")
moura_orchard_summary <- moura_orchard$mod_table[1L, ]
assert(close_enough(moura_orchard_summary$estimate, mf$estimate) &&
         close_enough(moura_orchard_summary$lowerCL, mf$confidence_interval_lower) &&
         close_enough(moura_orchard_summary$upperCL, mf$confidence_interval_upper) &&
         close_enough(moura_orchard_summary$lowerPR, mf$prediction_interval_lower) &&
         close_enough(moura_orchard_summary$upperPR, mf$prediction_interval_upper),
       "orchaRd Moura summary or prediction interval does not match the validated Reviewer 13 target.")
cat(sprintf("MOURA_ORCHARD_PI_VERIFIED %.16f %.16f\n", moura_orchard_summary$lowerPR, moura_orchard_summary$upperPR))
p_a_mf_mean <- orchaRd::orchard_plot(moura_orchard, xlab="Fisher's Z", k=TRUE, g=FALSE,
                                      k.pos="right", k.size=3.3, legend.pos="bottom.out",
                                      twig.size=.55, branch.size=1.2, trunk.size=.55) +
  labs(title="A. Moura BM — metafor", subtitle="Thick: pooled-mean 95% CI; thin: Reviewer 13 latent-effect 95% PI") +
  theme(plot.title=element_text(face="bold"), plot.subtitle=element_text(size=8.5))
moura_order <- c("Study variance", "Effect-size variance", "Species variance, non-phylogenetic", "Species variance, phylogenetic")
moura_profile <- moura_profile |>
  mutate(component=factor(component, levels=rev(moura_order)))
p_a_mf_var <- ggplot(moura_profile, aes(y=component,x=estimate)) +
  geom_errorbarh(aes(xmin=ci_lb,xmax=ci_ub),height=.18,linewidth=1.05,colour="#528B8B") +
  geom_point(size=3.1,colour="#528B8B") +
  scale_x_continuous(expand=expansion(mult=c(.02,.16))) +
  labs(title="Fitted variances",subtitle="Final-model profile-likelihood 95% CI",x="Variance",y=NULL) +
  theme_classic(base_size=10.5) + theme(plot.title=element_text(face="bold"))
p_a_brms_mean <- density_plot(bind_rows(
  tibble(value=latent$pooled_mean,parameter="Population mean"),
  tibble(value=latent$latent_true_effect,parameter="New latent true effect")
) |> mutate(parameter=factor(parameter, levels=c("Population mean", "New latent true effect"))),"brms","Fisher's Z")
brms_moura <- readRDS(file.path(root,"Rdata","tutorial_v2","moura2021_BM_brms.rds"))
assert(inherits(brms_moura,"brmsfit"),"Moura brms RDS unavailable.")
draws <- as_draws_df(brms_moura)
p_a_brms_var <- density_plot(tibble(
  value=c(draws$sd_study.id__Intercept^2,draws$sd_effect.size.id__Intercept^2,draws$sd_species.id__Intercept^2,draws$sd_species.id.phy__Intercept^2),
  parameter=factor(rep(moura_order,each=nrow(draws)), levels=moura_order)),"brms variance components","Variance")
panel_a <- (p_a_mf_mean / p_a_mf_var) | (p_a_brms_mean / p_a_brms_var)

# B. Validated Spain regional spatial-only comparison.
spain_dir <- file.path(root,"revision_checks","regional_cross_package_audit_outputs")
spain_mf <- read_csv(file.path(spain_dir,"metafor_spatial_only_result.csv"),show_col_types=FALSE)
spain_brms <- read_csv(file.path(spain_dir,"brms_output","brms_parameter_draws.csv"),show_col_types=FALSE)
assert(nrow(spain_mf)==1L && all(c("mean","ci_lb","ci_ub","iid_effect_variance","spatial_variance","rho_km") %in% names(spain_mf)),"Unexpected Spain metafor output.")
spain_tau_raw <- read_csv(file.path(spain_dir,"metafor_profile_ci_tau2.csv"),show_col_types=FALSE)
spain_rho_raw <- read_csv(file.path(spain_dir,"metafor_profile_ci_rho.csv"),show_col_types=FALSE)
spain_tau_profile <- spain_tau_raw |>
  filter(.data[[names(spain_tau_raw)[[1L]]]] == "tau^2") |>
  transmute(ci_lb=.data[["ci.lb"]],ci_ub=.data[["ci.ub"]])
spain_rho_profile <- spain_rho_raw |>
  filter(.data[[names(spain_rho_raw)[[1L]]]] == "rho") |>
  transmute(ci_lb=.data[["ci.lb"]],ci_ub=.data[["ci.ub"]])
assert(nrow(spain_tau_profile)==1L && nrow(spain_rho_profile)==1L &&
       all(is.finite(unlist(spain_tau_profile))) && all(is.finite(unlist(spain_rho_profile))),
       "Saved Spain spatial-variance and range profiles must be finite.")
spain_mf_fit <- readRDS(file.path(spain_dir, "metafor_spatial_only.rds"))
assert(inherits(spain_mf_fit, "rma.mv") && spain_mf_fit$k == 186L, "Unexpected Spain spatial-only metafor fit.")
spain_orchard <- orchaRd::mod_results(orchard_ready_fit(spain_mf_fit), mod="1", group="study_id")
spain_orchard_summary <- spain_orchard$mod_table[1L, ]
assert(close_enough(spain_orchard_summary$estimate, spain_mf$mean) &&
         close_enough(spain_orchard_summary$lowerCL, spain_mf$ci_lb) &&
         close_enough(spain_orchard_summary$upperCL, spain_mf$ci_ub),
       "orchaRd Spain mean or confidence interval does not match the validated spatial-only result.")
p_b_mf_mean <- orchaRd::orchard_plot(spain_orchard, xlab="Hedges' d", k=TRUE, g=FALSE,
                                      k.pos="right", k.size=3.3, legend.pos="bottom.out",
                                      twig.size=0, branch.size=1.2, trunk.size=.55) +
  labs(title="B. Spain regional spatial-only — metafor", subtitle="Pooled mean 95% CI; no spatial prediction interval") +
  theme(plot.title=element_text(face="bold"), plot.subtitle=element_text(size=8.5))
p_b_mf_variance_data <- tibble(
  item=factor(c("IID effect-size variance", "Spatial variance"), levels=c("Spatial variance", "IID effect-size variance")),
  value=c(spain_mf$iid_effect_variance,spain_mf$spatial_variance),
  ci_lb=c(NA_real_,spain_tau_profile$ci_lb), ci_ub=c(NA_real_,spain_tau_profile$ci_ub),
  interval_status=c("Point only: no saved or recoverable profile interval", "Constrained profile-likelihood 95% CI")
)
p_b_mf_variance <- ggplot(p_b_mf_variance_data,aes(y=item,x=value)) +
  geom_errorbarh(data=filter(p_b_mf_variance_data,is.finite(ci_lb)),aes(xmin=ci_lb,xmax=ci_ub),height=.18,linewidth=1.05,colour="#528B8B") +
  geom_point(size=3.1,colour="#528B8B") +
  scale_x_continuous(expand=expansion(mult=c(.02,.18))) +
  labs(title="Fitted variances",subtitle="Spatial: profile 95% CI; IID: point only",x="Variance",y=NULL) + theme_classic(base_size=10.5) + theme(plot.title=element_text(face="bold"))
p_b_mf_range_data <- tibble(item="Spatial range",value=spain_mf$rho_km,ci_lb=spain_rho_profile$ci_lb,ci_ub=spain_rho_profile$ci_ub)
p_b_mf_range <- ggplot(p_b_mf_range_data,aes(y=item,x=value)) +
  geom_errorbarh(aes(xmin=ci_lb,xmax=ci_ub),height=.18,linewidth=1.05,colour="#A44A3F") +
  geom_point(size=3.1,colour="#A44A3F") +
  scale_x_continuous(expand=expansion(mult=c(.02,.18))) + labs(title="Spatial range",subtitle="Profile-likelihood 95% CI",x="rho (km)",y=NULL) + theme_classic(base_size=10.5) + theme(plot.title=element_text(face="bold"))
p_b_mf_components <- p_b_mf_variance / p_b_mf_range
spain_order <- c("Pooled mean (Hedges' d)", "IID effect-size variance", "Spatial variance", "Spatial range (km)")
p_b_brms <- density_plot(spain_brms |> transmute(value,parameter=factor(parameter,levels=spain_order)),"brms","Parameter-specific scale")
panel_b <- (p_b_mf_mean / p_b_mf_components) | p_b_brms
main_figure <- panel_a / panel_b + plot_layout(heights=c(1.18,.82))
ggsave(file.path(fig_dir,"figure3_revised_cross_package.png"),main_figure,width=15,height=17,dpi=240,bg="white")
ggsave(file.path(fig_dir,"figure3_revised_cross_package.pdf"),main_figure,width=15,height=17,device=grDevices::cairo_pdf,bg="white")

# Supplementary Figure: primary cleaned Grau-Andres analysis only.
global_dir <- file.path(root,"revision_checks","reviewer18_influential_effects_outputs")
global <- read_csv(file.path(global_dir,"published_cleaned_primary_results.csv"),show_col_types=FALSE) |>
  mutate(model=factor(model,levels=c("unstructured_only","spatial_only","combined"),labels=c("Unstructured-only","Spatial-only","Combined")))
assert(nrow(global)==3L && all(global$n_effects==2355L) && all(global$n_studies==390L) && all(global$n_sites==380L),"Global data are not the primary cleaned dataset.")
cols <- c("Unstructured-only"="#528B8B","Spatial-only"="#CDAD00","Combined"="#A44A3F")
p_g_mean <- ggplot(global,aes(y=model,x=mean,colour=model)) + geom_vline(xintercept=0,linetype="dashed",colour="grey55") + geom_errorbarh(aes(xmin=ci_lb,xmax=ci_ub),height=.18,linewidth=.9) + geom_point(size=3) + scale_colour_manual(values=cols) + labs(title="Pooled Hedges' d with 95% CI",x="Hedges' d",y=NULL) + theme_classic(base_size=10.5) + theme(legend.position="none",plot.title=element_text(face="bold"))
global_profile_audit <- read_csv(file.path(profile_dir, "global_profile_boundary_audit.csv"), show_col_types=FALSE)
assert(nrow(global_profile_audit) == 9L, "Run reviewer13_global_profile_boundary_audit.R first.")
gvar <- global_profile_audit |>
  filter(parameter %in% c("sigma2", "tau2")) |>
  mutate(model=factor(model,levels=c("unstructured_only","spatial_only","combined"),labels=names(cols)),component=factor(component,levels=c("Study variance","Spatial variance","Effect-size variance")))
p_g_var <- ggplot(gvar,aes(y=model,x=estimate,colour=model)) +
  geom_errorbarh(data=filter(gvar, finite_two_sided),aes(xmin=profile_lower,xmax=profile_upper),height=.18,linewidth=1.05) +
  geom_point(size=3.1) + scale_colour_manual(values=cols) + facet_wrap(~component,scales="free_x",ncol=1) +
  labs(title="Variance components",subtitle="Bars: finite two-sided profile-likelihood 95% CI\nCombined spatial variance: point only (upper side unresolved)",x="Variance",y=NULL,colour=NULL) +
  theme_classic(base_size=10.5) + theme(plot.title=element_text(face="bold"),legend.position="bottom",strip.background=element_blank())
rho_profile <- global_profile_audit |>
  filter(parameter == "rho") |>
  mutate(model=factor(model,levels=c("unstructured_only","spatial_only","combined"),labels=names(cols)))
assert(nrow(rho_profile) == 2L && !any(rho_profile$finite_two_sided) &&
         sum(rho_profile$lower_is_search_boundary) == 2L && sum(rho_profile$upper_is_search_boundary) == 1L,
       "Global rho profiles must retain their boundary-status distinction.")
p_g_range <- ggplot(rho_profile,aes(y=model,x=estimate,colour=model)) +
  geom_point(size=3) + scale_colour_manual(values=cols) + scale_x_continuous(expand=expansion(mult=c(.02,.18))) +
  labs(title="Spatial ranges",subtitle="Point estimates only: spatial-only lower side unresolved; finite upper bound\nCombined model: both sides unresolved over searched range",x="Range rho (km)",y=NULL) +
  theme_classic(base_size=10.5) + theme(legend.position="none",plot.title=element_text(face="bold"))
supplementary <- p_g_mean / (p_g_var | p_g_range)
if (write_supplementary) {
  ggsave(file.path(fig_dir,"supplementary_cleaned_global_spatial_sensitivity.png"),supplementary,width=15,height=10.5,dpi=240,bg="white")
  ggsave(file.path(fig_dir,"supplementary_cleaned_global_spatial_sensitivity.pdf"),supplementary,width=15,height=10.5,device=grDevices::cairo_pdf,bg="white")
}
cat("REVIEWER13_FIGURE3_CHECKS_PASSED\n")
