#!/usr/bin/env Rscript

# Validate the affected tutorial examples in a fresh R session without sampling
# from, or optimising, any statistical model.

suppressPackageStartupMessages({
  library(ape)
  library(brms)
  library(metadat)
  library(metafor)
})

source_lines <- readLines("tutorial_v2.qmd", warn = FALSE)
source_text <- paste(source_lines, collapse = "\n")

stopifnot(
  !grepl("A = lim_vcv", source_text, fixed = TRUE),
  !grepl("phylo_eg2_tmb$fit$par", source_text, fixed = TRUE),
  grepl("phylo_eg2_tmb_ma$fit$par", source_text, fixed = TRUE)
)

dat_lim2014 <- metafor::escalc(
  measure = "ZCOR", ri = ri, ni = ni, data = metadat::dat.lim2014$o_o_unadj
)
dat_lim2014$phy <- dat_lim2014$species
dat_lim2014$id <- seq_len(nrow(dat_lim2014))

tre_lim2014 <- ape::compute.brlen(metadat::dat.lim2014$o_o_unadj_tree)
A <- ape::vcv.phylo(tre_lim2014, corr = TRUE)
V <- diag(dat_lim2014$vi)
rownames(V) <- colnames(V) <- dat_lim2014$id

lim_vcv_formula <- brms::bf(
  yi ~ 1 + (1 | species) + (1 | gr(phy, cov = A)) +
    (1 | gr(id, cov = V))
)
lim_vcv_priors <- c(
  brms::set_prior("normal(0, 1)", class = "Intercept"),
  brms::set_prior("exponential(1)", class = "sigma"),
  brms::set_prior("exponential(1)", class = "sd", group = "species"),
  brms::set_prior("exponential(1)", class = "sd", group = "phy"),
  brms::set_prior("constant(1)", class = "sd", group = "id")
)
vcv_stancode <- brms::make_stancode(
  lim_vcv_formula, data = dat_lim2014, data2 = list(A = A, V = V),
  prior = lim_vcv_priors, family = gaussian()
)

lim_se_formula <- brms::bf(
  yi | se(sqrt(vi), sigma = TRUE) ~ 1 + (1 | species) +
    (1 | gr(phy, cov = A))
)
lim_se_priors <- c(
  brms::set_prior("normal(0, 1)", class = "Intercept"),
  brms::set_prior("exponential(1)", class = "sigma"),
  brms::set_prior("exponential(1)", class = "sd", group = "species"),
  brms::set_prior("exponential(1)", class = "sd", group = "phy")
)
se_stancode <- brms::make_stancode(
  lim_se_formula, data = dat_lim2014, data2 = list(A = A),
  prior = lim_se_priors, family = gaussian()
)

dat_lim2014$environment <- factor(dat_lim2014$environment)
lim_mr_formula <- brms::bf(
  yi ~ 1 + environment + (1 | species) + (1 | gr(phy, cov = A)) +
    (1 | gr(id, cov = V))
)
lim_mr_priors <- c(
  brms::set_prior("normal(0, 1)", class = "Intercept"),
  brms::set_prior("normal(0, 1)", class = "b"),
  brms::set_prior("exponential(1)", class = "sigma"),
  brms::set_prior("exponential(1)", class = "sd", group = "species"),
  brms::set_prior("exponential(1)", class = "sd", group = "phy"),
  brms::set_prior("constant(1)", class = "sd", group = "id")
)
mr_stancode <- brms::make_stancode(
  lim_mr_formula, data = dat_lim2014, data2 = list(A = A, V = V),
  prior = lim_mr_priors, family = gaussian()
)

stopifnot(
  nchar(vcv_stancode) > 0L,
  nchar(se_stancode) > 0L,
  nchar(mr_stancode) > 0L
)

cat("TUTORIAL_VISIBLE_OBJECTS_PASSED\n")
