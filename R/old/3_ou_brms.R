library(brms)

lev <- levels(HubDen_2004$Phylo)
D <- as.matrix(D_HD)[lev, lev]

mds <- cmdscale(D, k = 10, eig = TRUE, add = TRUE)
coords <- as.data.frame(mds$points)
names(coords) <- paste0("Z", 1:ncol(coords))
coords$Phylo <- rownames(coords)


dat_gp <- merge(HubDen_2004, coords, by = "Phylo", all.x = TRUE)

test_mod <- bf(d | se(sqrt(var)) ~ 1 +
                   (1 | Reference) + (1 | effect_id) + (1 | Herbivore) +
                   gp(Z1, Z2, 
                      Z3, Z4, 
                      Z5, Z6, 
                      Z7, Z8, 
                      Z9, Z10,
                      cov = "exponential", gr = TRUE, scale = FALSE))

pri <- get_prior(test_mod,
                     data = dat_gp,
                     family = gaussian())


max_cores <- 10
num_chains <- 2
threads_per_chain <- floor(max_cores / num_chains)


fit_exp2 <- brm(test_mod,
  data = dat_gp,
  family = gaussian(),
  prior = pri,
  iter = 10000, 
  warmup = 7000, 
  chains = num_chains,
  backend = "cmdstanr",
  threads = threading(threads_per_chain),
  control = list(adapt_delta = 0.95)
  )

summary(fit_exp2)
summary(HD_spatial)
fit_exp_2alpha <- 1/(0.08)
HD_OU_brms <- exp(-fit_exp_2alpha * D_HD)



formula_dat1_OU <- bf(d | se(sqrt(var)) ~ 1 + 
                        (1 | Reference) + 
                        (1 | effect_id) + 
                        (1 | Herbivore) +  # non-phylogenetic random effect
                        (1 |a| gr(Phylo, cov = HD_OU_brms))
)

prior_dat1_OU <- get_prior(
  formula = formula_dat1_OU,
  data = HubDen_2004,
  data2 = list(HD_OU_brms = HD_OU_brms),
  family = gaussian()
)

OU_brms1 <- brm(
  formula = formula_dat1_OU,
  family = gaussian(),
  data = HubDen_2004,
  data2 = list(HD_OU_brms = HD_OU_brms),
  prior = prior_dat1_OU,
  iter = 8000, 
  warmup = 6000,  
  chains = 2,
  backend = "cmdstanr",
  threads = threading(threads_per_chain),
  control = list(adapt_delta = 0.95)
)

summary(OU_brms1)
summary(HD_OU)

#####----#####
D_org <- as.matrix(D_HD)[coords$Phylo, coords$Phylo, drop = FALSE]

Zmat  <- as.matrix(coords[ , grep("^Z", names(coords)), drop = FALSE])
D_rec <- as.matrix(dist(Zmat, method = "euclidean"))


D_org[1:10, 1:10]
D_rec[1:10, 1:10]

ix <- upper.tri(D_ord, diag = FALSE)

pearson  <- cor(D_ord[ix], D_rec[ix], method = "pearson")
spearman <- cor(D_ord[ix], D_rec[ix], method = "spearman")
rmse     <- sqrt(mean( (D_ord[ix] - D_rec[ix])^2 ))
stress1  <- sqrt( sum((D_ord[ix] - D_rec[ix])^2) / sum(D_ord[ix]^2) )

pearson; spearman; rmse; stress1
# [1] 0.9930093
# [1] 0.9303287
# [1] 0.06352045
# [1] 0.0837915

plot(D_ord[ix], D_rec[ix], pch = 16, cex = 0.6,
     xlab = "Original distance", ylab = "Reconstructed (MDS) distance")
abline(0, 1, lty = 2)
