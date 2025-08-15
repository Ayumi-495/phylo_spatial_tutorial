library(ggplot2)


# Figure 1 ----
# grid of (normalized) distances
d <- seq(0, 1, length.out = 501)

# parameters
rho <- 0.2          # spatial range parameter
alpha <- 1 / rho    # OU–exponential equivalence: alpha = 1/rho
bm_floor <- 1e-3    # tiny floor so BM-like curve doesn't hit zero (illustrative)

# kernel functions (illustrative BM; exact BM depends on shared path lengths on a tree)
k_linear      <- pmax(0, 1 - d / rho)                   # spatial linear with hard cutoff at d = rho
k_exponential <- exp(-d / rho)                          # spatial exponential
k_squaredexp  <- exp(-(d^2) / (rho^2))                  # spatial squared exponential (Gaussian)
k_ou          <- exp(-alpha * d)                        # phylogenetic OU; identical to exponential with alpha = 1/rho
k_bm_like     <- pmax(bm_floor, 1 - d)                  # BM-like linear decay without hard cutoff (illustrative)

# assemble tidy data
df <- tibble(
  distance = rep(d, 5),
  correlation = c(k_bm_like, k_ou, k_linear, k_exponential, k_squaredexp),
  model = factor(rep(c("BM",
                       "OU (α = 1/ρ)",
                       "Linear (spatial)",
                       "Exponential (spatial)",
                       "Squared exponential (spatial)"),
                     each = length(d)),
                 levels = c("BM",
                            "OU (α = 1/ρ)",
                            "Linear (spatial)",
                            "Exponential (spatial)",
                            "Squared exponential (spatial)"))
)

# plot
ggplot(df, aes(distance, correlation, color = model)) +
  geom_line(size = 1) +
  labs(x = "Distance (normalised)", y = "Correlation",
       color = NULL,
       title = "Comparison of correlation–distance relationships") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
