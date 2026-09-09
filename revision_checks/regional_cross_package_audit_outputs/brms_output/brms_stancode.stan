// generated with brms 2.22.0
functions {
  /* compute a latent Gaussian process with exponential kernel
   * Args:
   *   x: array of continuous predictor values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   *   zgp: vector of independent standard normal variables
   * Returns:
   *   a vector to be added to the linear predictor
   */
  vector gp_exponential(data array[] vector x, real sdgp, vector lscale, vector zgp) {
    int Dls = rows(lscale);
    int N = size(x);
    matrix[N, N] cov;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      cov = gp_exponential_cov(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      cov = gp_exponential_cov(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        cov = cov .* gp_exponential_cov(x[, d], 1, lscale[d]);
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-12;
    }
    return cholesky_decompose(cov) * zgp;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  // number of latent GP groups
  int<lower=1> Nsubgp_1;
  // indices of latent GP groups per observation
  array[N] int<lower=1> Jgp_1;
  array[Nsubgp_1] vector[Dgp_1] Xgp_1;  // covariates of the GP
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  array[Kgp_1] vector<lower=0>[1] lscale_1;  // GP length-scale parameters
  vector[Nsubgp_1] zgp_1;  // latent variables of the GP
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(Intercept | 0, 1);
  lprior += student_t_lpdf(sdgp_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  lprior += lognormal_lpdf(lscale_1[1][1] | log(50), 1);
  lprior += student_t_lpdf(sigma | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[Nsubgp_1] gp_pred_1 = gp_exponential(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + gp_pred_1[Jgp_1];
    target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2));
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zgp_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
