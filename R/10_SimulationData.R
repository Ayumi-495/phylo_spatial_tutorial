set.seed(357)

# True parameters (meta-analysis scale: Fisher's Z)
k.studies    <- 300
k.per.study  <- 20
k            <- k.studies * k.per.study

mu_z         <- atanh(0.20)
sigma2.u     <- 0.20     # between-study (non-spatial)
sigma2.s     <- 0.30     # within-study extra (non-spatial)
tau2         <- 0.40     # spatial variance
rho_within   <- 0.50     # within-study correlation in sampling error

# Spatial domain (study-level coords)
sd_x_km <- 800; sd_y_km <- 600
loc <- data.frame(
  study = seq_len(k.studies),
  x_km  = rnorm(k.studies, 0, sd_x_km),
  y_km  = rnorm(k.studies, 0, sd_y_km)
)

D_study  <- as.matrix(dist(loc[, c("x_km","y_km")]))
d_half   <- 800
rho_dec  <- log(2) / d_half
R_study  <- exp(-rho_dec * D_study)

# Spatial random effect: s_study ~ N(0, tau2 * R)
eps <- 1e-8
Sigma_sp <- tau2 * R_study + diag(eps, k.studies)
L <- chol(Sigma_sp)
s_study <- as.vector(t(L) %*% rnorm(k.studies))

# IDs
study <- rep(loc$study, each = k.per.study)
esid  <- ave(study, study, FUN = seq_along)
id    <- seq_len(k)

# Non-spatial random effects
u_study <- rnorm(k.studies, 0, sqrt(sigma2.u))
u_u <- u_study[study] # between-study mapped to effects
u_s <- rnorm(k, 0, sqrt(sigma2.s)) # within-study extra per effect

# Sampling sizes (lognormal) & variances on Z 
n_meanlog <- 3.7 # mean(log n)
n_sdlog <- 1.0 # sd(log n)

n_ij <- pmax(10L, round(rlnorm(k, meanlog = n_meanlog, sdlog = n_sdlog)))
vi   <- 1 / pmax(7, n_ij - 3) # Var(Z) ≈ 1/(n-3), guard small n

# Sampling error VCV: block-diagonal by study with rho_within
VCV <- matrix(0, nrow = k, ncol = k)
for (st in unique(study)) {
  idx <- which(study == st)
  S <- matrix(rho_within, nrow = length(idx), ncol = length(idx))
  diag(S) <- 1
  S <- S * sqrt(outer(vi[idx], vi[idx]))   # scale by sqrt(vi_i * vi_j)
  VCV[idx, idx] <- S
}
VCV[upper.tri(VCV)] <- t(VCV)[upper.tri(VCV)]
diag(VCV) <- vi

# Sampling error draw: mi ~ N(0, VCV)
mi <- MASS::mvrnorm(n = 1, mu = rep(0, k), Sigma = VCV)

# Observed effects (Z)
s  <- s_study[study]
yi <- mu_z + s + u_u + u_s + mi

# Combine 
dat <- data.frame(
  study = factor(study),
  esid  = esid,
  id = factor(id),
  yi = yi,
  vi = vi,
  x_km = loc$x_km[study],
  y_km  = loc$y_km[study],
  n = n_ij
)

head(dat, 30)
mean(dat$yi)
var(dat$yi)
mean(dat_vars$yi)
var(dat_means$yi)
dat_means <- aggregate(yi ~ study, data = dat, mean) 
var(dat_means$yi)
mean(dat$vi)

dat$const <- 1

system.time(
  fit <- rma.mv(
    yi = yi, 
    V = vi,
    random = list(
      ~ 1 | study,
      ~ 1 | id, 
    ~ x_km + y_km |const
    ), 
    struct = "SPEXP", 
    data = dat,
    sparse = TRUE,
    verbose = TRUE
    )
)
# user   system  elapsed 
# 12384.22   262.65 20069.32 
summary(fit)
# Multivariate Meta-Analysis Model (k = 6000; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc   
# -5610.7329  11221.4657  11231.4657  11264.9625  11231.4757   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed  factor 
# sigma^2.1  0.2205  0.4696    300     no   study 
# sigma^2.2  0.2879  0.5366   6000     no      id 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 300)
# 
# estim    sqrt  fixed 
# tau^2        0.2762  0.5255     no 
# rho        721.7237             no 
# 
# Test for Heterogeneity:
#   Q(df = 5999) = 280705.5665, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub     
# 0.5723  0.2083  2.7473  0.0060  0.1640  0.9806  ** 


dat$id <- factor(dat$id)
VCV <- diag(dat$vi, nrow = nrow(dat)) 
rownames(VCV)<- colnames(VCV)<- dat$id
VCV[1:5, 1:5]

dat$pos <- numFactor(dat$x_km, dat$y_km)

system.time(
  tmb <- glmmTMB(yi ~ 1 +
    equalto(0 + id | const, VCV) + 
    (1 | study) +
    exp(pos + 0 | const),     
  data = dat,
  REML = TRUE
  )
)
# user  system elapsed 
#812.421  11.121 834.193 

head(confint(tmb), 10)
#                                          2.5 %    97.5 %  Estimate
# (Intercept)                                              0.1527657 0.9918287 0.5722972
# Std.Dev.(Intercept)|study                                0.4155384 0.5306493 0.4695798
# Std.Dev.pos(-551.252000853363,-1370.3045899935)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(-218.706116577143,-1328.96626646076)|const.1 0.3433570 0.8042729 0.5255024
# Std.Dev.pos(378.931902660787,-1285.05806598531)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(-1073.8917378376,-1280.66200353426)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(44.5099736530509,-1220.60786402491)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(34.0059462884198,-1174.32729886188)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(104.067536928893,-1159.07462992095)|const.1  0.3433570 0.8042729 0.5255024
# Std.Dev.pos(-777.775431059205,-1136.9892610501)|const.1  0.3433570 0.8042729 0.5255024
# > 

sigma(tmb)^2
# 0.05065844
tmb_varcor <- VarCorr(tmb)$cond
# head(tmb_2_varcor, 10)

tmb$fit$par[[4]]



metafor_1 <- data.frame(model = "metafor", 
                        logLik = logLik(fit),
                        est = fit$b[[1]], 
                        se = fit$se[[1]], 
                        sigma2.u = fit$sigma2[1],  # among study variance
                        sigma2.m = fit$sigma2[2],  # within study variance
                        tau2 = fit$tau2,
                        rho = fit$rho)


glmmTMB_1 <- data.frame(model = "glmmTMB",
                        logLik = logLik(tmb)[1],
                        est = unlist(fixef(tmb))[[1]], #overall mean
                        se = as.numeric(sqrt(vcov(tmb)[[1]])), #overall mean SE
                        sigma2.u = (tmb_varcor$study_id[1]),  #among study variance estimate
                        sigma2.m = sigma(tmb)^2,  #within study variance estimate
                        tau2 = (tmb_varcor$const.1[1]),  #spatial variance estimate
                        rho = exp(tmb$fit$par[[4]])  # rho
)

output <- rbind(metafor_1, glmmTMB_1)
knitr::kable(output)
# |model   |    logLik|       est|        se|  sigma2.u|  sigma2.m|      tau2|      rho|
#   |:-------|---------:|---------:|---------:|---------:|---------:|---------:|--------:|
#   |metafor | -1659.075| 0.2742624| 0.2338043| 0.3048980| 0.2250742| 0.1622631| 2792.906|
#   |glmmTMB | -1663.334| 0.2742563| 0.2641737| 0.0929624| 0.0506584| 0.1622714| 2793.032|


 ## ----  # #
coords <- data.frame(
  study_id = factor(seq_len(n_studies)),
  x_km = x_km, y_km = y_km, lat = lat, lon = lon
)


enter <- st_sfc(st_point(c(10, 50)), crs = 4326) # 4326 = WGS84
center_proj <- st_transform(center, 3035)


x_m <- rnorm(n_studies, 0, 1200000) # 1200 km = 1,200,000 m
y_m <- rnorm(n_studies, 0, 900000)  # 900 km = 900,000 m
coords_proj <- st_coordinates(center_proj)[,1:2]
pts_proj <- cbind(coords_proj[1] + x_m, coords_proj[2] + y_m)

pts <- st_as_sf(data.frame(id = 1:n_studies),
                coords = pts_proj, crs = 3035)

pts_lonlat <- st_transform(pts, 4326)

head(st_coordinates(pts_lonlat))


