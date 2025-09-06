set.seed(42)
n_study <- 500
beta0_r <- 0.30
theta2 <- 7.5
tau2 <- 0.08
sigma2_study <- 0.10
sigma2_es <- 0.05

sd_x_km    <- 1200
sd_y_km    <-  900
center_lat <- 50
center_lon <- 10

n_meanlog <- 3.7   
n_sdlog   <- 1.0

# study location ----
x_km <- rnorm(n_study, 0, sd_x_km)
y_km <- rnorm(n_study, 0, sd_y_km)

km_per_deg_lat <- 111
km_per_deg_lon <- 111 * cos(center_lat * pi/180)
lat <- center_lat + y_km / km_per_deg_lat
lon <- center_lon + x_km / km_per_deg_lon

coords <- data.frame(
  study_id = factor(seq_len(n_study)),
  x_km = x_km, y_km = y_km, lat = lat, lon = lon
)

# distance matrix ----
D   <- as.matrix(dist(cbind(x_km, y_km)))  
phi <- exp(-theta2)
R   <- exp(-phi * D)
d_half <- log(2) / phi

# spatial random effect ----
eps  <- 1e-10
Sigma_sp <- tau2 * R + diag(eps, n_study)
L <- chol(Sigma_sp)
s <- as.vector(t(L) %*% rnorm(n_study))

# between study random effect -----
a <- rnorm(n_study, 0, sqrt(sigma2_study))

# each study have 10 datapoints ----
m_i <- rep(10, n_study)
beta0_z <- atanh(beta0_r)

dat <- vector("list", n_study)
for(i in seq_len(n_study)){
  k <- m_i[i]
  n_ij <- pmax(10, round(rlnorm(k, n_meanlog, n_sdlog)))
  vi <- 1/pmax(7, n_ij-3)
  u_ij <- rnorm(k, 0, sqrt(sigma2_es))
  
  z_true <- beta0_z + s[i] + a[i] + u_ij
  yi <- rnorm(k, mean = z_true, sd = sqrt(vi))
  r_obs <- tanh(yi)
  
  dat[[i]] <- data.frame(
    study_id = coords$study_id[i],
    x_km = coords$x_km[i],
    y_km = coords$y_km[i],
    lat = coords$lat[i],
    lon = coords$lon[i],
    n = n_ij,
    vi = vi,
    yi_z = yi,
    r = r_obs
  )
}
dat <- do.call(rbind, dat)
dat$es_id <- seq_len(nrow(dat))
dat$const <- 1
# write_csv(dat, here("data", "simulated_v1.csv"))
View(dat)

system.time(
  fit <- rma.mv(
    yi = yi_z, 
    V = vi,
    random = list(
      ~ 1 | study_id,
      ~ 1 | es_id, 
    ~ x_km + y_km |const
    ), 
    struct = "SPEXP", 
    data = dat,
    sparse = TRUE,
    verbose = TRUE
    )
)
summary(fit)

# Multivariate Meta-Analysis Model (k = 5000; method: REML)
# 
# logLik    Deviance         AIC         BIC        AICc   
# -1659.0752   3318.1505   3328.1505   3360.7354   3328.1625   
# 
# Variance Components:
#   
#   estim    sqrt  nlvls  fixed    factor 
# sigma^2.1  0.0930  0.3049    500     no  study_id 
# sigma^2.2  0.0507  0.2251   5000     no     es_id 
# 
# outer factor: const        (nlvls = 1)
# inner term:   ~x_km + y_km (nlvls = 500)
# 
# estim    sqrt  fixed 
# tau^2         0.1623  0.4028     no 
# rho        2792.9064             no 
# 
# Test for Heterogeneity:
#   Q(df = 4999) = 82763.9277, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub    
# 0.2743  0.2338  1.1730  0.2408  -0.1840  0.7325    
# 
# ---
dat$es_id <- factor(dat$es_id)
VCV <- diag(dat$vi, nrow = nrow(dat)) 
rownames(VCV)<- colnames(VCV)<- dat$es_id
VCV[1:5, 1:5]

dat$pos <- numFactor(dat$x_km, dat$y_km)

tmb <- system.time(
  glmmTMB(
  yi ~ 1 +
    equalto(0 + effect_id | const, VCV) + 
    (1 | Study) +
    exp(pos + 0 | const),     
  data = dat,
  REML = TRUE
)
)

### ---- ####
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


