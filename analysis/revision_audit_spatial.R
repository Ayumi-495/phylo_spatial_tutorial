#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(geosphere)
  library(metafor)
})

data_path <- file.path("data", "Roger_etal_2024", "Roger_etal_2024.csv")
dat <- read.csv(data_path)

keep <- with(dat,
  is.finite(d_Hedges) &
  is.finite(var_Hedges) &
  var_Hedges > 0 &
  is.finite(latitude) &
  is.finite(longitude)
)
dat <- dat[keep, , drop = FALSE]
dat$study_id <- factor(dat$study_id)
dat$effect_id <- factor(seq_len(nrow(dat)))

coord_key <- sprintf("%.8f_%.8f", dat$latitude, dat$longitude)
location_levels <- unique(coord_key)
dat$location_id <- factor(coord_key, levels = location_levels)
dat$const <- factor(1)

loc <- dat[match(location_levels, coord_key), c("longitude", "latitude")]
D_km <- geosphere::distm(
  as.matrix(loc[, c("longitude", "latitude")]),
  fun = geosphere::distHaversine
) / 1000
rownames(D_km) <- colnames(D_km) <- levels(dat$location_id)

upper_values <- function(x) x[upper.tri(x)]

# Quantify the distortion caused by using global Web Mercator coordinates.
earth_radius_m <- 6378137
lon_rad <- loc$longitude * pi / 180
lat_rad <- pmax(pmin(loc$latitude, 85.05112878), -85.05112878) * pi / 180
x_merc_km <- earth_radius_m * lon_rad / 1000
y_merc_km <- earth_radius_m * log(tan(pi / 4 + lat_rad / 2)) / 1000
D_merc_km <- as.matrix(dist(cbind(x_merc_km, y_merc_km)))
distance_ratio <- upper_values(D_merc_km) / upper_values(D_km)
distance_ratio <- distance_ratio[is.finite(distance_ratio) & upper_values(D_km) > 0]

cat("VERSIONS\n")
cat("R: ", R.version.string, "\n", sep = "")
cat("metafor: ", as.character(packageVersion("metafor")), "\n", sep = "")
cat("geosphere: ", as.character(packageVersion("geosphere")), "\n\n", sep = "")

cat("DATA AUDIT\n")
cat("effect sizes: ", nrow(dat), "\n", sep = "")
cat("studies: ", nlevels(dat$study_id), "\n", sep = "")
cat("unique locations: ", nlevels(dat$location_id), "\n", sep = "")
cat("Hedges d range: ", paste(range(dat$d_Hedges), collapse = " to "), "\n", sep = "")
cat("|d| > 5: ", sum(abs(dat$d_Hedges) > 5), "\n", sep = "")
cat("|d| > 10: ", sum(abs(dat$d_Hedges) > 10), "\n", sep = "")
cat("sampling variances <= 0: ", sum(dat$var_Hedges <= 0), "\n\n", sep = "")

cat("DISTANCE AUDIT\n")
cat("geodesic distance range (km): ", paste(range(upper_values(D_km)), collapse = " to "), "\n", sep = "")
cat(
  "Web Mercator / geodesic ratio (median, 95th percentile, max): ",
  paste(signif(c(
    median(distance_ratio),
    unname(quantile(distance_ratio, 0.95)),
    max(distance_ratio)
  ), 5), collapse = ", "),
  "\n\n",
  sep = ""
)

fit_models <- function(d, label) {
  cat("MODEL SET: ", label, "\n", sep = "")

  common <- list(
    yi = d$d_Hedges,
    V = d$var_Hedges,
    data = d,
    method = "REML",
    test = "t",
    sparse = TRUE,
    control = list(rel.tol = 1e-8)
  )

  unstructured <- do.call(
    rma.mv,
    c(common, list(
      random = list(~ 1 | study_id, ~ 1 | effect_id)
    ))
  )

  spatial_only <- do.call(
    rma.mv,
    c(common, list(
      random = list(~ 1 | study_id, ~ location_id | const),
      struct = "SPEXP",
      dist = list(location_id = D_km)
    ))
  )

  spatial_plus_unstructured <- do.call(
    rma.mv,
    c(common, list(
      random = list(
        ~ 1 | study_id,
        ~ 1 | effect_id,
        ~ location_id | const
      ),
      struct = "SPEXP",
      dist = list(location_id = D_km)
    ))
  )

  fits <- list(
    unstructured = unstructured,
    spatial_only = spatial_only,
    spatial_plus_unstructured = spatial_plus_unstructured
  )

  tab <- do.call(rbind, lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    data.frame(
      model = nm,
      k = fit$k,
      estimate = as.numeric(coef(fit)[1]),
      se = fit$se[1],
      ci_lb = fit$ci.lb[1],
      ci_ub = fit$ci.ub[1],
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit),
      sigma2 = paste(signif(fit$sigma2, 6), collapse = ";"),
      tau2 = if (is.null(fit$tau2)) NA_character_ else paste(signif(fit$tau2, 6), collapse = ";"),
      rho = if (is.null(fit$rho)) NA_character_ else paste(signif(fit$rho, 6), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  print(tab, row.names = FALSE)
  cat("\n")

  invisible(fits)
}

fits_full <- fit_models(dat, "all eligible effect sizes")

# A transparent diagnostic sensitivity analysis. This threshold is not a
# biological definition of implausibility and must not replace source-level
# verification of the effect-size calculations.
dat_sensitivity <- droplevels(dat[abs(dat$d_Hedges) <= 5, , drop = FALSE])
fits_sensitivity <- fit_models(dat_sensitivity, "diagnostic sensitivity: |Hedges d| <= 5")

cat("INTERPRETATION GUARDRAILS\n")
cat("metafor SPEXP uses cor(d) = exp(-d / rho); rho is a distance scale.\n")
cat("A spatial random effect and an unstructured effect-level random effect answer different questions.\n")
cat("Study/effect random effects do not encode non-diagonal sampling covariance.\n")
cat("The |d| <= 5 analysis is diagnostic only; extreme values require source-level verification.\n")
