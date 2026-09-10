# Reviewer Comment 18 audit: recover the source-publication influential-effect
# exclusions and fit corrected global spatial models to the published-cleaned
# spatially usable Grau-Andres data. This script never edits tutorial_v2.qmd.

suppressPackageStartupMessages({
  library(metafor)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript reviewer18_influential_effects.R <data_csv> <authors_code_R> <output_dir>")
}
data_csv <- args[[1L]]
authors_code <- args[[2L]]
out_dir <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(data_csv), file.exists(authors_code), dir.exists(out_dir))

dryad_expected <- c(
  data = "435c87cac364acc3bfb61e978b83c6bed06a5b08a2336c15e5f7427ed8674b45",
  code = "9a7aa3cbd14faeb59a8dd6773f52524db2a240637a40616a5cd8d42db80d015c"
)

sha256_file <- function(path) {
  exe <- Sys.which("sha256sum")
  if (nzchar(exe)) {
    return(strsplit(system2(exe, path, stdout = TRUE), "[[:space:]]+")[[1L]][1L])
  }
  exe <- Sys.which("shasum")
  if (!nzchar(exe)) stop("Neither sha256sum nor shasum is available")
  strsplit(system2(exe, c("-a", "256", path), stdout = TRUE), "[[:space:]]+")[[1L]][1L]
}

actual_hash <- c(data = sha256_file(data_csv), code = sha256_file(authors_code))
stopifnot(identical(unname(actual_hash), unname(dryad_expected)))
write.csv(data.frame(file = names(actual_hash), expected_sha256 = dryad_expected,
                     actual_sha256 = actual_hash,
                     exact_match = actual_hash == dryad_expected),
          file.path(out_dir, "dryad_file_hashes.csv"), row.names = FALSE)

# These are transcribed from the executable exclusion statements and adjacent
# Cook's-distance comments in the archived authors' code. Thresholds differed
# by response category; no new cutoff is introduced here.
excluded_keys <- data.frame(
  effect_uid = c("Ngugi_2022-2", "Gagnon_2015-2", "Moris_2017-1",
                 "Schwilk_1997-1", "Silveira_2016-4",
                 "Launonen_1999-1", "Ansley_2015-1"),
  study_id = c("Ngugi_2022", "Gagnon_2015", "Moris_2017",
               "Schwilk_1997", "Silveira_2016",
               "Launonen_1999", "Ansley_2015"),
  ES_num = c(2L, 2L, 1L, 1L, 4L, 1L, 1L),
  response_expected = c(rep("abundance", 3L), rep("diversity", 2L), rep("fitness", 2L)),
  cooks_threshold = c(rep(0.015, 5L), rep(0.06, 2L)),
  stringsAsFactors = FALSE
)

dat_all <- read.csv(data_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(dat_all) == 2363L,
          all(c("study_id", "ES_num", "latitude", "longitude", "response",
                "d_Hedges", "var_Hedges", "imputed", "source") %in% names(dat_all)))
dat_all$effect_uid <- paste(dat_all$study_id, dat_all$ES_num, sep = "-")
stopifnot(!anyDuplicated(dat_all$effect_uid))

hit <- match(excluded_keys$effect_uid, dat_all$effect_uid)
excluded <- cbind(excluded_keys, data.frame(
  present_in_dryad_csv = !is.na(hit),
  study_num = dat_all$study_num[hit],
  response = dat_all$response[hit],
  d_Hedges = dat_all$d_Hedges[hit],
  var_Hedges = dat_all$var_Hedges[hit],
  imputed = dat_all$imputed[hit],
  source = dat_all$source[hit],
  country = dat_all$country[hit],
  latitude = dat_all$latitude[hit],
  longitude = dat_all$longitude[hit],
  stringsAsFactors = FALSE
))
stopifnot(sum(excluded$present_in_dryad_csv) == 6L,
          identical(excluded$effect_uid[!excluded$present_in_dryad_csv], "Launonen_1999-1"),
          all(excluded$response[excluded$present_in_dryad_csv] ==
                excluded$response_expected[excluded$present_in_dryad_csv]))

spatial_ok_all <- complete.cases(dat_all[, c("latitude", "longitude")])
dat_spatial <- dat_all[spatial_ok_all, , drop = FALSE]
stopifnot(nrow(dat_spatial) == 2361L)
missing_coord <- dat_all[!spatial_ok_all,
                         c("effect_uid", "study_id", "study_num", "ES_num", "response",
                           "d_Hedges", "var_Hedges", "imputed", "source"), drop = FALSE]

excluded$spatially_usable <- excluded$effect_uid %in% dat_spatial$effect_uid
excluded$among_two_missing_coordinate_records <- excluded$effect_uid %in% missing_coord$effect_uid
abs_order_all <- dat_all$effect_uid[order(abs(dat_all$d_Hedges), decreasing = TRUE, na.last = NA)]
abs_order_spatial <- dat_spatial$effect_uid[order(abs(dat_spatial$d_Hedges), decreasing = TRUE, na.last = NA)]
excluded$absolute_effect_rank_all_dryad <- match(excluded$effect_uid, abs_order_all)
excluded$absolute_effect_rank_spatial <- match(excluded$effect_uid, abs_order_spatial)
excluded$within_top_7_absolute_spatial <- excluded$absolute_effect_rank_spatial <= 7L
write.csv(excluded, file.path(out_dir, "published_influential_effects.csv"), row.names = FALSE)
write.csv(missing_coord, file.path(out_dir, "missing_coordinate_records.csv"), row.names = FALSE)

top_abs <- dat_spatial[order(abs(dat_spatial$d_Hedges), decreasing = TRUE),
                       c("effect_uid", "study_id", "study_num", "ES_num", "response",
                         "d_Hedges", "var_Hedges", "imputed", "source"), drop = FALSE]
top_abs$absolute_effect_rank_spatial <- seq_len(nrow(top_abs))
write.csv(head(top_abs, 30L), file.path(out_dir, "top_30_absolute_spatial_effects.csv"), row.names = FALSE)

present_exclusions <- excluded$effect_uid[excluded$present_in_dryad_csv]
published_cleaned <- dat_spatial[!dat_spatial$effect_uid %in% present_exclusions, , drop = FALSE]
all_spatially_usable <- dat_spatial
stopifnot(nrow(published_cleaned) == nrow(all_spatially_usable) - 6L,
          !any(published_cleaned$effect_uid %in% present_exclusions))
write.csv(all_spatially_usable, file.path(out_dir, "all_spatially_usable.csv"), row.names = FALSE)
write.csv(published_cleaned, file.path(out_dir, "published_cleaned.csv"), row.names = FALSE)

count_dataset <- function(x, name) {
  data.frame(dataset = name, n_effects = nrow(x),
             n_studies = length(unique(x$study_id)),
             n_sites = length(unique(sprintf("%.8f_%.8f", x$latitude, x$longitude))),
             n_spain = sum(!is.na(x$country) & x$country == "Spain"),
             stringsAsFactors = FALSE)
}
dataset_counts <- rbind(count_dataset(all_spatially_usable, "all_spatially_usable"),
                        count_dataset(published_cleaned, "published_cleaned"))
write.csv(dataset_counts, file.path(out_dir, "dataset_counts.csv"), row.names = FALSE)

spain_effects <- all_spatially_usable$effect_uid[
  !is.na(all_spatially_usable$country) & all_spatially_usable$country == "Spain"]
spain_overlap <- excluded[excluded$effect_uid %in% spain_effects, , drop = FALSE]
write.csv(spain_overlap, file.path(out_dir, "spain_exclusion_overlap.csv"), row.names = FALSE)

prepare_model_data <- function(x) {
  x$effect_id <- factor(seq_len(nrow(x)), levels = seq_len(nrow(x)))
  x$study_id <- factor(x$study_id)
  x$site_key <- sprintf("%.8f_%.8f", x$latitude, x$longitude)
  site_levels <- sort(unique(x$site_key))
  x$site_id <- factor(x$site_key, levels = site_levels)
  x$const <- factor("all_sites")
  site_lookup <- unique(x[c("site_key", "latitude", "longitude")])
  site_lookup <- site_lookup[match(site_levels, site_lookup$site_key), , drop = FALSE]
  site_lookup$site_id <- site_levels
  coords <- as.matrix(site_lookup[c("longitude", "latitude")])
  distance_km <- geosphere::distm(coords, fun = geosphere::distGeo) / 1000
  rownames(distance_km) <- colnames(distance_km) <- site_levels
  stopifnot(identical(levels(x$site_id), site_levels),
            identical(rownames(distance_km), levels(x$site_id)),
            identical(colnames(distance_km), levels(x$site_id)),
            isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
            all(abs(diag(distance_km)) < 1e-10), all(distance_km >= 0),
            all(is.finite(distance_km)), all(x$var_Hedges > 0))
  list(dat = x, site_lookup = site_lookup, distance_km = distance_km)
}

p <- prepare_model_data(published_cleaned)
saveRDS(p, file.path(out_dir, "published_cleaned_prepared.rds"))
write.csv(p$site_lookup, file.path(out_dir, "published_cleaned_site_lookup.csv"), row.names = FALSE)
write.csv(p$distance_km, file.path(out_dir, "published_cleaned_distance_km.csv"), row.names = TRUE)

capture_result <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(force(expr), error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = unique(warnings))
}

fit_status <- function(fit, warnings) {
  if (!is.null(fit$converged)) return(if (isTRUE(fit$converged)) "converged" else "not_converged")
  if (!is.null(fit$optres$convergence)) {
    return(if (identical(fit$optres$convergence, 0L)) "converged" else "not_converged")
  }
  if (any(grepl("converg", warnings, ignore.case = TRUE))) return("optimizer_warning")
  "completed_no_explicit_optimizer_status"
}

save_fit <- function(name, captured, elapsed) {
  if (inherits(captured$value, "error")) {
    row <- data.frame(model = name, n_effects = nrow(p$dat),
                      n_studies = nlevels(p$dat$study_id), n_sites = nlevels(p$dat$site_id),
                      mean = NA, ci_lb = NA, ci_ub = NA, effect_variance = NA,
                      study_variance = NA, spatial_variance = NA, rho_km = NA,
                      logLik_REML = NA, AIC_REML = NA, convergence_status = "error",
                      elapsed_seconds = elapsed,
                      warnings = paste(c(conditionMessage(captured$value), captured$warnings), collapse = " | "))
  } else {
    fit <- captured$value
    saveRDS(fit, file.path(out_dir, paste0(name, ".rds")))
    row <- data.frame(model = name, n_effects = fit$k,
                      n_studies = nlevels(p$dat$study_id), n_sites = nlevels(p$dat$site_id),
                      mean = as.numeric(fit$b[1]), ci_lb = as.numeric(fit$ci.lb[1]),
                      ci_ub = as.numeric(fit$ci.ub[1]), effect_variance = fit$sigma2[1],
                      study_variance = if (name %in% c("unstructured_only", "combined")) fit$sigma2[2] else NA,
                      spatial_variance = if (name %in% c("spatial_only", "combined")) fit$tau2[1] else NA,
                      rho_km = if (name %in% c("spatial_only", "combined")) fit$rho[1] else NA,
                      logLik_REML = as.numeric(fit$fit.stats["ll", "REML"]),
                      AIC_REML = as.numeric(fit$fit.stats["AIC", "REML"]),
                      convergence_status = fit_status(fit, captured$warnings),
                      elapsed_seconds = elapsed,
                      warnings = paste(captured$warnings, collapse = " | "))
  }
  path <- file.path(out_dir, "published_cleaned_primary_results.csv")
  old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$model != name, , drop = FALSE]
  write.csv(rbind(old, row), path, row.names = FALSE)
  invisible(captured$value)
}

fit_one <- function(name, expr) {
  start <- proc.time()[["elapsed"]]
  captured <- capture_result(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  save_fit(name, captured, elapsed)
  message(name, " finished in ", round(elapsed, 1), " seconds")
}

# Assert matrix ordering immediately before every fit, as required.
assert_spatial_order <- function() {
  stopifnot(identical(rownames(p$distance_km), levels(p$dat$site_id)),
            identical(colnames(p$distance_km), levels(p$dat$site_id)))
}

assert_spatial_order()
fit_one("unstructured_only", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~1 | effect_id, ~1 | study_id),
  data = p$dat, method = "REML", test = "t", sparse = TRUE))

assert_spatial_order()
fit_one("spatial_only", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~1 | effect_id, ~site_id | const),
  struct = "SPEXP", dist = list(site_id = p$distance_km),
  data = p$dat, method = "REML", test = "t", sparse = TRUE))

assert_spatial_order()
fit_one("combined", rma.mv(
  yi = d_Hedges, V = var_Hedges,
  random = list(~1 | effect_id, ~1 | study_id, ~site_id | const),
  struct = "SPEXP", dist = list(site_id = p$distance_km),
  data = p$dat, method = "REML", test = "t", sparse = TRUE))

writeLines(c(
  paste("R:", R.version.string),
  paste("metafor:", as.character(packageVersion("metafor"))),
  paste("geosphere:", as.character(packageVersion("geosphere"))),
  paste("Dryad data SHA-256:", actual_hash[["data"]]),
  paste("Dryad code SHA-256:", actual_hash[["code"]]),
  "Distance: geosphere::distm(..., fun=geosphere::distGeo)/1000; WGS84 ellipsoidal geodesic kilometres.",
  "Published-cleaned means the 2361 spatially usable Dryad rows minus the six code-listed records present in the Dryad CSV.",
  "Launonen_1999-1 is named in the archived code but is already absent from the archived Dryad CSV; no replacement record was invented.",
  "Sampling V is diagonal var_Hedges; fixed effects are intercept-only; all models include iid effect-size heterogeneity.",
  "No additional iid site intercept was fitted; spatial outer group const has one level."
), file.path(out_dir, "audit_metadata.txt"))
