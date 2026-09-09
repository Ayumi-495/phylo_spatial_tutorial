# Correctness audit for the full Scholer et al. spatial example.
# This file deliberately performs no model fitting and does not edit tutorial_v2.qmd.

suppressPackageStartupMessages({
  library(dplyr)
  library(geosphere)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scholer_structure_audit.R <data_csv> <output_dir>")
}
data_csv <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(data_csv, stringsAsFactors = FALSE)
required <- c("X", "ref", "logit_survival", "se", "lat", "long")
stopifnot(all(required %in% names(dat)))

# The dataset has one row per effect size. Each unique recorded coordinate pair
# forms a location ID. This does not claim equal coordinates are the same exact
# field site; it records the coordinate resolution supplied by the dataset.
dat <- dat |>
  mutate(
    vi = se^2,
    effect_id = factor(sprintf("effect_%04d", seq_len(n()))),
    study_id = factor(ref),
    site_key = sprintf("%.8f_%.8f", lat, long)
  )
site_levels <- sort(unique(dat$site_key))
dat$site_id <- factor(dat$site_key, levels = site_levels)
dat$const <- factor("all_sites")

stopifnot(
  nrow(dat) == 949L,
  nlevels(dat$study_id) == 205L,
  nlevels(dat$site_id) == 454L,
  all(is.finite(dat$logit_survival)),
  all(is.finite(dat$se)), all(dat$se > 0),
  all(is.finite(dat$vi)), all(dat$vi > 0),
  all(is.finite(dat$lat)), all(is.finite(dat$long)),
  all(abs(dat$lat) <= 90), all(abs(dat$long) <= 180)
)

# Audit the historical group_by(lat, long) site construction: the numeric labels
# depend on grouping order, but its partition of rows should match site_key.
legacy <- dat |>
  group_by(lat, long) |>
  mutate(legacy_site_id = cur_group_id()) |>
  ungroup()
legacy_mapping <- legacy |>
  distinct(site_id, legacy_site_id)
legacy_one_to_one <- nrow(legacy_mapping) == nlevels(dat$site_id) &&
  n_distinct(legacy_mapping$site_id) == nrow(legacy_mapping) &&
  n_distinct(legacy_mapping$legacy_site_id) == nrow(legacy_mapping)
stopifnot(legacy_one_to_one)

site_lookup <- dat |>
  distinct(site_id, site_key, lat, long) |>
  arrange(site_id)
stopifnot(identical(as.character(site_lookup$site_id), levels(dat$site_id)))

# WGS84 ellipsoidal great-circle distance in kilometres, aligned explicitly
# to the factor levels that will be used in the spatial random term.
coords_lonlat <- as.matrix(site_lookup[c("long", "lat")])
distance_km <- geosphere::distm(coords_lonlat, fun = geosphere::distGeo) / 1000
rownames(distance_km) <- colnames(distance_km) <- levels(dat$site_id)
stopifnot(
  identical(rownames(distance_km), levels(dat$site_id)),
  identical(colnames(distance_km), levels(dat$site_id)),
  isTRUE(all.equal(distance_km, t(distance_km), tolerance = 1e-10)),
  all(abs(diag(distance_km)) < 1e-10), all(is.finite(distance_km)),
  all(distance_km >= 0)
)

study_locations <- dat |>
  group_by(study_id) |>
  summarise(n_effects = n(), n_sites = n_distinct(site_id), .groups = "drop") |>
  arrange(desc(n_sites), desc(n_effects))
shared_locations <- dat |>
  group_by(site_id) |>
  summarise(n_effects = n(), n_studies = n_distinct(study_id), .groups = "drop") |>
  arrange(desc(n_studies), desc(n_effects))

summary <- data.frame(
  effect_sizes = nrow(dat),
  references = nlevels(dat$study_id),
  recorded_coordinate_locations = nlevels(dat$site_id),
  studies_with_multiple_locations = sum(study_locations$n_sites > 1),
  locations_shared_across_studies = sum(shared_locations$n_studies > 1),
  rows_at_shared_locations = sum(shared_locations$n_effects[shared_locations$n_studies > 1]),
  maximum_locations_per_study = max(study_locations$n_sites),
  maximum_studies_per_location = max(shared_locations$n_studies),
  maximum_great_circle_distance_km = max(distance_km),
  historical_site_id_partition_correct = legacy_one_to_one,
  stringsAsFactors = FALSE
)

write.csv(summary, file.path(out_dir, "scholer_structure_summary.csv"), row.names = FALSE)
write.csv(study_locations, file.path(out_dir, "scholer_study_location_counts.csv"), row.names = FALSE)
write.csv(shared_locations, file.path(out_dir, "scholer_shared_location_counts.csv"), row.names = FALSE)
write.csv(site_lookup, file.path(out_dir, "scholer_site_lookup.csv"), row.names = FALSE)
write.csv(distance_km, file.path(out_dir, "scholer_distance_great_circle_km.csv"), row.names = TRUE)
saveRDS(list(dat = dat, site_lookup = site_lookup, distance_km = distance_km,
             study_locations = study_locations, shared_locations = shared_locations),
        file.path(out_dir, "scholer_prepared.rds"))

writeLines(c(
  "Existing QMD correctness check (no model run):",
  "1. The historical preparation code creates a 949 x 949 EPSG:3857 effect-level distance matrix, then a separate 454 x 454 site distance matrix.",
  "2. The displayed metafor call specifies ~ effect_id | const (949 inner levels) but supplies dist = list(site = dist_site), whose 454 row/column names refer to site IDs.",
  "3. The displayed model output instead reports inner term ~site_id (454 levels). It cannot be the output of the displayed call as written.",
  "4. confint(ma_sp2_exp) refers to an undefined object; the displayed fitted object is sp_exp_metafor_eg2.",
  "5. The prose says 204 studies, while the current data contain 205 unique ref values. The section header says Scholer et al. (2021), while its citation and linked paper say 2020.",
  "Conclusion: no existing numerical output in this section should be reused. The corrected future spatial term must be ~ site_id | const with dist = list(site_id = distance_km).",
  "No model fitting was performed by this audit."
), file.path(out_dir, "scholer_existing_qmd_consistency.txt"))

message("Scholer structure audit completed; no models were fitted.")
