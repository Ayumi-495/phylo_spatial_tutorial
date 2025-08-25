estimate_coords_from_dist <- function(
    dist_matrix, # A square matrix of pairwise distances (e.g., phylogenetic distances). Must be symmetric and non-negative.
    ref_points_df,#A data frame with reference points. Must have 3 columns: "name" (matching matrix names), "lon" (longitude), and "lat" (latitude).
    point_names = NULL, # Optional vector of point names. If NULL, names are taken from the rownames of dist_matrix.
    apply_scaling = TRUE, # If TRUE, the coordinates are scaled to match the distance between the first two reference points.
    max_scale_factor = 100, # Maximum allowed scale factor for coordinate stretching. Set NULL for no limit. Helps prevent extreme values.
    scale_mds = TRUE # If TRUE, standardize the MDS result (mean = 0, SD = 1 for each axis) before transformation. Recommended for stability.

) {
  # --- Input Validation ---
  if (!is.matrix(dist_matrix) || !is.numeric(dist_matrix)) {
    stop("`dist_matrix` must be a numeric matrix.")
  }
  if (nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("`dist_matrix` must be a square matrix.")
  }
  if (any(dist_matrix < 0)) {
    stop("Distances in `dist_matrix` cannot be negative.")
  }
  if (!all(diag(dist_matrix) == 0)) {
    warning("Diagonal elements of `dist_matrix` are not all zero. MDS assumes self-distance is zero.")
  }
  if (!is.data.frame(ref_points_df) || !all(c("name", "lon", "lat") %in% names(ref_points_df))) {
    stop("`ref_points_df` must be a data frame with columns 'name', 'lon', 'lat'.")
  }
  
  if (is.null(point_names)) {
    if (is.null(rownames(dist_matrix))) {
      point_names <- paste0("P", 1:nrow(dist_matrix))
      warning("`point_names` not provided and `dist_matrix` has no row names. Defaulting to P1, P2, ...")
    } else {
      point_names <- rownames(dist_matrix)
    }
  } else if (length(point_names) != nrow(dist_matrix)) {
    stop("Length of `point_names` must match the number of points in `dist_matrix`.")
  }
  
  rownames(dist_matrix) <- colnames(dist_matrix) <- point_names
  
  # --- 1. MDS ---
  mds_result <- cmdscale(d = dist_matrix, k = 2)
  if (scale_mds) {
    mds_result <- scale(mds_result)  # standardize
  }
  colnames(mds_result) <- c("MDS_Dim1", "MDS_Dim2")
  mds_df <- as.data.frame(mds_result)
  mds_df$name <- point_names
  
  # --- 2. Align to Reference Points ---
  mds_ref_points <- mds_df[match(ref_points_df$name, mds_df$name), ]
  if (any(is.na(mds_ref_points))) {
    stop("One or more reference point names in `ref_points_df` not found in `dist_matrix` names.")
  }
  
  mds_final_coords <- mds_df
  mds_final_coords$lon_estimated <- NA
  mds_final_coords$lat_estimated <- NA
  
  if (nrow(ref_points_df) >= 1) {
    # Translation
    ref1 <- ref_points_df[1, ]
    ref1_mds <- mds_ref_points[1, c("MDS_Dim1", "MDS_Dim2")]
    translation <- c(ref1$lon - ref1_mds$MDS_Dim1, ref1$lat - ref1_mds$MDS_Dim2)
    
    mds_translated <- mds_df
    mds_translated$MDS_Dim1 <- mds_df$MDS_Dim1 + translation[1]
    mds_translated$MDS_Dim2 <- mds_df$MDS_Dim2 + translation[2]
    
    if (nrow(ref_points_df) >= 2 && apply_scaling) {
      ref2 <- ref_points_df[2, ]
      ref2_mds <- mds_translated[mds_translated$name == ref2$name, ]
      
      vec_mds <- c(ref2_mds$MDS_Dim1 - ref1$lon, ref2_mds$MDS_Dim2 - ref1$lat)
      vec_true <- c(ref2$lon - ref1$lon, ref2$lat - ref1$lat)
      
      len_mds <- sqrt(sum(vec_mds^2))
      len_true <- sqrt(sum(vec_true^2))
      
      scale_factor <- if (len_mds == 0) 1 else len_true / len_mds
      if (!is.null(max_scale_factor)) {
        scale_factor <- min(scale_factor, max_scale_factor)
      }
      
      angle_mds <- atan2(vec_mds[2], vec_mds[1])
      angle_true <- atan2(vec_true[2], vec_true[1])
      rotation_angle <- angle_true - angle_mds
      
      rot_matrix <- matrix(c(cos(rotation_angle), -sin(rotation_angle),
                             sin(rotation_angle),  cos(rotation_angle)), nrow = 2, byrow = TRUE)
      
      centered <- as.matrix(mds_translated[, c("MDS_Dim1", "MDS_Dim2")]) - matrix(rep(c(ref1$lon, ref1$lat), each = nrow(mds_translated)), ncol = 2, byrow = FALSE)
      rotated_scaled <- (centered %*% rot_matrix) * scale_factor
      
      mds_final_coords$lon_estimated <- rotated_scaled[, 1] + ref1$lon
      mds_final_coords$lat_estimated <- rotated_scaled[, 2] + ref1$lat
    } else {
      # translation only
      mds_final_coords$lon_estimated <- mds_translated$MDS_Dim1
      mds_final_coords$lat_estimated <- mds_translated$MDS_Dim2
    }
  }
  
  return(mds_final_coords[, c("name", "lon_estimated", "lat_estimated")])
}

ref_points <- data.frame(
  name = c("Homo_sapiens_ott770315", "Megaptera_novaeangliae_ott226198"),
  lon = c(0, 1),  
  lat = c(0, 0) 
)

test <- estimate_coords_from_dist(
  dist_matrix = D,
  ref_points_df = ref_points,
  apply_scaling = TRUE,
  scale_mds = FALSE
)

summary(test)
head(test, 20)
plot(test$lon_estimated, test$lat_estimated, asp = 1, pch = 20,
     xlab = "Estimated longitude", ylab = "Estimated latitude")

# merge the dataset with coordinate distance ----
setdiff(unique(dat$species.id), test$name)
test$species.id <- test$name
dat_with_coords <- merge(dat, test[, c("species.id", "lon_estimated", "lat_estimated")],
                         by = "species.id", all.x = TRUE)

# gaussian process model in brms ----
formula_a <-  brms::bf(yi | se(sqrt(vi)) ~ 1 + 
                         gp(lon_estimated, lat_estimated, cov = "exp_quad") +
                         (1 | study.id) +
                         (1 | effect.size.id) + 
                         (1 | species.id)
)

prior_ma_a <- get_prior(
  formula = formula_1,
  data = dat_with_coords,
  family = gaussian()
)

system.time(
  fit_brm1 <- brm(
    formula_1,
    data = dat_with_coords,
    chains = 2, 
    iter = 10000, 
    warmup = 8000,
    prior = prior_ma_a,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.95)
  ))
# user    system   elapsed 
# 
summary(fit_brm1)


test2 <- estimate_coords_from_dist(
  dist_matrix = D,
  ref_points_df = ref_points,
  apply_scaling = TRUE,
  scale_mds = TRUE
)

test3 <- estimate_coords_from_dist(
  dist_matrix = D,
  ref_points_df = ref_points,
  apply_scaling = FALSE,
  scale_mds = TRUE
)

summary(test3)

plot(test3$lon_estimated, test3$lat_estimated, asp = 1, pch = 20,
     xlab = "Estimated longitude", ylab = "Estimated latitude")

# merge the dataset with coordinate distance ----
setdiff(unique(dat$species.id), test3$name)
test3$species.id <- test3$name
dat_with_coords <- merge(dat, test3[, c("species.id", "lon_estimated", "lat_estimated")],
                         by = "species.id", all.x = TRUE)

# gaussian process model in brms ----
formula_b <-  brms::bf(yi | se(sqrt(vi)) ~ 1 + 
                         gp(lon_estimated, lat_estimated, cov = "exp_quad") +
                         (1 | study.id) +
                         (1 | effect.size.id) + 
                         (1 | species.id)
)

prior_ma_b <- get_prior(
  formula = formula_b,
  data = dat_with_coords,
  family = gaussian()
)

system.time(
  fit_brm2 <- brm(
    formula_b,
    data = dat_with_coords,
    chains = 2, 
    iter = 2000, 
    warmup = 1000,
    prior = prior_ma_b,
    # backend = "cmdstanr",
    control = list(adapt_delta = 0.95)
  ))
# user    system   elapsed 
# 
summary(fit_brm2)
