pacman:: p_load(tidyverse, brms, metafor, geiger, mvtnorm)

dat <- dat.moura2021$dat
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)
tree <- dat.moura2021$tree
tree <- compute.brlen(tree, method="Grafen", power=1)

summary(dat)
tree <- dat.moura2021$tree
tree <- compute.brlen(tree, method="Grafen", power=1)

A_BM <- vcv(tree, corr=TRUE)
A_BM[1:5, 1:5]
I <- matrix(1, nrow = 341, ncol = 341)

D <- I - A_BM # make distance matrix
D[1:20, 1:20]


# estimate coordinates from phylogenetic distance matrix using MDS ----
estimate_coords_from_dist <- function(dist_matrix, ref_points_df, point_names = NULL) {
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
  if (nrow(ref_points_df) < 1) {
    stop("At least one reference point must be provided in `ref_points_df`.")
  }
  # Assign point names from dist_matrix if not provided
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
  # Ensure dist_matrix has names for matching
  rownames(dist_matrix) <- colnames(dist_matrix) <- point_names
  # --- 1. Perform Classical Multidimensional Scaling (MDS) ---
  mds_result <- cmdscale(d = dist_matrix, k = 2) # k=2 for 2D coordinates (lon/lat)
  colnames(mds_result) <- c("MDS_Dim1", "MDS_Dim2")
  mds_df <- as.data.frame(mds_result)
  mds_df$name <- point_names
  # --- 2. Align MDS Output to Reference Points ---
  num_ref_points <- nrow(ref_points_df)
  # Extract MDS coordinates for reference points
  mds_ref_points <- mds_df[match(ref_points_df$name, mds_df$name), ]
  if (any(is.na(mds_ref_points))) {
    stop("One or more reference point names in `ref_points_df` not found in `dist_matrix` names.")
  }
  # Get true coordinates for reference points
  true_ref_points <- ref_points_df
  # Initialize translated and final coordinates
  mds_final_coords <- mds_df
  mds_final_coords$lon_estimated <- NA
  mds_final_coords$lat_estimated <- NA
  if (num_ref_points >= 1) {
    # Always translate using the first reference point
    ref1_mds_x <- mds_ref_points$MDS_Dim1[1]
    ref1_mds_y <- mds_ref_points$MDS_Dim2[1]
    ref1_true_lon <- true_ref_points$lon[1]
    ref1_true_lat <- true_ref_points$lat[1]
    # Translation vector
    translation_x <- ref1_true_lon - ref1_mds_x
    translation_y <- ref1_true_lat - ref1_mds_y
    mds_translated_x <- mds_df$MDS_Dim1 + translation_x
    mds_translated_y <- mds_df$MDS_Dim2 + translation_y
    if (num_ref_points >= 2) {
      # Use second reference point for rotation and scaling
      ref2_mds_x_translated <- mds_translated_x[mds_df$name == true_ref_points$name[2]]
      ref2_mds_y_translated <- mds_translated_y[mds_df$name == true_ref_points$name[2]]
      ref2_true_lon <- true_ref_points$lon[2]
      ref2_true_lat <- true_ref_points$lat[2]
      # Vector from ref1 to ref2 in translated MDS space
      vec_mds_x <- ref2_mds_x_translated - ref1_true_lon
      vec_mds_y <- ref2_mds_y_translated - ref1_true_lat
      # Vector from ref1 to ref2 in true geographical space
      vec_true_x <- ref2_true_lon - ref1_true_lon
      vec_true_y <- ref2_true_lat - ref1_true_lat
      # Calculate lengths of the vectors
      len_mds <- sqrt(vec_mds_x^2 + vec_mds_y^2)
      len_true <- sqrt(vec_true_x^2 + vec_true_y^2)
      # Calculate scaling factor
      scale_factor <- if (len_mds == 0) 1 else len_true / len_mds
      # Calculate angles
      angle_mds <- atan2(vec_mds_y, vec_mds_x)
      angle_true <- atan2(vec_true_y, vec_true_x)
      # Calculate rotation angle needed
      rotation_angle <- angle_true - angle_mds
      # Create rotation matrix
      rotation_matrix <- matrix(c(cos(rotation_angle), -sin(rotation_angle),
                                  sin(rotation_angle),  cos(rotation_angle)), nrow = 2, byrow = TRUE)
      # Shift points so ref1 is at the origin for rotation and scaling
      points_centered_x <- mds_translated_x - ref1_true_lon
      points_centered_y <- mds_translated_y - ref1_true_lat
      points_centered <- as.matrix(cbind(points_centered_x, points_centered_y))
      # Apply scaling and rotation
      rotated_scaled_points <- (points_centered %*% rotation_matrix) * scale_factor
      # Shift back by adding ref1's true coordinates
      mds_final_coords$lon_estimated <- rotated_scaled_points[,1] + ref1_true_lon
      mds_final_coords$lat_estimated <- rotated_scaled_points[,2] + ref1_true_lat
    } else {
      # Only one reference point, only translation applied
      mds_final_coords$lon_estimated <- mds_translated_x
      mds_final_coords$lat_estimated <- mds_translated_y
    }
  }
  return(mds_final_coords[, c("name", "lon_estimated", "lat_estimated")])
}

# try the function with Moura et al. 2021 dataset ----
ref_points <- data.frame(
  name = c("Homo_sapiens_ott770315", "Megaptera_novaeangliae_ott226198"),
  lon = c(0, 1),  
  lat = c(0, 0) 
)

estimated_coords <- estimate_coords_from_dist(D, ref_points)
names(estimated_coords)
# [1] "name"          "lon_estimated" "lat_estimated" "species.id" 
head(estimated_coords, 20)
D[1:5, 1:5]
# merge the dataset with coordinate distance ----
setdiff(unique(dat$species.id), estimated_coords$name)
estimated_coords$species.id <- estimated_coords$name
dat_with_coords <- merge(dat, estimated_coords[, c("species.id", "lon_estimated", "lat_estimated")],
                         by = "species.id", all.x = TRUE)

# gaussian process model in brms ----
formula_1 <-  brms::bf(yi | se(sqrt(vi)) ~ 1 + 
                     gp(lon_estimated, lat_estimated, cov = "exp_quad") +
                     (1 | study.id) +
                     (1 | effect.size.id) + 
                     (1 | species.id)
                     )

prior_ma1 <- get_prior(
  formula = formula_1,
  data = dat_with_coords,
  family = gaussian()
)

system.time(
  fit_brm1 <- brm(
               formula_1,
               data = dat_with_coords,
               chains = 2, 
               iter = 3000, 
               warmup = 2000,
               prior = prior_ma1,
               # backend = "cmdstanr",
               control = list(adapt_delta = 0.95)
                  ))
# user    system   elapsed 
# 
summary(fit_brm1)
