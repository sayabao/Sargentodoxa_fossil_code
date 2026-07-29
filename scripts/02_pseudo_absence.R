# ============================================================================
# 02_pseudo_absence.R - Pseudo-absence generation (two schemes)
# ============================================================================
# Functionality:
#   - Generate pseudo-absences for two schemes via biomod2 built-in functions
#   - Strategy: random, global land extent
#   - Inner buffer: 100 km (lon/lat: dist.min = pa_dist_min degrees)
#
# Two schemes (read from config.R):
#   1. Scheme 1: 1:1 ratio (nb.absences = equal, nb.rep = 1)
#   2. Scheme 2: ~1:3 ratio (nb.absences = 1000, nb.rep = 1)
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  02 Pseudo Absence Generation - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Read data
# ============================================================================

# Read presence SpatVector
presence_file <- file.path(output_dirs$data_prep, "sargentodoxa_presence_spat.rds")
if (!file.exists(presence_file)) {
  stop("Error: Presence data not found. Run 01_data_preparation.R first.")
}

presence_spat <- readRDS(presence_file)
n_presence <- nrow(presence_spat)

cat("  Loading presence data:", presence_file, "\n")
cat("  Presence points:", n_presence, "\n")

# Read modern climate data
cat("\n  Loading modern climate data from:", modern_clim_dir, "\n")
clim_files <- list.files(modern_clim_dir, pattern = modern_clim_pattern, full.names = TRUE)

if (length(clim_files) == 0) {
  stop("Error: No climate files found in ", modern_clim_dir)
}

cat("  Found", length(clim_files), "climate files\n")

# Read climate rasters and order them
modern_clim <- rast(clim_files)
modern_clim <- modern_clim[[mixedsort(names(modern_clim))]]
names(modern_clim) <- paste0("Bio", 1:nlyr(modern_clim))

# Build a "complete environmental values" mask for pseudo-absences, not just the first layer
# Prefer manual_var_select so variables used later are all non-NA
mask_vars <- if (!is.null(manual_var_select) && length(manual_var_select) > 0) {
  intersect(manual_var_select, names(modern_clim))
} else {
  names(modern_clim)
}

if (length(mask_vars) == 0) {
  stop("Error: No valid variables available to build pseudo-absence mask.")
}

cat("  Building valid-environment mask with", length(mask_vars), "variables:\n")
cat("   ", paste(mask_vars, collapse = ", "), "\n")

mask_stack <- modern_clim[[mask_vars]]
valid_mask <- app(mask_stack, fun = function(x) {
  if (all(!is.na(x))) 1 else NA
})

# Use the region where all key variables are non-NA as the sampling mask
land_mask <- valid_mask
cat("  Using complete-case environmental mask for pseudo-absence sampling\n")

# ============================================================================
# 2. Check and set CRS
# ============================================================================

# Ensure presence_spat has a CRS
if (is.na(crs(presence_spat)) || crs(presence_spat) == "") {
  crs(presence_spat) <- crs(land_mask)
  cat("  Set presence CRS to match climate data:", crs(land_mask), "\n")
}

# Ensure consistency
cat("  Presence CRS:", crs(presence_spat), "\n")
cat("  Land mask CRS:", crs(land_mask), "\n")

# ============================================================================
# 2.1 KDE helper (compatible with older R)
# ============================================================================

build_gaussian_kernel <- function(sigma_cells) {
  radius_cells <- max(1, ceiling(3 * sigma_cells))
  offsets <- -radius_cells:radius_cells
  kernel <- outer(offsets, offsets, function(x, y) exp(-(x^2 + y^2) / (2 * sigma_cells^2)))
  kernel / sum(kernel)
}

build_kde_weight_raster <- function(presence_spat, land_mask, max_dist_deg) {
  buffer_deg <- max_dist_deg + (3 * kde_bandwidth_km / 111)
  presence_extent <- ext(presence_spat)
  kde_extent <- ext(
    xmin(presence_extent) - buffer_deg,
    xmax(presence_extent) + buffer_deg,
    ymin(presence_extent) - buffer_deg,
    ymax(presence_extent) + buffer_deg
  )
  
  xmin(kde_extent) <- max(xmin(kde_extent), -180)
  xmax(kde_extent) <- min(xmax(kde_extent), 180)
  ymin(kde_extent) <- max(ymin(kde_extent), -90)
  ymax(kde_extent) <- min(ymax(kde_extent), 90)
  
  land_mask_crop <- crop(land_mask, kde_extent)
  native_res_km <- estimate_raster_resolution_km(land_mask_crop, presence_spat)$mean_km
  
  kde_base_raster <- land_mask_crop
  aggregation_factor <- 1
  
  if (!is.null(kde_compute_resolution_km) && is.finite(kde_compute_resolution_km) &&
      kde_compute_resolution_km > native_res_km) {
    aggregation_factor <- max(1, round(kde_compute_resolution_km / native_res_km))
    kde_base_raster <- aggregate(land_mask_crop, fact = aggregation_factor, fun = "mean", na.rm = TRUE)
    kde_base_raster[!is.finite(kde_base_raster)] <- NA
    kde_base_raster[!is.na(kde_base_raster)] <- 1
  }
  
  res_km <- estimate_raster_resolution_km(kde_base_raster, presence_spat)$mean_km
  sigma_cells <- max(kde_bandwidth_km / res_km, 1)
  kernel <- build_gaussian_kernel(sigma_cells)
  
  pres_raster <- rasterize(presence_spat, kde_base_raster, field = 1, fun = "sum", background = 0)
  kde_raster <- focal(pres_raster, w = kernel, fun = "sum", na.policy = "omit", fillvalue = 0)
  kde_raster[is.na(kde_base_raster)] <- NA
  
  dist_to_presence <- distance(kde_base_raster, y = presence_spat)
  kde_raster[dist_to_presence > (max_dist_deg * 111000)] <- NA
  
  list(
    weights = kde_raster,
    native_resolution_km = native_res_km,
    resolution_km = res_km,
    sigma_cells = sigma_cells,
    kde_extent = kde_extent,
    buffer_deg = buffer_deg,
    aggregation_factor = aggregation_factor
  )
}

filter_by_presence_distance <- function(candidate_xy, presence_spat, min_dist_meters) {
  if (nrow(candidate_xy) == 0) {
    return(candidate_xy[FALSE, , drop = FALSE])
  }
  
  candidate_sf <- sf::st_as_sf(
    data.frame(Longitude = candidate_xy[, 1], Latitude = candidate_xy[, 2]),
    coords = c("Longitude", "Latitude"),
    crs = target_crs
  )
  
  presence_xy <- terra::crds(presence_spat)
  presence_sf <- sf::st_as_sf(
    data.frame(Longitude = presence_xy[, 1], Latitude = presence_xy[, 2]),
    coords = c("Longitude", "Latitude"),
    crs = target_crs
  )
  
  dist_mat <- sf::st_distance(candidate_sf, presence_sf)
  min_dist <- apply(dist_mat, 1, min)
  keep_idx <- which(as.numeric(min_dist) >= min_dist_meters)
  
  candidate_xy[keep_idx, , drop = FALSE]
}

sample_kde_pseudo <- function(n_points, presence_spat, land_mask, kde_info, min_dist_meters) {
  weight_vals <- terra::values(kde_info$weights, mat = FALSE)
  available_cells <- which(!is.na(weight_vals))
  available_weights <- weight_vals[available_cells]
  available_weights[!is.finite(available_weights) | available_weights < 0] <- 0
  
  if (sum(available_weights) <= 0) {
    available_weights <- rep(1, length(available_cells))
  }
  
  selected_cells <- integer(0)
  iter <- 1
  
  while (length(selected_cells) < n_points && length(available_cells) > 0 && iter <= kde_max_iter) {
    remaining_need <- n_points - length(selected_cells)
    draw_n <- min(length(available_cells), max(remaining_need * kde_candidate_multiplier, remaining_need))
    
    sampled_idx <- sample(seq_along(available_cells), size = draw_n, replace = FALSE, prob = available_weights)
    sampled_cells <- available_cells[sampled_idx]
    sampled_xy <- terra::xyFromCell(kde_info$weights, sampled_cells)
    
    valid_xy <- filter_by_presence_distance(sampled_xy, presence_spat, min_dist_meters)
    
    if (nrow(valid_xy) > 0) {
      valid_cells <- terra::cellFromXY(kde_info$weights, valid_xy)
      selected_cells <- unique(c(selected_cells, valid_cells))
    }
    
    available_cells <- available_cells[-sampled_idx]
    available_weights <- available_weights[-sampled_idx]
    iter <- iter + 1
  }
  
  if (length(selected_cells) < n_points) {
    warning("KDE strategy could only generate ", length(selected_cells), " pseudo-absence points out of requested ", n_points)
  }
  
  selected_cells <- selected_cells[seq_len(min(length(selected_cells), n_points))]
  selected_xy <- terra::xyFromCell(kde_info$weights, selected_cells)
  
  if (kde_jitter_within_cell && nrow(selected_xy) > 0) {
    cell_res <- terra::res(kde_info$weights)
    jitter_x <- stats::runif(nrow(selected_xy), min = -0.5 * cell_res[1], max = 0.5 * cell_res[1])
    jitter_y <- stats::runif(nrow(selected_xy), min = -0.5 * cell_res[2], max = 0.5 * cell_res[2])
    selected_xy[, 1] <- selected_xy[, 1] + jitter_x
    selected_xy[, 2] <- selected_xy[, 2] + jitter_y
    
    jitter_df <- data.frame(Longitude = selected_xy[, 1], Latitude = selected_xy[, 2])
    jitter_vect <- vect(jitter_df, geom = c("Longitude", "Latitude"), crs = crs(land_mask))
    
    ext_vals <- ext(kde_info$weights)
    coords <- terra::crds(jitter_vect)
    inside_extent <- coords[, 1] >= xmin(ext_vals) & coords[, 1] <= xmax(ext_vals) &
      coords[, 2] >= ymin(ext_vals) & coords[, 2] <= ymax(ext_vals)
    
    if (any(!inside_extent)) {
      selected_xy[!inside_extent, ] <- terra::xyFromCell(kde_info$weights, selected_cells[!inside_extent])
    }
  }
  
  if (nrow(selected_xy) == 0) {
    stop("Error: KDE strategy failed to generate any pseudo-absence points. Try reducing pa_dist_min or increasing kde_candidate_multiplier.")
  }
  
  vect(
    data.frame(Longitude = selected_xy[, 1], Latitude = selected_xy[, 2]),
    geom = c("Longitude", "Latitude"),
    crs = crs(land_mask)
  )
}

# ============================================================================
# 3. Geographic filter on the sampling extent (drop polar regions)
# ============================================================================

cat("\n  Filtering land mask by latitude limits...\n")
cat("  Target Range: ", pa_lat_limit_south, " to ", pa_lat_limit_north, "\n")

# Get land_mask coordinates
# Logical mask: 1 inside the extent, NA outside
# Use terra::init or a simple coordinate check
lat_mask <- init(land_mask, "y") # raster whose values equal latitude

# Set regions outside the extent to NA
land_mask_filtered <- land_mask
land_mask_filtered[lat_mask > pa_lat_limit_north | lat_mask < pa_lat_limit_south] <- NA

# Update the variable for later functions
# Replace the original mask with the filtered one
land_mask <- land_mask_filtered

# Compute the filtered land-area fraction (optional, for logging)
non_na_count <- global(land_mask, "notNA")$notNA
cat("  Land mask updated. Available pixels for sampling:", non_na_count, "\n")

# Extent-visualization function: clip to the joint presence + pseudo-absence extent with a buffer
plot_buffer_deg <- 8

get_plot_extent <- function(presence_spat, pseudo_list, buffer_deg = 8) {
  all_extents <- c(list(ext(presence_spat)), lapply(pseudo_list, ext))
  
  xmin_val <- min(sapply(all_extents, xmin))
  xmax_val <- max(sapply(all_extents, xmax))
  ymin_val <- min(sapply(all_extents, ymin))
  ymax_val <- max(sapply(all_extents, ymax))
  
  plot_extent <- ext(
    xmin_val - buffer_deg,
    xmax_val + buffer_deg,
    ymin_val - buffer_deg,
    ymax_val + buffer_deg
  )
  
  xmin(plot_extent) <- max(xmin(plot_extent), -180)
  xmax(plot_extent) <- min(xmax(plot_extent), 180)
  ymin(plot_extent) <- max(ymin(plot_extent), -90)
  ymax(plot_extent) <- min(ymax(plot_extent), 90)
  
  plot_extent
}

# Core-parameter guard (use optimal defaults if config.R is not updated)
if (!exists("pa_strategy")) pa_strategy <- "disk"
if (!exists("pa_dist_min")) pa_dist_min <- 0.5   # degrees
if (!exists("pa_dist_max")) pa_dist_max <- 20    # degrees (critical: limits to ~2200 km)

cat("\n  PA Strategy Parameters:\n")
cat("  - Strategy:", pa_strategy, "\n")
cat("  - Min Distance:", pa_dist_min, "degrees (~", round(pa_dist_min * 111), "km)\n")
cat("  - Max Distance:", pa_dist_max, "degrees (~", round(pa_dist_max * 111), "km)\n")

# Unit conversion (biomod2 requires metres for distances under lon/lat)
min_dist_meters <- pa_dist_min * 111000
max_dist_meters <- pa_dist_max * 111000

kde_info <- NULL
if (pa_strategy == "kde") {
  cat("  - KDE Bandwidth:", kde_bandwidth_km, "km\n")
  cat("  - KDE Compute Resolution:", if (is.null(kde_compute_resolution_km)) "native" else paste0(kde_compute_resolution_km, " km"), "\n")
  cat("  - KDE Candidate Multiplier:", kde_candidate_multiplier, "\n")
  cat("  - KDE Max Iter:", kde_max_iter, "\n")
  cat("  - KDE Jitter Within Cell:", kde_jitter_within_cell, "\n")
  
  kde_info <- build_kde_weight_raster(presence_spat, land_mask, pa_dist_max)
  cat("  - Native Raster Resolution:", round(kde_info$native_resolution_km, 2), "km\n")
  cat("  - KDE Raster Resolution:", round(kde_info$resolution_km, 2), "km\n")
  cat("  - KDE Aggregation Factor:", kde_info$aggregation_factor, "\n")
  cat("  - KDE Sigma (cells):", round(kde_info$sigma_cells, 2), "\n")
  cat("  - KDE Crop Buffer:", round(kde_info$buffer_deg, 2), "degrees\n")
  cat("  - KDE Extent: ",
      round(xmin(kde_info$kde_extent), 2), " to ", round(xmax(kde_info$kde_extent), 2), ", ",
      round(ymin(kde_info$kde_extent), 2), " to ", round(ymax(kde_info$kde_extent), 2), "\n", sep = "")
}

# ============================================================================
# 4. Define the two schemes (read from config.R)
# ============================================================================

scheme1 <- pa_scheme1

# Scheme 2
scheme2 <- pa_scheme2

cat("\n  Scheme definitions (from config.R):\n")
cat("  - Scheme 1: 1:1 ratio\n")
cat("  - Scheme 2: 1000 points (1:3 ratio)\n")

# ============================================================================
# 5. Pseudo-absence generation function
# ============================================================================

generate_pseudo <- function(scheme, presence_spat, land_mask, kde_info = NULL) {
  cat("\n  Generating", scheme$name, "...\n")

  # Seed the RNG so pseudo-absence generation is reproducible across re-runs.
  # Without this, sample() and runif() jitter below draw from an unseeded
  # stream and every run gives different points.
  set.seed(seed_val)
  
  # Preprocess the count vector
  if (is.numeric(scheme$nb.absences)) {
    nb.vec <- scheme$nb.absences
    if (length(nb.vec) < scheme$nb.rep) nb.vec <- rep(nb.vec, scheme$nb.rep)
  } else {
    nb.vec <- rep(scheme$nb.absences, scheme$nb.rep)
  }
  
  results <- list()
  
  for (rep in 1:scheme$nb.rep) {
    curr_nb <- nb.vec[rep]
    # Handle the "equal" string
    if (is.character(curr_nb) && curr_nb == "equal") {
      curr_nb_val <- nrow(presence_spat)
    } else {
      curr_nb_val <- as.numeric(curr_nb)
    }
    
    cat("    Replicate", rep, "- Requesting", curr_nb_val, "points...\n")
    
    if (pa_strategy == "kde") {
      if (is.null(kde_info)) {
        stop("Error: KDE info not initialized.")
      }
      
      res_spat <- sample_kde_pseudo(
        n_points = curr_nb_val,
        presence_spat = presence_spat,
        land_mask = land_mask,
        kde_info = kde_info,
        min_dist_meters = min_dist_meters
      )
    } else {
      pa_args <- list(
        resp.var = presence_spat,
        expl.var = land_mask,
        nb.absences = curr_nb_val,
        nb.rep = 1,
        strategy = pa_strategy
      )
      
      if (pa_strategy == "disk") {
        pa_args$dist.min <- min_dist_meters
        pa_args$dist.max <- max_dist_meters
      } else if (pa_strategy == "sre") {
        pa_args$sre.quant <- sre_quant
      }
      
      pa_obj <- do.call(bm_PseudoAbsences, pa_args)
      pa_df <- as.data.frame(pa_obj)
      all_cols <- names(pa_df)
      x_col <- intersect(all_cols, c("x", "xy.x", "Longitude"))[1]
      y_col <- intersect(all_cols, c("y", "xy.y", "Latitude"))[1]
      resp_col <- all_cols[3]
      pa_select_col <- grep("PA", all_cols, value = TRUE)[1]
      pa_coords <- pa_df[pa_df[[pa_select_col]] == TRUE & is.na(pa_df[[resp_col]]), c(x_col, y_col)]
      
      if (nrow(pa_coords) == 0) {
        stop("Error: No pseudo-absence points were generated. Check dist.min or land_mask.")
      }
      
      res_spat <- vect(pa_coords, geom = c(x_col, y_col), crs = crs(land_mask))
    }
    
    res_spat$presence <- 0
    res_spat$rep <- rep
    
    results[[rep]] <- res_spat
    cat("      Successfully generated:", nrow(res_spat), "points\n")
  }
  return(results)
}

# ============================================================================
# 6. Generate both schemes
# ============================================================================

cat("\n  Generating pseudo absence points...\n")
cat("  Strategy:", pa_strategy, "\n")
cat("  Distance minimum:", pa_dist_min, "degree (≈", round(pa_dist_min * 111), "km)\n")

# Generate scheme 1
scheme1_results <- generate_pseudo(scheme1, presence_spat, land_mask, kde_info)

# Generate scheme 2
scheme2_results <- generate_pseudo(scheme2, presence_spat, land_mask, kde_info)

# ============================================================================
# 7. Save data
# ============================================================================

dir.create(output_dirs$pseudo_absence, showWarnings = FALSE, recursive = TRUE)

# Save scheme 1
cat("\n  Saving Scheme 1 (1:1)...\n")
pseudo1 <- scheme1_results[[1]]
pseudo1$pseudoabsence <- 0
pseudo1_csv <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1to1.csv")
write.csv(as.data.frame(pseudo1, geom = "XY"), pseudo1_csv, row.names = FALSE)
cat("    Saved:", pseudo1_csv, "\n")
pseudo1_rds <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1to1_spat.rds")
saveRDS(pseudo1, pseudo1_rds)

# Save scheme 2
cat("\n  Saving Scheme 2 (1000 points)...\n")
pseudo2 <- scheme2_results[[1]]
pseudo2$pseudoabsence <- 0
pseudo2_csv <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1000.csv")
write.csv(as.data.frame(pseudo2, geom = "XY"), pseudo2_csv, row.names = FALSE)
cat("    Saved:", pseudo2_csv, "\n")
pseudo2_rds <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1000_spat.rds")
saveRDS(pseudo2, pseudo2_rds)

# ============================================================================
# 8. Visualization check
# ============================================================================

cat("\n  Generating visualization plot...\n")

pdf_file <- file.path(output_dirs$pseudo_absence, "pseudoabsence_check.pdf")
pdf(pdf_file, width = 12, height = 10)

# Scheme 1 plot extent
plot_extent_1 <- get_plot_extent(presence_spat, list(pseudo1), plot_buffer_deg)
land_mask_plot_1 <- crop(land_mask, plot_extent_1)
cat("  Visualization extent - Scheme 1: ",
    round(xmin(plot_extent_1), 2), " to ", round(xmax(plot_extent_1), 2), ", ",
    round(ymin(plot_extent_1), 2), " to ", round(ymax(plot_extent_1), 2), "\n", sep = "")

# Page 1: scheme 1
plot(land_mask_plot_1, main = "Scheme 1: 1:1 Ratio (Zoomed View)",
     col = gray(0.7, alpha = 0.5))
points(presence_spat, pch = 21, col = "black", bg = "green", cex = 0.6)
points(pseudo1, pch = 23, col = "red", bg = "orange", cex = 0.8)
legend("bottomright",
       legend = c("Presence", "Pseudo Absence"),
       pch = c(21, 23),
       col = c("black", "red"),
       pt.bg = c("green", "orange"),
       bg = "white",
       cex = 0.8)

# Scheme 2 plot extent
plot_extent_2 <- get_plot_extent(presence_spat, list(pseudo2), plot_buffer_deg)
land_mask_plot_2 <- crop(land_mask, plot_extent_2)
cat("  Visualization extent - Scheme 2: ",
    round(xmin(plot_extent_2), 2), " to ", round(xmax(plot_extent_2), 2), ", ",
    round(ymin(plot_extent_2), 2), " to ", round(ymax(plot_extent_2), 2), "\n", sep = "")

# Page 2: scheme 2
plot(land_mask_plot_2, main = "Scheme 2: 1000 Points (1:3 Ratio) (Zoomed View)",
     col = gray(0.7, alpha = 0.5))
points(presence_spat, pch = 21, col = "black", bg = "green", cex = 0.6)
points(pseudo2, pch = 23, col = "red", bg = "orange", cex = 0.8)
legend("bottomright",
       legend = c("Presence", "Pseudo Absence"),
       pch = c(21, 23),
       col = c("black", "red"),
       pt.bg = c("green", "orange"),
       bg = "white",
       cex = 0.8)

dev.off()

cat("  Saved visualization plot:", pdf_file, "\n")

# ============================================================================
# 9. Save metadata
# ============================================================================

metadata <- list(
  presence_count = n_presence,
  strategy = pa_strategy,
  seed = seed_val,
  dist_min = pa_dist_min,  # degrees
  dist_min_km = round(pa_dist_min * 111),  # ~111 km
  dist_max = pa_dist_max,
  dist_max_km = round(pa_dist_max * 111),
  kde = if (pa_strategy == "kde") list(
    bandwidth_km = kde_bandwidth_km,
    compute_resolution_km = kde_compute_resolution_km,
    jitter_within_cell = kde_jitter_within_cell,
    candidate_multiplier = kde_candidate_multiplier,
    max_iter = kde_max_iter,
    native_resolution_km = kde_info$native_resolution_km,
    raster_resolution_km = kde_info$resolution_km,
    aggregation_factor = kde_info$aggregation_factor,
    sigma_cells = kde_info$sigma_cells,
    crop_buffer_deg = kde_info$buffer_deg,
    kde_extent = c(
      xmin = xmin(kde_info$kde_extent),
      xmax = xmax(kde_info$kde_extent),
      ymin = ymin(kde_info$kde_extent),
      ymax = ymax(kde_info$kde_extent)
    )
  ) else NULL,
  scheme1 = list(
    name = scheme1$name,
    nb.absences = scheme1$nb.absences,
    nb.rep = scheme1$nb.rep,
    pseudoabsence_count = nrow(pseudo1)
  ),
  scheme2 = list(
    name = scheme2$name,
    nb.absences = scheme2$nb.absences,
    nb.rep = scheme2$nb.rep,
    pseudoabsence_count = nrow(pseudo2)
  )
)

metadata_file <- file.path(output_dirs$pseudo_absence, "pseudoabsence_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  02 Pseudo Absence Generation - Completed\n")
cat("=================================================================\n")
cat("\n")
# Free memory
gc()
