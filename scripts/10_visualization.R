# ============================================================================
# 10_visualization.R - Result visualization
# ============================================================================
# Functionality:
#   - Visualize paleoclimate projection results
#   - Generate static figures and animation
#   - Compare fossil records with model predictions
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  10 Visualization - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Set output directory
# ============================================================================

viz_dir <- output_dirs$visualization
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(viz_dir, "static"), showWarnings = FALSE)
dir.create(file.path(viz_dir, "animation"), showWarnings = FALSE)
dir.create(file.path(viz_dir, "comparison"), showWarnings = FALSE)

cat("  Output directory:", viz_dir, "\n")

# ============================================================================
# 2. Load paleoclimate projection results
# ============================================================================

paleo_proj_dir <- output_dirs$paleoclimate

# Check output directory
if (!dir.exists(paleo_proj_dir)) {
  stop("Error: Paleoclimate projection output not found. Run 09 first.")
}

# Get all tif files
tif_files <- list.files(
  file.path(paleo_proj_dir, "projections"),
  pattern = "\\.tif$",
  full.names = TRUE
)

if (length(tif_files) == 0) {
  stop("Error: No projection TIF files found.")
}

# One period usually yields several ensemble rasters (e.g. EMwmeanByTSS / AUCroc / AUCprg / BOYCE).
# Keep one representative raster per period to avoid duplicate names causing layer-name clashes.
extract_period_name <- function(file_path, species_name) {
  pattern_to_remove <- paste0("_", species_name, ".*")
  gsub(pattern_to_remove, "", basename(file_path))
}

extract_projection_metric <- function(file_path) {
  metric <- sub(".*_EMwmeanBy([A-Za-z0-9]+)_.*", "\\1", basename(file_path))
  if (identical(metric, basename(file_path))) {
    return(NA_character_)
  }
  metric
}

select_projection_files <- function(files, species_name,
                                    metric_priority = c("TSS", "AUCroc", "AUCprg", "BOYCE")) {
  period_names_all <- vapply(files, extract_period_name, FUN.VALUE = character(1), species_name = species_name)
  metric_names_all <- vapply(files, extract_projection_metric, FUN.VALUE = character(1))
  
  selected_files <- character(0)
  
  for (period in unique(period_names_all)) {
    idx <- which(period_names_all == period)
    period_files <- files[idx]
    period_metrics <- metric_names_all[idx]
    
    chosen_idx <- NA_integer_
    for (metric_name in metric_priority) {
      hit <- which(period_metrics == metric_name)
      if (length(hit) > 0) {
        chosen_idx <- hit[1]
        break
      }
    }
    
    if (is.na(chosen_idx)) {
      chosen_idx <- 1
      cat("  Warning: No preferred metric found for", period,
          "- using", basename(period_files[chosen_idx]), "\n")
    }
    
    selected_files <- c(selected_files, period_files[chosen_idx])
  }
  
  selected_files
}

metric_priority <- unique(c(
  visualization_projection_metric,
  setdiff(c("TSS", "AUCroc", "AUCprg", "BOYCE"), visualization_projection_metric)
))

cat("  Visualization metric priority:", paste(metric_priority, collapse = " > "), "\n")
tif_files <- select_projection_files(
  tif_files,
  species_name = species_name,
  metric_priority = metric_priority
)
cat("  Using", length(tif_files), "projection files after metric filtering (one per period)\n")

format_ma_label <- function(ma_value) {
  if (!is.finite(ma_value)) {
    return(NA_character_)
  }
  if (ma_value < 1) {
    return(paste0(round(ma_value * 1000, 3), " ka"))
  }
  paste0(round(ma_value, 3), " Ma")
}

extract_period_ka <- function(period_name) {
  has_ka <- grepl("_(\\d+)k$", period_name)
  if (!has_ka) {
    return(NA_real_)
  }
  as.numeric(sub(".*_(\\d+)k$", "\\1", period_name))
}

build_period_metadata <- function(period_names, simulations_file_path = NULL) {
  period_labels <- setNames(period_names, period_names)
  period_ma_numeric <- setNames(seq_along(period_names), period_names)
  
  parsed_ka <- vapply(period_names, extract_period_ka, FUN.VALUE = numeric(1))
  parsed_idx <- which(is.finite(parsed_ka))
  if (length(parsed_idx) > 0) {
    period_labels[parsed_idx] <- paste0(parsed_ka[parsed_idx], " ka")
    period_ma_numeric[parsed_idx] <- parsed_ka[parsed_idx] / 1000
  }
  
  used_simulation_file <- FALSE
  if (!is.null(simulations_file_path) && file.exists(simulations_file_path)) {
    sim_info <- read.csv(simulations_file_path, stringsAsFactors = FALSE)
    if (all(c("Simulation.Name", "Ma") %in% names(sim_info))) {
      matched_periods <- intersect(period_names, sim_info$Simulation.Name)
      if (length(matched_periods) > 0) {
        matched_info <- sim_info[match(matched_periods, sim_info$Simulation.Name), , drop = FALSE]
        period_ma_numeric[matched_periods] <- matched_info$Ma
        period_labels[matched_periods] <- vapply(matched_info$Ma, format_ma_label, FUN.VALUE = character(1))
        used_simulation_file <- TRUE
      }
    }
  }
  
  list(
    period_labels = period_labels,
    period_ma_numeric = period_ma_numeric,
    used_simulation_file = used_simulation_file
  )
}

# --- Dynamically extract period names (no hardcoded "_Sargentodoxa") ---
period_names <- vapply(tif_files, extract_period_name, FUN.VALUE = character(1), species_name = species_name)
period_meta <- build_period_metadata(
  period_names = period_names,
  simulations_file_path = if (exists("simulations_file")) simulations_file else NULL
)
period_labels <- period_meta$period_labels
period_ma_numeric <- period_meta$period_ma_numeric

order_idx <- order(period_ma_numeric[period_names], decreasing = TRUE, na.last = TRUE)
tif_files <- tif_files[order_idx]
period_names <- period_names[order_idx]

if (isTRUE(period_meta$used_simulation_file)) {
  cat("  Ordered periods using simulations_file metadata\n")
} else if (all(is.finite(period_ma_numeric[period_names]))) {
  cat("  Ordered periods by parsed ages from period names\n")
} else {
  cat("  Warning: Could not fully infer period ages; using fallback order\n")
}

cat("  Found", length(tif_files), "projection files\n")

normalize_prediction_raster <- function(x) {
  xmax <- terra::global(x, "max", na.rm = TRUE)[1, 1]
  if (is.finite(xmax) && xmax > 1) {
    x <- x / 1000
  }
  x
}

suitability_palette <- c(
  "#f7f7f7", "#e9f2c9", "#cde48b", "#94c95b", "#4ea34a", "#1f6b2d"
)

# Load projections as a raster stack
proj_stack <- rast(tif_files)
proj_stack <- normalize_prediction_raster(proj_stack)

# --- Key fix: manually name proj_stack layers here ---
names(proj_stack) <- period_names
# --- end fix ---

cat("  Loaded projection stack:", nrow(proj_stack), "x", ncol(proj_stack), 
    "x", nlyr(proj_stack), "\n")

# (Optional) verify the names are set correctly
cat("  Verifying layer names. First 5 names in stack are:\n")
print(head(names(proj_stack)))


# Load the statistics summary
stats_dir <- file.path(paleo_proj_dir, "stats")
if (dir.exists(stats_dir)) {
  stats_files <- list.files(stats_dir, pattern = "summary\\.csv$", full.names = TRUE)
  stats_all <- lapply(stats_files, read.csv)
  names(stats_all) <- gsub("_summary\\.csv$", "", basename(stats_files))
  cat("  Loaded", length(stats_files), "stat files\n")
}

# ============================================================================
# 3. Load modern projection (for comparison)
# ============================================================================

modern_proj_file <- file.path(output_dirs$projection, "projection_modern_map.pdf")
if (file.exists(modern_proj_file)) {
  cat("  Modern projection found:", modern_proj_file, "\n")
}

# ============================================================================
# 4. Load fossil records (for validation)
# ============================================================================

# Check the fossil record file
fossil_data <- NULL
if (file.exists(fossil_record_file)) {
  fossil_data <- read.csv(fossil_record_file)
  # Check required columns: paleoLong, paleoLat, age
  required_cols <- c("paleoLong", "paleoLat", "age")
  if (all(required_cols %in% names(fossil_data))) {
    cat("  Loaded fossil records:", nrow(fossil_data), "records\n")
    cat("  Columns:", paste(required_cols, collapse = ", "), "\n")
    cat("  Available fossil ages:", paste(sort(unique(fossil_data$age)), collapse = ", "), "\n")
  } else {
    fossil_data <- NULL
    cat("  Fossil record file missing required columns:", required_cols, "\n")
  }
}

# ============================================================================
# 5. Visualization functions (full support for real Ma ages)
# ============================================================================

# 5.1 Single-period distribution map
plot_single_period <- function(r, period_name, output_file, 
                               period_label = period_name, # <-- new: accept real Ma label
                               add_fossil = FALSE, fossil_data = NULL,
                               threshold = 0.5, ...) {
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "value"
  
  theme_opts <- theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm")
  )
  
  p <- ggplot() +
    geom_tile(data = r_df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradientn(
      name = "Habitat suitability (0-1)",
      colors = suitability_palette,
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(crs = "EPSG:4326") +
    labs(
      title = paste(species_name, "-", period_label), # <-- key: title uses real Ma label
      subtitle = "Species Distribution Model Projection",
      x = "Longitude", y = "Latitude"
    ) +
    theme_opts
  
  if (add_fossil && !is.null(fossil_data)) {
    p <- p + geom_point(data = fossil_data, aes(x = lon, y = lat), color = "red", size = 2, shape = 1)
  }
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  return(p)
}

# 5.2 Multi-period comparison (memory-optimized)
plot_multi_period_comparison <- function(proj_stack, periods, output_file, 
                                         period_labels, ncol = 4, 
                                         downsample_factor = 3) { # <-- new: downsample parameter
  plot_data <- list()
  
  cat("    Extracting and downsampling raster layers (factor=", downsample_factor, ")...\n", sep="")
  
  for (i in seq_along(periods)) {
    if (periods[i] %in% names(proj_stack)) {
      r <- proj_stack[[periods[i]]]
      
      # --- Core memory optimization: lower spatial resolution ---
      # With factor = 3, data volume drops to 1/9!
      if (downsample_factor > 1) {
        r <- terra::aggregate(r, fact = downsample_factor, fun = mean, na.rm = TRUE)
      }
      # ----------------------------------------
      
      r_df <- as.data.frame(r, xy = TRUE)
      names(r_df)[3] <- "value"
      
      label <- period_labels[periods[i]]
      r_df$period_label <- factor(label, levels = period_labels[periods]) 
      plot_data[[i]] <- r_df
    }
  }
  
  if (length(plot_data) == 0) stop("No valid periods found")
  
  cat("    Merging data and rendering ggplot (this may still take a moment)...\n")
  combined_df <- do.call(rbind, plot_data)
  
  # Scale image height by panel count (~2.5 inch per row)
  n_panels <- length(periods)
  n_rows <- ceiling(n_panels / ncol)
  dynamic_height <- max(10, n_rows * 2.5) 
  
  p <- ggplot() +
    geom_tile(data = combined_df, aes(x = x, y = y, fill = value)) +
    facet_wrap(~ period_label, ncol = ncol) +
    scale_fill_gradientn(
      name = "Habitat suitability (0-1)",
      colors = suitability_palette,
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(crs = "EPSG:4326") +
    theme_minimal() + 
    theme(legend.position = "bottom", strip.background = element_rect(fill = "lightgray"))
  
  ggsave(output_file, p, width = 14, height = dynamic_height, dpi = 300)
  
  # --- Actively free memory ---
  rm(combined_df, plot_data)
  gc()
  # ----------------------
  
  return(p)
}

# 5.3 Habitat area trend (numeric real-age X axis)
plot_area_trend <- function(stats_all, period_order, period_ma_numeric, output_file) { # <-- new: numeric age dictionary
  
  valid_periods <- intersect(period_order, names(stats_all))
  ordered_stats <- stats_all[valid_periods]
  
  area_data <- data.frame(
    period = valid_periods,
    mean = sapply(ordered_stats, function(x) x$mean),
    sd = sapply(ordered_stats, function(x) x$sd)
  )
  
  # --- Key: map X axis to real numeric ages ---
  area_data$Ma <- period_ma_numeric[valid_periods]
  
  p <- ggplot(area_data, aes(x = Ma, y = mean)) +  # <-- use numeric Ma as X axis
    geom_line(color = "steelblue", linewidth = 1) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, fill = "steelblue") +
    geom_point(color = "darkblue", size = 2) +
    scale_x_reverse() + # --- key: reverse X axis (97.2 left, 0 right) ---
    labs(
      title = paste0(species_name, " - Habitat Suitability Trend"),
      subtitle = "Trend over geological time",
      x = "Age (Ma)", # <-- axis renamed to Age (Ma)
      y = "Mean Suitability Score"
    ) +
    theme_minimal()
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  return(p)
}

# 5.4 Single-period fossil comparison map
plot_fossil_comparison <- function(proj_stack, period_name, fossil_data, 
                                   period_ma, output_file, threshold = 0.5) {
  # Get the raster for the corresponding period
  if (!period_name %in% names(proj_stack)) {
    cat("    Warning: period", period_name, "not found in proj_stack\n")
    return(NULL)
  }
  r <- proj_stack[[period_name]]
  
  # Filter fossil records for this period
  fossil_subset <- fossil_data[fossil_data$age == period_ma, ]
  if (nrow(fossil_subset) == 0) {
    cat("    Warning: no fossil records for period", period_name, "at", period_ma, "Ma\n")
    return(NULL)
  }
  
  # Convert to SpatVector (using paleoLong, paleoLat)
  fossil_vect <- vect(fossil_subset, geom = c("paleoLong", "paleoLat"), crs = "EPSG:4326")
  
  # Convert to data frame for ggplot
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "value"
  
  # Get the period label
  period_label <- if (!is.null(period_labels[[period_name]])) {
    period_labels[[period_name]]
  } else {
    period_name
  }
  
  p <- ggplot() +
    geom_tile(data = r_df, aes(x = x, y = y, fill = value)) +
    geom_spatvector(data = fossil_vect, color = "red", size = 2, shape = 1) +
    scale_fill_gradientn(
      name = "Habitat suitability (0-1)",
      colors = suitability_palette,
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(crs = "EPSG:4326") +
    labs(
      title = paste0(species_name, " - ", period_label),
      subtitle = paste0("Fossil records at ", period_ma, " Ma (red circles)"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  cat("    Created:", basename(output_file), "\n")
  return(p)
}

# 5.5 Multi-period fossil comparison panel (final version - precise fossil matching)
plot_fossil_multi_comparison <- function(proj_stack, period_names, fossil_data,
                                         period_ma_numeric, output_file, 
                                         ncol = 5, downsample_factor = 4) {
  
  # --- 1. Data preparation ---
  
  # Prepare two lists: rasters and fossil data
  raster_data_list <- list()
  fossil_data_list <- list()
  
  # Keep only periods with fossil records
  valid_periods <- c()
  for (p in period_names) {
    ma_val <- period_ma_numeric[p]
    if (!is.na(ma_val) && any(fossil_data$age == ma_val, na.rm = TRUE)) {
      valid_periods <- c(valid_periods, p)
    }
  }
  
  if (length(valid_periods) == 0) {
    cat("    Warning: no periods with fossil records found\n")
    return(NULL)
  }
  
  cat("    Plotting", length(valid_periods), "periods with fossil records\n")
  
  # Loop over each valid period
  for (i in seq_along(valid_periods)) {
    p_name <- valid_periods[i]
    ma_val <- period_ma_numeric[p_name]
    label <- period_labels[p_name]
    
    # Prepare raster data
    r <- proj_stack[[p_name]]
    if (downsample_factor > 1) {
      r <- aggregate(r, fact = downsample_factor)
    }
    r_df <- as.data.frame(r, xy = TRUE)
    names(r_df)[3] <- "value"
    r_df$period_label <- factor(label, levels = period_labels[valid_periods])
    raster_data_list[[i]] <- r_df
    
    # --- Key: prepare fossil data for the corresponding period ---
    fossil_subset <- fossil_data[fossil_data$age == ma_val, ]
    if(nrow(fossil_subset) > 0) {
      fossil_subset$period_label <- factor(label, levels = period_labels[valid_periods])
      fossil_data_list[[i]] <- fossil_subset
    }
  }
  
  # Merge into two large data frames
  combined_raster_df <- do.call(rbind, raster_data_list)
  combined_fossil_df <- do.call(rbind, fossil_data_list)
  
  # --- 2. Plotting ---
  
  p <- ggplot() +
    # Layer 1: downsampled raster base map
    geom_tile(data = combined_raster_df, aes(x = x, y = y, fill = value)) +
    
    # Layer 2: precisely matched fossil points
    geom_point(data = combined_fossil_df, aes(x = paleoLong, y = paleoLat), 
               color = "red", size = 2, shape = 1, stroke = 1.2) +
    
    # Layer 3: facet panels
    facet_wrap(~ period_label, ncol = ncol) +
    
    # Layer 4: styling
    scale_fill_gradientn(
      name = "Habitat suitability (0-1)",
      colors = suitability_palette,
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(crs = "EPSG:4326") +
    labs(
      title = paste0(species_name, " - Prediction vs Fossil Records"),
      subtitle = "Multiple paleoclimate periods with fossil occurrences",
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "lightgray"),
      plot.background = element_blank(),
      panel.background = element_blank()
    )
  
  # Scale image height by panel count
  dynamic_height <- max(8, 2.5 * ceiling(length(valid_periods) / ncol))
  
  ggsave(output_file, p, width = 16, height = dynamic_height, dpi = 300)
  cat("    Created:", basename(output_file), "\n")
  return(p)
}

# 5.6 Fossil-point suitability bar chart
plot_fossil_suitability_bars <- function(proj_stack, fossil_data, period_ma_numeric,
                                         period_labels, output_file) {
  fossil_records <- list()
  
  for (period_name in names(period_ma_numeric)) {
    if (!period_name %in% names(proj_stack)) next
    
    ma_val <- period_ma_numeric[[period_name]]
    fossil_subset <- fossil_data[fossil_data$age == ma_val, , drop = FALSE]
    if (nrow(fossil_subset) == 0) next
    
    fossil_vect <- vect(fossil_subset, geom = c("paleoLong", "paleoLat"), crs = "EPSG:4326")
    extracted_vals <- terra::extract(proj_stack[[period_name]], fossil_vect)[, 2]
    
    fossil_subset$suitability <- extracted_vals
    fossil_subset$period_name <- period_name
    fossil_subset$period_label <- unname(period_labels[[period_name]])
    fossil_subset$fossil_id <- paste0(fossil_subset$period_label, "_", seq_len(nrow(fossil_subset)))
    fossil_records[[length(fossil_records) + 1]] <- fossil_subset
  }
  
  if (length(fossil_records) == 0) {
    return(NULL)
  }
  
  fossil_plot_df <- do.call(rbind, fossil_records)
  fossil_plot_df <- fossil_plot_df[is.finite(fossil_plot_df$suitability), , drop = FALSE]
  if (nrow(fossil_plot_df) == 0) {
    return(NULL)
  }
  
  fossil_plot_df <- fossil_plot_df[order(fossil_plot_df$age, decreasing = TRUE), , drop = FALSE]
  fossil_plot_df$fossil_id <- factor(fossil_plot_df$fossil_id, levels = fossil_plot_df$fossil_id)
  
  p <- ggplot(fossil_plot_df, aes(x = fossil_id, y = suitability, fill = period_label)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0(species_name, " - Fossil Point Suitability"),
      subtitle = "Predicted suitability extracted at fossil occurrence coordinates",
      x = "Fossil record",
      y = "Habitat suitability (0-1)",
      fill = "Period"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  ggsave(output_file, p, width = 12, height = max(6, 0.28 * nrow(fossil_plot_df)), dpi = 300)
  cat("    Created:", basename(output_file), "\n")
  return(p)
}

# 5.7 Fossil vs background suitability distribution
plot_fossil_vs_background <- function(proj_stack, fossil_data, period_ma_numeric,
                                      period_labels, output_file,
                                      background_n = 5000, seed_value = 123) {
  comparison_list <- list()
  
  for (period_name in names(period_ma_numeric)) {
    if (!period_name %in% names(proj_stack)) next
    
    ma_val <- period_ma_numeric[[period_name]]
    fossil_subset <- fossil_data[fossil_data$age == ma_val, , drop = FALSE]
    if (nrow(fossil_subset) == 0) next
    
    r <- proj_stack[[period_name]]
    fossil_vect <- vect(fossil_subset, geom = c("paleoLong", "paleoLat"), crs = "EPSG:4326")
    fossil_vals <- terra::extract(r, fossil_vect)[, 2]
    fossil_vals <- fossil_vals[is.finite(fossil_vals)]
    if (length(fossil_vals) == 0) next
    
    background_vals <- terra::values(r, mat = FALSE)
    background_vals <- background_vals[is.finite(background_vals)]
    if (length(background_vals) == 0) next
    
    set.seed(seed_value)
    background_vals <- sample(
      background_vals,
      size = min(background_n, length(background_vals)),
      replace = FALSE
    )
    
    period_label <- unname(period_labels[[period_name]])
    comparison_list[[length(comparison_list) + 1]] <- data.frame(
      period = period_label,
      group = "Background",
      suitability = background_vals,
      stringsAsFactors = FALSE
    )
    comparison_list[[length(comparison_list) + 1]] <- data.frame(
      period = period_label,
      group = "Fossil",
      suitability = fossil_vals,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(comparison_list) == 0) {
    return(NULL)
  }
  
  comparison_df <- do.call(rbind, comparison_list)
  if (nrow(comparison_df) == 0) {
    return(NULL)
  }
  
  period_order <- unname(period_labels[names(period_ma_numeric)])
  period_order <- unique(period_order[period_order %in% comparison_df$period])
  comparison_df$period <- factor(comparison_df$period, levels = period_order)
  
  background_df <- comparison_df[comparison_df$group == "Background", , drop = FALSE]
  fossil_df <- comparison_df[comparison_df$group == "Fossil", , drop = FALSE]
  
  fossil_counts <- aggregate(suitability ~ period, data = fossil_df, FUN = length)
  names(fossil_counts)[2] <- "n_fossil"
  fossil_counts$label <- paste0("n = ", fossil_counts$n_fossil)
  
  fossil_point_layer <- geom_point(
    data = fossil_df,
    aes(x = period, y = suitability, colour = group),
    size = 2.4,
    alpha = 0.9,
    position = position_nudge(x = 0)
  )
  
  p <- ggplot() +
    geom_violin(
      data = background_df,
      aes(x = period, y = suitability, fill = group),
      width = 0.9,
      alpha = 0.5,
      colour = NA,
      trim = TRUE
    ) +
    geom_boxplot(
      data = background_df,
      aes(x = period, y = suitability),
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      colour = "grey30"
    ) +
    fossil_point_layer +
    geom_text(
      data = fossil_counts,
      aes(x = period, y = 0.98, label = label),
      inherit.aes = FALSE,
      size = 3.5,
      vjust = 1
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("Background" = "grey70")) +
    scale_colour_manual(values = c("Fossil" = "#c43c39")) +
    labs(
      title = paste0(species_name, " - Fossil vs Background Suitability"),
      subtitle = "Background suitability distribution with fossil-point values overlaid",
      x = "Period",
      y = "Habitat suitability (0-1)",
      colour = NULL,
      fill = "Background"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(output_file, p, width = max(10, 1.6 * length(period_order)), height = 6.5, dpi = 300)
  cat("    Created:", basename(output_file), "\n")
  return(p)
}

# ============================================================================
# 6. Generate visualizations
# ============================================================================

cat("\n  Generating visualizations...\n")

# 6.1 Single-period distribution maps
cat("  Generating high-resolution static maps for all periods...\n")
for (period in period_names) {
  output_file <- file.path(viz_dir, "static", paste0(period, "_distribution.png"))
  if (!file.exists(output_file)) {
    if (period %in% names(proj_stack)) {
      r <- proj_stack[[period]]
      # --- Pass the corresponding period age label ---
      plot_single_period(r, period, output_file, period_label = period_labels[period])
      cat("    Created:", basename(output_file), "\n")
    }
  }
}
cat("  All static maps are up to date.\n")

# 6.2 Multi-period comparison (all periods)
if (length(period_names) >= 4) {
  output_file <- file.path(viz_dir, "comparison", "multi_period_comparison.png")
  
  cat("  Generating multi-period comparison map...\n")
  
  # downsample_factor = 3 or 4 saves a lot of memory
  plot_multi_period_comparison(proj_stack, period_names[1:length(period_names)], 
                               output_file, ncol = 4, 
                               period_labels = period_labels,
                               downsample_factor = 3) # <-- tune to balance speed vs clarity
  
  cat("    Created: multi_period_comparison.png\n")
}

# 6.3 Habitat area trend
if (exists("stats_all")) {
  output_file <- file.path(viz_dir, "static", "area_trend.png")
  # --- Pass the extracted numeric age dictionary ---
  plot_area_trend(stats_all, period_names, period_ma_numeric, output_file)
  cat("    Created: area_trend.png\n")
}

# 6.4 Compare with fossil records
if (!is.null(fossil_data) && "paleoLong" %in% names(fossil_data) && "paleoLat" %in% names(fossil_data)) {
  cat("\n  Generating fossil comparison plots...\n")
  
  # 6.4.1 Single-period fossil comparison
  for (period in period_names) {
    ma_val <- period_ma_numeric[period]
    if (is.na(ma_val)) next
    
    # Check whether this period has fossil records
    fossil_subset <- fossil_data[fossil_data$age == ma_val, ]
    if (nrow(fossil_subset) == 0) next
    
    output_file <- file.path(viz_dir, "comparison", paste0("fossil_comparison_", period, ".png"))
    if (!file.exists(output_file)) {
      plot_fossil_comparison(proj_stack, period, fossil_data, ma_val, output_file)
    }
  }
  
  # 6.4.2 Multi-period fossil comparison panel
  output_file_multi <- file.path(viz_dir, "comparison", "fossil_multi_period_comparison.png")
  if (!file.exists(output_file_multi)) {
    plot_fossil_multi_comparison(proj_stack, period_names, fossil_data,
                                 period_ma_numeric, output_file_multi,
                                 ncol = 3, downsample_factor = 4)
  }
  
  output_file_bars <- file.path(viz_dir, "comparison", "fossil_point_suitability.png")
  if (!file.exists(output_file_bars)) {
    plot_fossil_suitability_bars(
      proj_stack = proj_stack,
      fossil_data = fossil_data,
      period_ma_numeric = period_ma_numeric,
      period_labels = period_labels,
      output_file = output_file_bars
    )
  }
  
  output_file_bg <- file.path(viz_dir, "comparison", "fossil_vs_background_suitability.png")
  if (!file.exists(output_file_bg)) {
    plot_fossil_vs_background(
      proj_stack = proj_stack,
      fossil_data = fossil_data,
      period_ma_numeric = period_ma_numeric,
      period_labels = period_labels,
      output_file = output_file_bg
    )
  }
}

# ============================================================================
# 7. Generate GIF animation (downscaled from high-res files)
# ============================================================================
cat("\n  Generating animation from high-resolution files...\n")

if (require(magick)) {
  
  # 1. Collect all high-res PNG paths in the right order
  high_res_pngs <- file.path(viz_dir, "static", paste0(period_names, "_distribution.png"))
  high_res_pngs <- high_res_pngs[file.exists(high_res_pngs)]
  
  if (length(high_res_pngs) > 0) {
    
    cat("  Found", length(high_res_pngs), "images to process.\n")
    cat("  Resizing images for animation (this may take a moment)...\n")
    
    # 2. Read each high-res file, downscale, and store in a list
    #    image_scale "800x" scales width proportionally to 800 px
    images_resized <- image_read(high_res_pngs) %>%
      image_scale("600x")
    
    # Print new size for verification
    cat("  Images resized to:\n")
    print(image_info(images_resized)[, 1:3])
    
    # 3. Combine the downscaled images into an animation
    animation <- image_animate(images_resized, fps = 2) # fps=2 means 2 frames per second
    
    # 4. Save the GIF
    output_file <- file.path(viz_dir, "animation", "paleoclimate_timeline.gif")
    image_write(animation, output_file)
    cat("    Successfully created animation:", basename(output_file), "\n")
    
    # 5. Free memory (important!)
    rm(images_resized, animation)
    gc()
    
  } else {
    cat("  Warning: No static PNG files found to create animation.\n")
  }
} else {
  cat("  Note: magick package not available. Skipping animation generation.\n")
  cat("  To enable animation, run: install.packages('magick')\n")
}

# ============================================================================
# 8. Save metadata
# ============================================================================

viz_metadata <- list(
  created = Sys.time(),
  n_periods = length(period_names),
  periods = period_names,
  output_dir = viz_dir,
  static_dir = file.path(viz_dir, "static"),
  animation_dir = file.path(viz_dir, "animation"),
  comparison_dir = file.path(viz_dir, "comparison")
)

saveRDS(viz_metadata, file.path(viz_dir, "visualization_metadata.rds"))

cat("\n")
cat("=================================================================\n")
cat("  10 Visualization - Completed\n")
cat("=================================================================\n")
cat("  Output directory:", viz_dir, "\n")
cat("\n")
