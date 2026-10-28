# ============================================================================
# 04_data_formatting.R - Data formatting for biomod2
# ============================================================================
# Purpose:
#   - Read presence points, pseudo-absence points, and selected climate variables
#   - Merge data and split into train/test sets
#   - Format into biomod2-compatible format
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  04 Data Formatting - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Read data
# ============================================================================

# Read presence points as SpatVector
presence_file <- file.path(output_dirs$data_prep, "sargentodoxa_presence_spat.rds")
if (!file.exists(presence_file)) {
  stop("Error: Presence data not found. Run 01_data_preparation.R first.")
}

presence_spat <- readRDS(presence_file)
n_presence <- nrow(presence_spat)

cat("  Loading presence data:", presence_file, "\n")
cat("  Presence points:", n_presence, "\n")

# Select pseudo-absence file based on config
cat("\n  Pseudo-absence scheme:", pa_scheme_to_use, "\n")

if (pa_scheme_to_use == "scheme_1to1") {
  pseudo_file <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1to1_spat.rds")
  cat("  Loading 1:1 pseudo-absence data...\n")
} else if (pa_scheme_to_use == "scheme_1000") {
  pseudo_file <- file.path(output_dirs$pseudo_absence, "sargentodoxa_pseudoabsence_1000_spat.rds")
  cat("  Loading 1000 points pseudo-absence data...\n")
} else {
  stop("Error: Unknown pa_scheme_to_use: ", pa_scheme_to_use)
}

if (!file.exists(pseudo_file)) {
  stop("Error: Pseudo-absence data not found: ", pseudo_file)
}

pseudo_spat <- readRDS(pseudo_file)
n_pseudo <- nrow(pseudo_spat)

cat("  Loaded:", pseudo_file, "\n")
cat("  Pseudo-absence points:", n_pseudo, "\n")

# Read selected variables list
var_file <- file.path(output_dirs$variable_select, "variables_selected.csv")
if (!file.exists(var_file)) {
  stop("Error: Selected variables not found. Run 03_variable_selection.R first.")
}

selected_vars <- read.csv(var_file, stringsAsFactors = FALSE)
selected_vars <- selected_vars$variable

cat("\n  Loading selected variables:", var_file, "\n")
cat("  Selected variables:", length(selected_vars), "\n")
cat("  ", paste(selected_vars, collapse = ", "), "\n")

# Read modern climate data
cat("\n  Loading modern climate data from:", modern_clim_dir, "\n")
clim_files <- list.files(modern_clim_dir, pattern = modern_clim_pattern, full.names = TRUE)

if (length(clim_files) == 0) {
  stop("Error: No climate files found in ", modern_clim_dir)
}

modern_clim <- rast(clim_files)
modern_clim <- modern_clim[[mixedsort(names(modern_clim))]]

names(modern_clim) <- paste0("Bio", 1:nlyr(modern_clim)) 
# Keep only selected variables
modern_clim_selected <- modern_clim[[selected_vars]]

cat("  Climate data loaded:", nlyr(modern_clim_selected), "variables\n")

# ============================================================================
# 2. Merge full dataset (Scheme B: no random train/test split)
# ============================================================================

cat("\n  Merging full dataset for spatial CV...\n")

presence_data_clean <- presence_spat[, FALSE]
presence_data_clean$presence <- 1

pseudo_data_clean <- pseudo_spat[, FALSE]
pseudo_data_clean$presence <- 0

all_points <- rbind(presence_data_clean, pseudo_data_clean)

cat("  Total points:", nrow(all_points), "\n")
cat("  - Presence:", n_presence, "\n")
cat("  - Pseudo-absence:", n_pseudo, "\n")

# ============================================================================
# 3. Filter points with missing predictor values, then build custom spatial CV folds
# ============================================================================

cat("\n  Filtering points with missing predictor values...\n")

all_point_env <- terra::extract(modern_clim_selected, all_points)
all_point_env <- as.data.frame(all_point_env)
if ("ID" %in% names(all_point_env)) {
  all_point_env$ID <- NULL
}

valid_env_rows <- complete.cases(all_point_env)
removed_na_points <- sum(!valid_env_rows)

if (removed_na_points > 0) {
  cat("  - Removing", removed_na_points, "points with NA predictor values\n")
  all_points <- all_points[valid_env_rows, ]
} else {
  cat("  - No points removed for predictor NAs\n")
}

n_presence_filtered <- sum(all_points$presence == 1)
n_pseudo_filtered <- sum(all_points$presence == 0)

cat("  - Remaining total points:", nrow(all_points), "\n")
cat("  - Remaining presence:", n_presence_filtered, "\n")
cat("  - Remaining pseudo-absence:", n_pseudo_filtered, "\n")

# ============================================================================
# 4. Build custom spatial CV folds
# ============================================================================

cat("\n  Building custom spatial CV table...\n")
cat("  - Method:", spatial_cv_method, "\n")
cat("  - Folds:", cv_k, "\n")
cat("  - Block size:", spatial_cv_block_size_km, "km\n")
cat("  - Projection CRS:", spatial_cv_projection_crs, "\n")

build_spatial_cv_table <- function(points_spat,
                                   k_folds,
                                   n_reps,
                                   block_size_km,
                                   projection_crs,
                                   retry_enabled = TRUE,
                                   retry_factor = 1.5,
                                   max_tries = 4) {
  if (k_folds < 2) {
    stop("Error: spatial CV requires at least 2 folds.")
  }
  if (n_reps < 1) {
    stop("Error: spatial CV requires at least 1 repetition.")
  }
  
  current_block_size_km <- block_size_km
  attempt <- 1
  
  while (attempt <= max_tries) {
    points_proj <- terra::project(points_spat, projection_crs)
    coords_proj <- terra::crds(points_proj)
    block_size_m <- current_block_size_km * 1000
    
    x_origin <- min(coords_proj[, 1], na.rm = TRUE)
    y_origin <- min(coords_proj[, 2], na.rm = TRUE)
    
    cv_table <- as.data.frame(matrix(FALSE, nrow = nrow(points_spat), ncol = k_folds * n_reps))
    fold_summary_list <- vector("list", n_reps)
    block_lookup_list <- vector("list", n_reps)
    fold_ids_list <- vector("list", n_reps)
    valid_folds <- TRUE
    run_counter <- 1L
    
    for (rep_id in seq_len(n_reps)) {
      shift_fraction <- (rep_id - 1) / n_reps
      x_origin_rep <- x_origin + (shift_fraction * block_size_m)
      y_origin_rep <- y_origin + (shift_fraction * block_size_m)
      
      block_col <- floor((coords_proj[, 1] - x_origin_rep) / block_size_m)
      block_row <- floor((coords_proj[, 2] - y_origin_rep) / block_size_m)
      block_id <- paste(block_col, block_row, sep = "_")
      
      unique_blocks <- unique(block_id)
      block_lookup <- data.frame(
        rep = rep_id,
        block_id = unique_blocks,
        stringsAsFactors = FALSE
      )
      block_lookup$col <- as.integer(vapply(strsplit(block_lookup$block_id, "_"), `[`, 1, FUN.VALUE = character(1)))
      block_lookup$row <- as.integer(vapply(strsplit(block_lookup$block_id, "_"), `[`, 2, FUN.VALUE = character(1)))
      block_lookup$fold_id <- ((block_lookup$col + block_lookup$row) %% k_folds) + 1L
      block_lookup$x_shift_fraction <- shift_fraction
      block_lookup$y_shift_fraction <- shift_fraction
      
      point_folds <- block_lookup$fold_id[match(block_id, block_lookup$block_id)]
      fold_ids_list[[rep_id]] <- point_folds
      block_lookup_list[[rep_id]] <- block_lookup
      
      fold_summary_rep <- dplyr::bind_rows(lapply(seq_len(k_folds), function(fold_id_value) {
        idx <- which(point_folds == fold_id_value)
        data.frame(
          rep = rep_id,
          fold = fold_id_value,
          run = run_counter + fold_id_value - 1L,
          x_shift_fraction = shift_fraction,
          y_shift_fraction = shift_fraction,
          n_total = length(idx),
          n_presence = sum(points_spat$presence[idx] == 1),
          n_pseudo = sum(points_spat$presence[idx] == 0)
        )
      }))
      
      if (!(all(fold_summary_rep$n_presence > 0) &&
            all(fold_summary_rep$n_pseudo > 0) &&
            all(fold_summary_rep$n_total > 0))) {
        valid_folds <- FALSE
        break
      }
      
      for (fold_id_value in seq_len(k_folds)) {
        cv_table[[run_counter]] <- point_folds != fold_id_value
        colnames(cv_table)[run_counter] <- paste0("_allData_RUN", run_counter)
        run_counter <- run_counter + 1L
      }
      
      fold_summary_list[[rep_id]] <- fold_summary_rep
    }
    
    if (valid_folds) {
      cv_table$`_allData_allRun` <- TRUE
      fold_summary <- dplyr::bind_rows(fold_summary_list)
      block_lookup <- dplyr::bind_rows(block_lookup_list)
      
      return(list(
        cv_table = cv_table,
        fold_ids = fold_ids_list,
        fold_summary = fold_summary,
        block_lookup = block_lookup,
        block_size_km = current_block_size_km,
        projection_crs = projection_crs,
        method = spatial_cv_method,
        attempt = attempt,
        n_reps = n_reps
      ))
    }
    
    if (!retry_enabled || attempt == max_tries) {
      stop(
        "Error: Spatial CV folds are imbalanced after ", attempt,
        " attempt(s). Try reducing cv_k or spatial_cv_block_size_km."
      )
    }
    
    current_block_size_km <- current_block_size_km * retry_factor
    attempt <- attempt + 1
  }
}

spatial_cv <- build_spatial_cv_table(
  points_spat = all_points,
  k_folds = cv_k,
  n_reps = cv_nb_rep,
  block_size_km = spatial_cv_block_size_km,
  projection_crs = spatial_cv_projection_crs,
  retry_enabled = spatial_cv_auto_retry,
  retry_factor = spatial_cv_retry_factor,
  max_tries = spatial_cv_max_tries
)

cv_table <- spatial_cv$cv_table
all_points$cv_fold <- spatial_cv$fold_ids[[1]]

cat("  Spatial CV table built successfully\n")
cat("  - Final block size:", spatial_cv$block_size_km, "km\n")
cat("  - Repetitions:", spatial_cv$n_reps, "\n")
cat("  - Total runs:", cv_k * spatial_cv$n_reps, "\n")
cat("  - Attempts:", spatial_cv$attempt, "\n")
for (i in seq_len(nrow(spatial_cv$fold_summary))) {
  cat(
    "  - Rep", spatial_cv$fold_summary$rep[i],
    "Fold", spatial_cv$fold_summary$fold[i],
    "(RUN", spatial_cv$fold_summary$run[i], "):",
    spatial_cv$fold_summary$n_total[i], "points (presence =",
    spatial_cv$fold_summary$n_presence[i], ", pseudo =",
    spatial_cv$fold_summary$n_pseudo[i], ")\n"
  )
}

# ============================================================================
# 5. Format for biomod2
# ============================================================================

cat("\n  Formatting full dataset for biomod2...\n")

models_dir <- "./Models"
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

formatted_data <- BIOMOD_FormatingData(
  resp.name = species_name,
  resp.var = all_points$presence,
  expl.var = modern_clim_selected,
  resp.xy = terra::crds(all_points),
  dir.name = models_dir,
  PA.nb.rep = 0
)

cat("  Data formatted successfully\n")
cat("  Formatted data object saved in:", models_dir, "\n")

# ============================================================================
# 6. Save results
# ============================================================================

cat("\n  Saving results...\n")

dir.create(output_dirs$formatted, showWarnings = FALSE, recursive = TRUE)

formatted_file <- file.path(output_dirs$formatted, "biomod2_formatted_data.RData")
save(formatted_data, file = formatted_file)
cat("  Saved formatted data:", formatted_file, "\n")

all_points_file <- file.path(output_dirs$formatted, "all_points.csv")
write.csv(as.data.frame(all_points, geom = "XY"), all_points_file, row.names = FALSE)
cat("  Saved all points:", all_points_file, "\n")

cv_table_file <- file.path(output_dirs$formatted, "spatial_cv_table.csv")
write.csv(cv_table, cv_table_file, row.names = FALSE)
cat("  Saved spatial CV table:", cv_table_file, "\n")

cv_object_file <- file.path(output_dirs$formatted, "spatial_cv_info.rds")
saveRDS(spatial_cv, cv_object_file)
cat("  Saved spatial CV info:", cv_object_file, "\n")

fold_summary_file <- file.path(output_dirs$formatted, "spatial_cv_fold_summary.csv")
write.csv(spatial_cv$fold_summary, fold_summary_file, row.names = FALSE)
cat("  Saved fold summary:", fold_summary_file, "\n")

# ============================================================================
# 7. Visualization
# ============================================================================

cat("\n  Generating visualization...\n")

pdf_file <- file.path(output_dirs$formatted, "spatial_cv_summary.pdf")
pdf(pdf_file, width = 12, height = 10)

fold_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

for (rep_id in seq_len(spatial_cv$n_reps)) {
  all_points$cv_fold <- spatial_cv$fold_ids[[rep_id]]
  point_colors <- fold_colors[((all_points$cv_fold - 1) %% length(fold_colors)) + 1]
  
  plot(
    modern_clim_selected[[1]],
    main = paste(
      "All Points by Spatial Fold - REP", rep_id,
      "(block =", spatial_cv$block_size_km, "km)"
    ),
    col = gray(0.7, alpha = 0.5)
  )
  points(all_points[all_points$presence == 0, ], pch = 1, col = point_colors[all_points$presence == 0], cex = 0.6)
  points(all_points[all_points$presence == 1, ], pch = 16, col = point_colors[all_points$presence == 1], cex = 0.7)
  legend(
    "bottomright",
    legend = paste("Fold", seq_len(cv_k)),
    col = fold_colors[seq_len(cv_k)],
    pch = 16,
    bg = "white",
    cex = 0.8
  )
}

plot.new()
text(0.5, 0.92, "Spatial CV Summary", cex = 1.5, font = 2)
text(0.5, 0.80, paste("Total points:", nrow(all_points)), cex = 1.2)
text(0.5, 0.72, paste("Presence points:", n_presence_filtered), cex = 1.1)
text(0.5, 0.64, paste("Pseudo-absence points:", n_pseudo_filtered), cex = 1.1)
text(0.5, 0.56, paste("Variables:", length(selected_vars)), cex = 1.1)
text(0.5, 0.48, paste("Spatial CV method:", spatial_cv$method), cex = 1.1)
text(0.5, 0.40, paste("Folds per repetition:", cv_k), cex = 1.1)
text(0.5, 0.32, paste("Repetitions:", spatial_cv$n_reps, " (total runs =", cv_k * spatial_cv$n_reps, ")"), cex = 1.1)
text(0.5, 0.24, paste("Final block size:", spatial_cv$block_size_km, "km"), cex = 1.1)
text(0.5, 0.16, paste("Projection CRS:", spatial_cv$projection_crs), cex = 1.1)
text(0.5, 0.08, paste("Filtered NA points:", removed_na_points, "; Attempts:", spatial_cv$attempt), cex = 1.1)

plot.new()
text(0.5, 0.94, "Spatial CV Fold Balance", cex = 1.4, font = 2)
rep_levels <- unique(spatial_cv$fold_summary$rep)
rep_y <- seq(0.82, 0.18, length.out = length(rep_levels))
for (idx in seq_along(rep_levels)) {
  rep_id <- rep_levels[idx]
  rep_rows <- spatial_cv$fold_summary[spatial_cv$fold_summary$rep == rep_id, , drop = FALSE]
  rep_text <- paste(
    paste0(
      "F", rep_rows$fold,
      ": n=", rep_rows$n_total,
      " (P=", rep_rows$n_presence,
      ", A=", rep_rows$n_pseudo, ")"
    ),
    collapse = "   "
  )
  text(0.5, rep_y[idx], paste("REP", rep_id, "-", rep_text), cex = 0.9)
}

dev.off()

cat("  Saved visualization:", pdf_file, "\n")

# ============================================================================
# 8. Save metadata
# ============================================================================

metadata <- list(
  pa_scheme = pa_scheme_to_use,
  n_presence = n_presence_filtered,
  n_pseudo = n_pseudo_filtered,
  n_total = nrow(all_points),
  removed_na_points = removed_na_points,
  n_variables = length(selected_vars),
  variables = selected_vars,
  seed = seed_val,
  cv_strategy = cv_strategy,
  cv_nb_rep = spatial_cv$n_reps,
  spatial_cv_method = spatial_cv$method,
  spatial_cv_block_size_km = spatial_cv$block_size_km,
  spatial_cv_projection_crs = spatial_cv$projection_crs,
  cv_k = cv_k,
  cv_table_file = cv_table_file,
  cv_summary_file = fold_summary_file
)

metadata_file <- file.path(output_dirs$formatted, "data_formatting_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  04 Data Formatting - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Summary:\n")
cat("  - Pseudo-absence scheme:", pa_scheme_to_use, "\n")
cat("  - Total points:", nrow(all_points), "\n")
cat("  - Removed NA predictor points:", removed_na_points, "\n")
cat("  - Spatial CV folds:", cv_k, "\n")
cat("  - Final spatial block size:", spatial_cv$block_size_km, "km\n")
cat("  - Variables:", length(selected_vars), "\n")
cat("\n")

gc()
