# ============================================================================
# 08_projection.R - Model Projection
# ============================================================================
# Functionality:
#   - Project model onto modern climate data (model validation)
#   - Prepare for subsequent paleoclimate projection
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  08 Projection - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Load data
# ============================================================================

# Load model training results (myBiomodModel)
model_results_file <- file.path(output_dirs$models, "model_results.RData")
if (!file.exists(model_results_file)) {
  stop("Error: Single model results not found. Run 05_model_training.R first.")
}
load(model_results_file)
cat("  Loaded single model results\n")

# Load ensemble model (ensemble_model)
ensemble_file <- file.path(output_dirs$ensemble, "ensemble_model.RData")
if (!file.exists(ensemble_file)) {
  stop("Error: Ensemble model not found. Run 07_ensemble_modeling.R first.")
}
load(ensemble_file)
cat("  Loaded ensemble model:", ensemble_file, "\n")

# Load selected variables
variables_file <- file.path(output_dirs$variable_select, "variables_selected.csv")
if (file.exists(variables_file)) {
  selected_vars <- read.csv(variables_file)$variable
  cat("  Selected variables:", paste(selected_vars, collapse = ", "), "\n")
} else {
  selected_vars <- NULL
  cat("  Variables not found, using all available\n")
}

training_metadata_file <- file.path(output_dirs$models, "model_training_metadata.rds")
training_metadata <- NULL
if (file.exists(training_metadata_file)) {
  training_metadata <- readRDS(training_metadata_file)
}

evaluation_summary_file <- file.path(output_dirs$models, "model_evaluation_summary.csv")
evaluation_summary <- NULL
if (file.exists(evaluation_summary_file)) {
  evaluation_summary <- read.csv(evaluation_summary_file, stringsAsFactors = FALSE)
}

# Load modern climate data (uses standardized function from config.R)
cat("\n  Loading and formatting modern climate data...\n")

modern_clim <- load_standard_climate(clim_dir = modern_clim_dir, pattern = modern_clim_pattern)

cat("  Loaded modern climate variables:", paste(names(modern_clim), collapse = ", "), "\n")

# Keep only selected variables
if (!is.null(selected_vars)) {
  available_vars <- names(modern_clim)[names(modern_clim) %in% selected_vars]
  if (length(available_vars) > 0) {
    modern_clim <- modern_clim[[available_vars]]
    cat("  Using", length(available_vars), "selected variables for projection:\n")
    cat("  ", paste(available_vars, collapse = ", "), "\n")
  } else {
    stop("Error: Selected variables do not match modern climate layers.")
  }
}

data_file_temp <- file.path(output_dirs$formatted, "biomod2_formatted_data.RData") 
if (!file.exists(data_file_temp)) {
  # If not saved separately, try extracting from model_results.RData
  # or from the myBiomodModel object (if linked)
  stop("Error: myBiomodData.RData not found. Needed for coordinates.")
}
load(data_file_temp) # after loading, object name is biomod2_formatted_data.RData
cat("  Loaded formatted data (biomod2_formatted_data.RData)\n")

normalize_prediction_raster <- function(x) {
  xmax <- terra::global(x, "max", na.rm = TRUE)[1, 1]
  if (is.finite(xmax) && xmax > 1) {
    x <- x / 1000
  }
  x
}

get_best_algo <- function(summary_df, metric = "TSS", allowed_algos = NULL) {
  if (is.null(summary_df)) {
    return(NA_character_)
  }
  metric_col <- paste0("mean_", metric)
  if (!metric_col %in% names(summary_df)) {
    return(NA_character_)
  }
  candidate_df <- summary_df
  if (!is.null(allowed_algos)) {
    candidate_df <- candidate_df[candidate_df$algo %in% allowed_algos, , drop = FALSE]
  }
  valid_idx <- which(is.finite(candidate_df[[metric_col]]))
  if (length(valid_idx) == 0) {
    return(NA_character_)
  }
  candidate_df$algo[valid_idx][which.max(candidate_df[[metric_col]][valid_idx])]
}

plot_projection_map <- function(r, title_text, land_vect, presence_xy = NULL,
                                palette, show_points = TRUE) {
  r <- normalize_prediction_raster(r)
  plot(
    r,
    main = title_text,
    col = palette,
    zlim = c(0, 1),
    range = c(0, 1),
    plg = list(title = "Habitat suitability")
  )
  plot(
    land_vect,
    add = TRUE,
    border = "grey30",
    lwd = 0.4,
    col = NA
  )
  if (show_points && !is.null(presence_xy) && nrow(presence_xy) > 0) {
    points(presence_xy$x, presence_xy$y, col = adjustcolor(presence_color, alpha.f = 0.75),
           pch = presence_pch, cex = 0.35)
  }
}

extract_presence_background_values <- function(r, presence_xy, background_n = 5000,
                                               seed_value = 123) {
  r <- normalize_prediction_raster(r)
  presence_vals <- numeric(0)
  if (!is.null(presence_xy) && nrow(presence_xy) > 0) {
    presence_vect <- terra::vect(presence_xy, geom = c("x", "y"), crs = terra::crs(r))
    extracted_presence <- terra::extract(r, presence_vect)[, 2]
    presence_vals <- extracted_presence[is.finite(extracted_presence)]
  }
  
  all_vals <- terra::values(r, mat = FALSE)
  all_vals <- all_vals[is.finite(all_vals)]
  if (length(all_vals) == 0) {
    return(list(presence = numeric(0), background = numeric(0)))
  }
  
  set.seed(seed_value)
  background_size <- min(background_n, length(all_vals))
  background_vals <- sample(all_vals, size = background_size, replace = FALSE)
  
  list(presence = presence_vals, background = background_vals)
}

plot_presence_background_comparison <- function(r, title_text, presence_xy,
                                                background_n = 5000,
                                                seed_value = 123) {
  comparison_vals <- extract_presence_background_values(
    r = r,
    presence_xy = presence_xy,
    background_n = background_n,
    seed_value = seed_value
  )
  
  if (length(comparison_vals$presence) == 0 || length(comparison_vals$background) == 0) {
    plot.new()
    text(0.5, 0.5, paste("No valid comparison values for", title_text))
    return(invisible(NULL))
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  hist(
    comparison_vals$background,
    breaks = seq(0, 1, by = 0.02),
    freq = FALSE,
    main = paste("Presence vs Background:", title_text),
    xlab = "Suitability (0-1)",
    col = rgb(0.6, 0.6, 0.6, 0.6),
    border = "white",
    xlim = c(0, 1)
  )
  hist(
    comparison_vals$presence,
    breaks = seq(0, 1, by = 0.02),
    freq = FALSE,
    add = TRUE,
    col = rgb(0.85, 0.2, 0.2, 0.45),
    border = "white"
  )
  legend(
    "topright",
    legend = c(
      paste0("Background (n=", length(comparison_vals$background), ")"),
      paste0("Presence (n=", length(comparison_vals$presence), ")")
    ),
    fill = c(rgb(0.6, 0.6, 0.6, 0.6), rgb(0.85, 0.2, 0.2, 0.45)),
    bg = "white",
    cex = 0.85
  )
  
  boxplot(
    list(
      Background = comparison_vals$background,
      Presence = comparison_vals$presence
    ),
    col = c("grey75", "tomato"),
    main = "Distribution summary",
    ylab = "Suitability (0-1)",
    ylim = c(0, 1)
  )
}
# ============================================================================
# 2. Project to modern climate (using BIOMOD_Projection properly)
# ============================================================================

cat("\n  Projecting to modern climate...\n")

# 2.1 Prepare single-model list (exact match required)
selected_file <- file.path(output_dirs$evaluation, "selected_models.rds")
selected_models <- readRDS(selected_file)
all_built_models <- get_built_models(myBiomodModel)

models_to_project <- as.character(unlist(lapply(selected_models, function(algo) {
  pattern <- paste0(algo, "$")
  grep(pattern, all_built_models, value = TRUE)
}), use.names = FALSE))
names(models_to_project) <- NULL

best_tss_algo <- if (!is.null(training_metadata) && !is.null(training_metadata$best_model_tss)) {
  if (training_metadata$best_model_tss %in% selected_models) training_metadata$best_model_tss else
    get_best_algo(evaluation_summary, "TSS", allowed_algos = selected_models)
} else {
  get_best_algo(evaluation_summary, "TSS", allowed_algos = selected_models)
}
best_boyce_algo <- get_best_algo(evaluation_summary, "BOYCE", allowed_algos = selected_models)

cat("  Best single model by TSS:", best_tss_algo, "\n")
cat("  Best single model by BOYCE:", best_boyce_algo, "\n")

# 2.2 Run single-model projection (BIOMOD_Projection)
proj_name <- paste0("modern_", projection_name)

# Pre-create the projection output directory; newer terra no longer auto-creates
# it, so biomod2's writeRaster fails with "[writeRaster] path does not exist".
dir.create(file.path(getwd(), "Models", myBiomodModel@sp.name,
                     paste0("proj_", proj_name), "individual_projections"),
           recursive = TRUE, showWarnings = FALSE)

cat("  Running Single Models Projection (", proj_name, ")...\n")
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModel,
  proj.name = proj_name,
  new.env = modern_clim,
  models.chosen = models_to_project,
  metric.binary = NULL, # set to NULL to skip binarization
  build.clamping.mask = FALSE,
  output.format = output_format,
  nb.cpu = 1,
  keep.in.memory = FALSE,     # FALSE = drop from memory after each map, keep only .tif on disk
  do.stack = FALSE,           # FALSE = don't stack dozens of maps into one huge object
  overwrite = TRUE
)
cat("  Single model projection completed\n")

# 2.3 Ensemble forecasting (BIOMOD_EnsembleForecasting)
# Build ensemble from single-model projections (bm.proj) - most stable approach
cat("\n  Running Ensemble Models Forecasting...\n")
myBiomodEF <- BIOMOD_EnsembleForecasting(
  bm.em = ensemble_model,
  bm.proj = myBiomodProj,
  models.chosen = 'all',
  metric.binary = NULL, #metric_binary, # <--- keep "TSS" here to binarize only a few ensemble models
  nb.cpu = 1,
  keep.in.memory = FALSE,
  do.stack = FALSE,
  output.format = ".tif"
)
cat("  Ensemble forecasting completed\n")

# ============================================================================
# 3. Extract projection statistics (memory-safe approach)
# ============================================================================

cat("\n  Extracting projection results...\n")

# Note: biomod2 saves files under "species_name/proj_proj_name" by default
sp_name <- myBiomodModel@sp.name
proj_dir <- file.path(getwd(), "Models", sp_name, paste0("proj_", proj_name), "individual_projections")

if (dir.exists(proj_dir)) {
  proj_files <- list.files(proj_dir,  pattern = "\\.tif$", full.names = TRUE)
  
  if (length(proj_files) > 0) {
    cat("  Found", length(proj_files), "projection files\n")
    
    # Use terra::global for memory-safe statistics
    proj_stats <- data.frame(
      file = basename(proj_files),
      min = sapply(proj_files, function(f) terra::global(rast(f), "min", na.rm = TRUE)[1, 1]),
      max = sapply(proj_files, function(f) terra::global(rast(f), "max", na.rm = TRUE)[1, 1]),
      mean = sapply(proj_files, function(f) terra::global(rast(f), "mean", na.rm = TRUE)[1, 1]),
      na_count = sapply(proj_files, function(f) terra::global(is.na(rast(f)), "sum", na.rm = TRUE)[1, 1])
    )
    
    cat("\n  Projection statistics:\n")
    print(proj_stats)
  }
} else {
  cat("  Warning: Projection output directory not found at", proj_dir, "\n")
}

# ============================================================================
# 4. Save projection results
# ============================================================================

cat("\n  Saving projection objects...\n")
dir.create(output_dirs$projection, showWarnings = FALSE, recursive = TRUE)

proj_single_file <- file.path(output_dirs$projection, paste0("myBiomodProj_", projection_name, ".RData"))
save(myBiomodProj, file = proj_single_file)

proj_ens_file <- file.path(output_dirs$projection, paste0("myBiomodEF_", projection_name, ".RData"))
save(myBiomodEF, file = proj_ens_file)

if (exists("proj_stats")) {
  stats_file <- file.path(output_dirs$projection, paste0("projection_", projection_name, "_summary.csv"))
  write.csv(proj_stats, stats_file, row.names = FALSE)
}
cat("  Results saved to", output_dirs$projection, "\n")

# ============================================================================
# 5. Visualization - modern climate projection map
# ============================================================================

cat("\n  Generating visualization...\n")
# ================== Generate land boundary (based on climate data NA) ==================
cat("  Creating land mask from modern climate...\n")

# Build mask from first climate variable (all variables share the same NA distribution)
land_mask <- !is.na(modern_clim[[1]])

# Convert to vector boundary (dissolve=TRUE merges into single land polygon)
land_vect <- terra::as.polygons(land_mask, dissolve = TRUE)

cat("  Land mask created\n")

pdf_file <- file.path(output_dirs$projection, "projection_modern_map.pdf")
pdf(pdf_file, width = 14, height = 10)

# Safely extract presence points from model object for overlay
myResp <- get_formal_data(myBiomodModel, "resp.var")
myXY <- formatted_data@coord
# Keep only true presence points (value = 1)
pres_xy <- myXY[myResp == 1 & !is.na(myResp), ]

if (dir.exists(proj_dir)) {
  proj_ensemble_hist <- NULL
  proj_single_tss_hist <- NULL
  proj_single_boyce_hist <- NULL
  
  # --- 1. Plot ensemble (EMwmean) continuous probability map ---
  # Exclude binary maps, keep continuous probability maps
  ensemble_files <- list.files(proj_dir, pattern = "EMwmean", full.names = TRUE)
  prob_ens_files <- ensemble_files[!grepl("_TSSbin", ensemble_files)]
  
  if (length(prob_ens_files) > 0) {
    proj_ensemble <- normalize_prediction_raster(rast(prob_ens_files[1]))
    proj_ensemble_hist <- proj_ensemble
    plot_projection_map(
      proj_ensemble,
      "Modern suitability: ensemble (EMwmean)",
      land_vect,
      pres_xy,
      suitability_palette,
      show_points = TRUE
    )
    legend("bottomright", legend = "Presence", col = presence_color, pch = presence_pch, bg = "white")
  }
  
  # --- 2. Plot best TSS single model ---
  if (!is.na(best_tss_algo)) {
    single_files_tss <- list.files(proj_dir, pattern = best_tss_algo, full.names = TRUE)
    prob_single_tss <- single_files_tss[!grepl("_TSSbin", single_files_tss)]
    
    if (length(prob_single_tss) > 0) {
      proj_single_tss <- normalize_prediction_raster(rast(prob_single_tss[1]))
      proj_single_tss_hist <- proj_single_tss
      plot_projection_map(
        proj_single_tss,
        sprintf("Modern suitability: best TSS model (%s)", best_tss_algo),
        land_vect,
        pres_xy,
        suitability_palette,
        show_points = TRUE
      )
    }
  }
  
  # --- 3. Plot best BOYCE single model ---
  if (!is.na(best_boyce_algo)) {
    single_files_boyce <- list.files(proj_dir, pattern = best_boyce_algo, full.names = TRUE)
    prob_single_boyce <- single_files_boyce[!grepl("_TSSbin", single_files_boyce)]
    
    if (length(prob_single_boyce) > 0) {
      proj_single_boyce <- normalize_prediction_raster(rast(prob_single_boyce[1]))
      proj_single_boyce_hist <- proj_single_boyce
      plot_projection_map(
        proj_single_boyce,
        sprintf("Modern suitability: best BOYCE model (%s)", best_boyce_algo),
        land_vect,
        pres_xy,
        suitability_palette,
        show_points = TRUE
      )
    }
  }
  
  # --- 4. Side-by-side comparison: best TSS vs best BOYCE ---
  if (!is.na(best_tss_algo) && !is.na(best_boyce_algo)) {
    tss_files <- list.files(proj_dir, pattern = best_tss_algo, full.names = TRUE)
    boyce_files <- list.files(proj_dir, pattern = best_boyce_algo, full.names = TRUE)
    tss_prob <- tss_files[!grepl("_TSSbin", tss_files)]
    boyce_prob <- boyce_files[!grepl("_TSSbin", boyce_files)]
    
    if (length(tss_prob) > 0 && length(boyce_prob) > 0) {
      par(mfrow = c(1, 2), mar = c(3, 3, 3, 5))
      plot_projection_map(
        normalize_prediction_raster(rast(tss_prob[1])),
        sprintf("Best TSS: %s", best_tss_algo),
        land_vect,
        pres_xy,
        suitability_palette,
        show_points = FALSE
      )
      plot_projection_map(
        normalize_prediction_raster(rast(boyce_prob[1])),
        sprintf("Best BOYCE: %s", best_boyce_algo),
        land_vect,
        pres_xy,
        suitability_palette,
        show_points = FALSE
      )
      par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    }
  }
  
  # --- 5. Projection value histograms ---
  hist_targets <- list(
    proj_ensemble_hist,
    proj_single_tss_hist,
    proj_single_boyce_hist
  )
  names(hist_targets) <- c(
    "Ensemble (EMwmean)",
    paste0("Best TSS (", best_tss_algo, ")"),
    paste0("Best BOYCE (", best_boyce_algo, ")")
  )
  hist_targets <- hist_targets[!vapply(hist_targets, is.null, logical(1))]
  
  if (length(hist_targets) > 0) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(length(hist_targets), 1), mar = c(4, 4, 3, 1))
    
    for (hist_name in names(hist_targets)) {
      hist_values <- terra::values(hist_targets[[hist_name]], mat = FALSE)
      hist_values <- hist_values[is.finite(hist_values)]
      hist(
        hist_values,
        breaks = seq(0, 1, by = 0.02),
        main = paste("Suitability histogram:", hist_name),
        xlab = "Suitability (0-1)",
        col = "grey70",
        border = "white",
        xlim = c(0, 1)
      )
    }
    
    par(old_par)
  }
  
  # --- 6. Presence vs background comparison ---
  if (length(hist_targets) > 0) {
    for (hist_name in names(hist_targets)) {
      plot_presence_background_comparison(
        r = hist_targets[[hist_name]],
        title_text = hist_name,
        presence_xy = pres_xy,
        background_n = 5000,
        seed_value = seed_val
      )
    }
  }
  
} else {
  plot.new()
  text(0.5, 0.5, "Projection maps not found.")
}

dev.off()
cat("  Saved visualization:", pdf_file, "\n")

# ============================================================================
# 6. Save metadata
# ============================================================================

metadata <- list(
  projection_type = "modern",
  climate_data = modern_clim_dir,
  variables_used = if (!is.null(selected_vars)) selected_vars else "all",
  output_dir = proj_dir,
  created_at = Sys.time()
)

metadata_file <- file.path(output_dirs$projection, "projection_modern_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  08 Projection - Completed\n")
cat("=================================================================\n")
cat("\n")

# Clean up memory
gc()
