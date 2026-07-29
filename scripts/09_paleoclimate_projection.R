# ============================================================================
# 09_paleoclimate_projection.R - Paleoclimate Projection
# ============================================================================
# Functionality:
#   - Project model onto paleoclimate data
#   - Support per-period runs (controlled by projection_periods in config)
#   - Output flattened directory structure
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  09 Paleoclimate Projection - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Load model and data
# ============================================================================

# Load single-model training results
model_results_file <- file.path(output_dirs$models, "model_results.RData")
if (!file.exists(model_results_file)) {
  stop("Error: Single model results not found. Run 05_model_training.R first.")
}
load(model_results_file)
cat("  Loaded single model results\n")

# Load ensemble model
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

# Load formatted data (for coordinates etc.)
data_file_temp <- file.path(output_dirs$formatted, "biomod2_formatted_data.RData")
if (!file.exists(data_file_temp)) {
  stop("Error: Formatted data not found.")
}
load(data_file_temp)
cat("  Loaded formatted data\n")

# ============================================================================
# 2. Get list of periods to project
# ============================================================================

cat("\n  Scanning available paleoclimate periods...\n")

# Use config function to get period list
all_periods <- list_paleo_periods()
cat("  Found", length(all_periods), "periods in", current_past_clim, "\n")

# Filter by projection_periods in config
if (is.null(projection_periods)) {
  # Default: project all periods
  periods_to_project <- all_periods
  cat("  Projection periods: ALL (", length(periods_to_project), ")\n")
} else if (is.character(projection_periods)) {
  # Specified period list
  periods_to_project <- intersect(projection_periods, all_periods)
  cat("  Projection periods: CUSTOM (", length(periods_to_project), ")\n")
} else if (is.list(projection_periods) && !is.null(projection_periods$time_range)) {
  # Filter by time range (extract time from period names)
  time_range <- projection_periods$time_range
  # Extract numeric value from period name (unit: ka)
  period_years <- as.numeric(sub(".*_([0-9]+)k", "\\1", all_periods))
  periods_to_project <- all_periods[period_years >= time_range[1] & period_years <= time_range[2]]
  cat("  Projection periods: TIME RANGE", time_range[1], "-", time_range[2], "ka (", length(periods_to_project), ")\n")
} else {
  periods_to_project <- all_periods
  cat("  Projection periods: ALL (fallback)\n")
}

cat("  Periods to project:", paste(head(periods_to_project, 5), collapse = ", "))
if (length(periods_to_project) > 5) cat(" ... ")
cat("\n")

if (length(periods_to_project) == 0) {
  stop("Error: No periods to project. Check projection_periods in config.R")
}

# ============================================================================
# 3. Create output directories
# ============================================================================

# Main output directory
output_paleo_base <- "./output/09_paleoclimate"
dir.create(output_paleo_base, showWarnings = FALSE, recursive = TRUE)

# Projection results directory
output_proj_dir <- file.path(output_paleo_base, "projections")
dir.create(output_proj_dir, showWarnings = FALSE, recursive = TRUE)

# Statistics directory
output_stats_dir <- file.path(output_paleo_base, "stats")
dir.create(output_stats_dir, showWarnings = FALSE, recursive = TRUE)

cat("  Output base:", output_paleo_base, "\n")

# ============================================================================
# 4. Get single-model list (for projection)
# ============================================================================

selected_file <- file.path(output_dirs$evaluation, "selected_models.rds")
selected_models <- readRDS(selected_file)
all_built_models <- get_built_models(myBiomodModel)

models_to_project <- as.character(unlist(lapply(selected_models, function(algo) {
  pattern <- paste0(algo, "$")
  grep(pattern, all_built_models, value = TRUE)
}), use.names = FALSE))
names(models_to_project) <- NULL

cat("  Models to project:", length(models_to_project), "\n")

# ============================================================================
# 5. Project per period
# ============================================================================

cat("\n")
cat("=================================================================\n")
cat("  Starting projection for", length(periods_to_project), "periods\n")
cat("=================================================================\n")

# Store statistics
all_stats <- data.frame()

# Set projection name prefix
proj_prefix <- paste0("paleo_", projection_name)

# Progress tracking
n_periods <- length(periods_to_project)

for (i in seq_along(periods_to_project)) {
  period <- periods_to_project[i]

  cat("\n  [", i, "/", n_periods, "] Processing:", period, "\n")

  # ----- 5.1 Load climate data for this period -----
  load_error <- FALSE

  tryCatch({
    paleo_clim <- load_paleo_climate(period)

    # Keep only selected variables
    if (!is.null(selected_vars)) {
      available_vars <- names(paleo_clim)[names(paleo_clim) %in% selected_vars]
      if (length(available_vars) > 0) {
        paleo_clim <- paleo_clim[[available_vars]]
      } else {
        cat("    Warning: Selected variables not found in paleoclimate data\n")
      }
    }

    cat("    Loaded", nlyr(paleo_clim), "climate variables\n")

  }, error = function(e) {
    cat("    Error loading paleoclimate:", e$message, "\n")
    load_error <<- TRUE
  })

  # Skip periods that failed to load
  if (load_error) {
    cat("    Skipping this period due to error\n")
    next
  }

  # ----- 5.2 Single-model projection -----
  proj_name <- paste0(proj_prefix, "_", period)

  # Pre-create the projection output directory; newer terra no longer
  # auto-creates it, so biomod2's writeRaster fails with "path does not exist".
  dir.create(file.path(getwd(), "Models", myBiomodModel@sp.name,
                       paste0("proj_", proj_name), "individual_projections"),
             recursive = TRUE, showWarnings = FALSE)

  cat("    Running single model projection...")

  myBiomodProj <- BIOMOD_Projection(
    bm.mod = myBiomodModel,
    proj.name = proj_name,
    new.env = paleo_clim,
    models.chosen = models_to_project,
    metric.binary = NULL,
    build.clamping.mask = TRUE,
    output.format = output_format,
    nb.cpu = 2,
    keep.in.memory = FALSE,
    do.stack = FALSE,
    overwrite = TRUE
  )

  cat(" done\n")

  # ----- 5.3 Ensemble forecast -----
  cat("    Running ensemble forecasting...")

  myBiomodEF <- BIOMOD_EnsembleForecasting(
    bm.em = ensemble_model,
    bm.proj = myBiomodProj,
    models.chosen = 'all',
    metric.binary = NULL,
    nb.cpu = 2,
    keep.in.memory = FALSE,
    do.stack = FALSE,
    output.format = output_format
  )

  cat(" done\n")

  # ----- 5.4 Extract and save statistics -----
  cat("    Extracting statistics...")

  sp_name <- myBiomodModel@sp.name
  proj_dir <- file.path(getwd(), "Models", sp_name, paste0("proj_", proj_name), "individual_projections")

  if (dir.exists(proj_dir)) {
    proj_files <- list.files(proj_dir, pattern = "\\.tif$", full.names = TRUE)

    # Extract EMwmean result (recommended)
    em_files <- proj_files[grepl("EMwmean", proj_files) & !grepl("_TSSbin", proj_files)]

    if (length(em_files) > 0) {
      # Copy to output directory
      for (f in em_files) {
        new_name <- paste0(period, "_", basename(f))
        file.copy(f, file.path(output_proj_dir, new_name), overwrite = TRUE)
      }

      # Compute statistics
      r <- rast(em_files[1])
      stats <- data.frame(
        period = period,
        min = global(r, "min", na.rm = TRUE)[1, 1],
        max = global(r, "max", na.rm = TRUE)[1, 1],
        mean = global(r, "mean", na.rm = TRUE)[1, 1],
        sd = global(r, "sd", na.rm = TRUE)[1, 1],
        na_count = global(is.na(r), "sum", na.rm = TRUE)[1, 1]
      )

      # Save per-period statistics
      stats_file <- file.path(output_stats_dir, paste0(period, "_summary.csv"))
      write.csv(stats, stats_file, row.names = FALSE)

      # Append to summary
      all_stats <- rbind(all_stats, stats)
    }
  }

  cat(" done\n")

  # ----- 5.5 Memory cleanup -----
  if (i %% gc_trigger == 0) {
    gc()
    cat("    [Garbage collected]\n")
  }

  # Remove current period's model objects to free memory
  rm(myBiomodProj, myBiomodEF, paleo_clim)
}

# ============================================================================
# 6. Save summary statistics
# ============================================================================

cat("\n  Collecting all available statistics files...\n")

# 1. List all summary CSV files in the stats directory
stat_files <- list.files(output_stats_dir, pattern = "_summary.csv$", full.names = TRUE)

if (length(stat_files) > 0) {
  cat("  Found", length(stat_files), "individual period stats. Merging...\n")

  # 2. Read and merge all files
  all_stats <- do.call(rbind, lapply(stat_files, function(f) {
    read.csv(f, stringsAsFactors = FALSE)
  }))

  # 3. Sort by period name (optional)
  # all_stats <- all_stats[order(all_stats$period), ]

  # 4. Save the final merged summary
  summary_file <- file.path(output_paleo_base, "summary_all_periods.csv")
  write.csv(all_stats, summary_file, row.names = FALSE)
  cat("  Saved final summary:", summary_file, "\n")

  cat("\n  Projection summary:\n")
  print(head(all_stats, 10))
} else {
  cat("  Warning: No statistics files found in", output_stats_dir, "\n")
}

# ============================================================================
# 7. Metadata
# ============================================================================

metadata <- list(
  projection_type = "paleoclimate",
  climate_data_source = current_past_clim,
  climate_data_dir = past_clim_dir,
  variables_used = if (!is.null(selected_vars)) selected_vars else "all",
  periods_projected = periods_to_project,
  n_periods = length(periods_to_project),
  projection_name = projection_name,
  output_dir = output_paleo_base,
  created_at = Sys.time()
)

metadata_file <- file.path(output_paleo_base, "projection_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  09 Paleoclimate Projection - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Projected", length(periods_to_project), "periods\n")
cat("  Output:", output_paleo_base, "\n")
cat("\n")

# Cleanup
gc()
