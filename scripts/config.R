# ============================================================================
# config.R
# Global configuration for the Sargentodoxa SDM project
#
# Species: Sargentodoxa cuneata
# Main pipeline:
#   1. Occurrence cleaning and spatial thinning
#   2. Pseudo-absence generation
#   3. Environmental variable selection
#   4. Spatial cross-validation
#   5. Single-model training, filtering, and ensemble
#   6. Modern / paleoclimate projection and visualization
#
# Tips:
#   - Prefer changing parameter values; avoid editing helper functions below.
#   - Sections follow the analysis order; review top to bottom.
#   - When switching climate data source or projection version, check path config first.
#   - Each section tags parameters as:
#       Common: likely to change on first run
#       Advanced: tune only with a clear modeling goal
#       Leave default: keep unless debugging or refactoring
# ============================================================================

# ==================== 1. Path configuration ====================
# This repo is de-servered:
#   - Local inputs (occurrence / fossil / TEYF metadata) ship with the repo via relative paths.
#   - Server-side climate rasters (WorldClim / TEYF) are configured via env vars; defaults
#     are empty. Local copies are not bundled with the repo.
#   - Plotting scripts in scripts/visualization/ infer project_root themselves and do not
#     depend on this section; it mainly serves pipeline scripts 01-09.
# Convention: run pipeline scripts from the repo root, e.g. Rscript scripts/01_data_preparation.R

# Infer repo root (works with both source() and Rscript, independent of cwd).
# config.R lives at scripts/config.R, so the root is one level up.
project_dir <- tryCatch({
  ofiles <- vapply(sys.frames(), function(f) f$ofile, character(1))
  ofile <- ofiles[length(ofiles)]
  if (!is.na(ofile) && nzchar(ofile)) {
    normalizePath(file.path(dirname(ofile), ".."), mustWork = FALSE)
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    fa <- args[grep("^--file=", args)]
    if (length(fa)) {
      normalizePath(file.path(dirname(sub("^--file=", "", fa[1])), ".."), mustWork = FALSE)
    } else {
      normalizePath(".", mustWork = FALSE)
    }
  }
}, error = function(e) normalizePath(".", mustWork = FALSE))

data_dir <- file.path(project_dir, "data")

# Common: occurrence data (CVH / GBIF)
# To switch data source, usually just change sp_file below.
sp_file_cvh <- file.path(project_dir, "data/occ_data_cvh/sargentodoxa_occ_cvh_counties_filtered.csv")
sp_file_gbif <- file.path(project_dir, "data/occ_data_cvh/sargentodoxa_occ_filtered.csv")

Sys.setenv(
  SARGENTO_MODERN_CLIM_DIR = "/home/cedrus/data/shz/data/wc2.1_5m_bio",
  SARGENTO_TEWT_CLIM_DIR   = "/home/cedrus/data/shz/data/tewt_0_3600ka/masked_landsea",
  SARGENTO_TEYF_CLIM_DIR   = "/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min"
)
# Leave default: modern climate directory (WorldClim 2.1, 5 arc-min)
# Set the SARGENTO_MODERN_CLIM_DIR env var to your climate data directory.
# Not bundled; set the env var to rerun steps 01-08 locally.
modern_clim_dir <- Sys.getenv(
  "SARGENTO_MODERN_CLIM_DIR",
  unset = ""
)
modern_clim_pattern <- "*.tif$"

# Common: paleoclimate data source
# To switch between the two paleoclimate datasets, usually just change current_past_clim.
# 1. TEWT 0-3600 ka (set the SARGENTO_TEWT_CLIM_DIR env var to your climate data directory)
past_clim_dir_tewt <- Sys.getenv(
  "SARGENTO_TEWT_CLIM_DIR",
  unset = ""
)
past_clim_subdir_tewt <- TRUE  # TRUE = one subdir per period

# 2. TEYF (set the SARGENTO_TEYF_CLIM_DIR env var to your climate data directory)
past_clim_dir_teyf <- Sys.getenv(
  "SARGENTO_TEYF_CLIM_DIR",
  unset = ""
)
past_clim_subdir_teyf <- TRUE  # TRUE = one subdir per period

# Common: currently active paleoclimate source
# Options: "tewt" / "teyf"
current_past_clim <- "teyf"

# Auto-select the active paleoclimate directory based on current_past_clim
if (current_past_clim == "tewt") {
  past_clim_dir <- past_clim_dir_tewt
  past_clim_subdir <- past_clim_subdir_tewt
} else {
  past_clim_dir <- past_clim_dir_teyf
  past_clim_subdir <- past_clim_subdir_teyf
}

# Advanced: land boundary (optional, visualization only)
# Set to NULL if no shapefile is available
land_file <- NULL

# Common: fossil records
# Used to compare paleoclimate projections against fossil evidence.
fossil_dir <- file.path(project_dir, "data/fossil_record")
# fossil_record_file must contain at least: paleoLong, paleoLat, age, simulation_name
fossil_record_file <- file.path(fossil_dir, "FossilRecord-paleocor-v2.1.csv")

# Advanced: paleoclimate period ordering file
# Bundled with the repo (data/teyf_simulations.csv).
teyf_simulations_file <- file.path(project_dir, "data/teyf_simulations.csv")
simulations_file <- teyf_simulations_file

# ==================== 2. Species and input data ====================

species_name <- "Sargentodoxa"
resp_name <- "Sargentodoxa"

# Common: default occurrence data source
# On a first run this is usually the first input to confirm.
sp_file <- sp_file_cvh

# ==================== 3. Compatibility aliases ====================
# This section mainly supports older scripts; new scripts that use modern_clim_dir /
# past_clim_dir uniformly usually need no changes here.

# Leave default: modern climate directory alias
wc2_clim_dir <- modern_clim_dir

# Leave default: paleoclimate directory alias
tewt_dir <- past_clim_dir

# ==================== 4. Output directories ====================

output_base <- "./output"
dir.create(output_base, showWarnings = FALSE)

# Leave default: per-step output directories
# Usually no change needed unless restructuring project outputs.
output_dirs <- list(
  data_prep = "./output/01_data_prep",
  pseudo_absence = "./output/02_pseudo_absence",
  variable_select = "./output/03_variable_select",
  formatted = "./output/04_formatted",
  models = "./output/05_models",
  evaluation = "./output/06_evaluation",
  ensemble = "./output/07_ensemble",
  projection = "./output/08_projection",
  paleoclimate = "./output/09_paleoclimate",
  visualization = "./output/10_visualization"

)

# Auto-create output directories on startup
for (dir_name in output_dirs) {
  dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
}

# ==================== 5. Spatial thinning ====================

# Common: spatially thin presence points in `01_data_preparation`?
# Keep TRUE for data with obvious sampling aggregation.
enable_spatial_thinning <- TRUE

# Common: thinning method
# `hybrid` is usually the safest first choice.
# "cell"     = keep 1 point per climate raster cell
# "distance" = thin by minimum distance only
# "hybrid"   = cell dedup, then distance constraint
thinning_method <- "hybrid"

# Advanced: use modern climate raster as the grid for cell thinning?
# Keep TRUE to align occurrences with predictor rasters.
thinning_use_predictor_cell <- TRUE

# Common: distance threshold source
# "auto"   = derive from modern climate raster resolution
# "manual" = use thinning_distance_km directly
thinning_distance_mode <- "manual"

# Advanced: auto threshold formula
# modern raster resolution (km) * thinning_distance_factor
thinning_distance_factor <- 1.0

# Common: manual minimum retained distance (km)
# Core parameter controlling thinning intensity.
thinning_distance_km <- 20

# Advanced: when points conflict, which order decides "who stays"
# "random" = greedy random order (recommended; reduces input-order bias)
# "first"  = keep by original input order
thinning_keep <- "random"

# Leave default: thinning random seed (for reproducibility)
thinning_seed <- 42

# Advanced: nearest-neighbor distance thresholds reported in diagnostics (km)
# Affects diagnostic tables/plots only, not the thinning result.
thinning_diagnostic_thresholds_km <- c(10, 20, 30)

# Advanced: x-axis upper limit for nearest-neighbor distance histogram (km)
thinning_histogram_xlim_km <- 100

# ==================== 6. Pseudo-absence ====================

# Common: pseudo-absence generation strategy
# `kde` is fine for a first run.
# "kde"    = weighted sampling by occurrence density
# "disk"   = sample within a distance band around presence
# "random" = random sampling
# "sre"    = Surface Range Envelope
pa_strategy <- "kde"

# Advanced: KDE base params (only when pa_strategy = "kde")
# Controls kernel density estimation details; usually unchanged on first use.
kde_sigma <- "Scott"  # bandwidth: "Scott", "Stoyan", "Diggle", "likelihood", "geometry"
kde_window <- "hull"  # window type: "hull" (convex hull), "rect" (rectangle)

# Common: custom KDE implementation params (for older R versions)
# Bandwidth (km): smoothing scale spreading occurrence hotspots
# Larger = smoother, more dispersed pseudo-absence; smaller = closer to occurrence clusters.
kde_bandwidth_km <- 1200

# Advanced: KDE weight-map computation resolution (km)
# Larger = faster but coarser; NULL uses the native climate raster resolution
kde_compute_resolution_km <- 30

# Advanced: add random jitter within KDE coarse cells to reduce blockiness
kde_jitter_within_cell <- T

# Advanced: candidate oversampling factor
# Larger = easier to satisfy min-distance constraints, but slower
kde_candidate_multiplier <- 10

# Advanced: maximum resampling iterations
kde_max_iter <- 14

# Advanced: SRE params (only when pa_strategy = "sre")
sre_quant <- 0.1      # half-quantile

# Common: distance limits between pseudo-absence and presence (degrees)
# These directly control how near/far pseudo-absences sit from known occurrences.
# Mainly affects `disk`; in the current implementation `kde` also uses them as constraints.
pa_dist_min <- 2  # 1 degree ~ 111 km
pa_dist_max <- 60 # max distance (degrees) from presence for KDE/disk pseudo-absence

# Advanced: pseudo-absence sampling latitude limits
# Exclude high latitudes to avoid unreasonable polar samples
pa_lat_limit_north <- 75   # drop north of 75N (e.g. near the North Pole)
pa_lat_limit_south <- -60  # drop south of 60S (e.g. Antarctica)

# Common: pseudo-absence scheme 1: 1:1 ratio
# Pseudo-absence count equals presence count; good for a baseline test.
# Note: the count is not hardcoded; `nb.absences = "equal"` auto-matches presence count.
pa_scheme1 <- list(
  name = "scheme_1to1",
  nb.absences = "equal",
  nb.rep = 1,
  desc = "1:1 ratio (equal to presence)"
)

# Common: pseudo-absence scheme 2: ~1:3 ratio
# Try this when you want stronger background constraint and discrimination.
# Note: the count is a manual example, not a universal recommendation;
# adjust manually based on presence count and study design.
pa_scheme2 <- list(
  name = "scheme_1000",
  nb.absences = 1000,
  nb.rep = 1,
  desc = "1:10 ratio (1000 points)"
)

# ==================== 7. Data formatting ====================

# Common: data scheme for downstream modeling
# The main pipeline keeps two schemes for clarity and maintenance.
# "scheme_1to1" = equal to presence
# "scheme_1000" = fixed-count pseudo-absence points
pa_scheme_to_use <- "scheme_1to1"

# ==================== 8. Variable selection ====================

# Common: correlation filter threshold (initial screening)
# Lower = stricter; higher = keeps more variables.
corr_threshold <- 0.8

# Common: VIF threshold (final screening)
# Controls multicollinearity; a core variable-selection parameter.
vif_threshold <- 12

# Common: run stepwise VIF removal?
use_vif_step <- TRUE

# Common: manually preset candidate variables
# NULL = rely entirely on auto "correlation + VIF" selection
# BioClim variable meanings:
# If you have clear ecological judgement, restrict candidates here first.
#   Bio1 = Annual mean temp, Bio2 = Mean diurnal range, Bio3 = Isothermality, Bio4 = Temperature seasonality
#   Bio5 = Max temp of warmest month, Bio6 = Min temp of coldest month, Bio7 = Temp annual range
#   Bio8 = Mean temp of wettest quarter, Bio9 = Mean temp of driest quarter, Bio10 = Mean temp of warmest quarter, Bio11 = Mean temp of coldest quarter
#   Bio12 = Annual precipitation, Bio13 = Precip of wettest month, Bio14 = Precip of driest month, Bio15 = Precip seasonality
#   Bio16 = Precip of wettest quarter, Bio17 = Precip of driest quarter, Bio18 = Precip of warmest quarter, Bio19 = Precip of coldest quarter
manual_var_select <- c(#'Bio1','Bio2'
                       'Bio3', #'Bio5','Bio6','Bio7','Bio4',
                       'Bio8','Bio9', 'Bio10', 'Bio11',
                       #'Bio12',#'Bio13','Bio14',
                       'Bio15', 'Bio18','Bio19'#'Bio16','Bio17',
                       )
# ==================== 9. Data split and cross-validation ====================

# Common: enable cross-validation?
# Keep TRUE for formal modeling.
use_cv <- TRUE

# Common: cross-validation strategy
# This pipeline targets spatial CV, so `user.defined` is recommended.
# "kfold"        = standard k-fold CV
# "strat"        = stratified sampling
# "random"       = random split
# "user.defined" = custom spatial folds from `04_data_formatting.R`
cv_strategy <- "user.defined"

cv_nb_rep <- 3        # Common: repeats; under user.defined = number of shifted grids
cv_k <- 5             # Common: folds per grid
cv_perc <- 0.8        # Leave default: used only by non-user.defined strategies

# Common: spatial CV method
# "checkerboard" = checkerboard assignment on a spatial grid
spatial_cv_method <- "checkerboard"

# Common: spatial block size (km)
# Should clearly exceed the thinning distance, else train/test separation is weak.
# A key parameter to understand and tune in spatial CV.
spatial_cv_block_size_km <- 500

# Advanced: projection CRS for spatial blocking (units must be meters)
spatial_cv_projection_crs <- "EPSG:3857"

# Advanced: if some folds lack presence/pseudo-absence, auto-grow block size and retry
spatial_cv_auto_retry <- TRUE
spatial_cv_retry_factor <- 1.5
spatial_cv_max_tries <- 4

# ==================== 10. Model algorithms ====================

# Common: single-model algorithm list
# On first contact, prioritize a stable pipeline over many algorithms.
# Options: "ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF", "RFd"
model_list <- c( "CTA", "FDA", "GAM", "GLM", "MARS", "MAXNET", "RF")
# model_list <- c("MAXNET")
# If some algorithms are unavailable in your environment, trim the list temporarily
# model_list <- c("GAM", "GBM", "GLM", "MAXNET", "RF")

# Common: quick-test config (2 models only)
# Use to check the pipeline runs end to end.
model_list_test <- c("MAXNET", "RF")

# ==================== 11. Model training ====================

# Common: parameter optimization strategy
# "bigboss" = biomod2 recommended params
# "user"    = user-defined params
# "none"    = package defaults
opt_strategy <- "bigboss"

# Common: model evaluation metrics
# Affect single-model evaluation, downstream filtering, and ensembling.
metric_eval <- c("TSS", "AUCroc", "AUCprg", "BOYCE")

# Advanced: variable-importance permutation repeats
# Larger = more stable but slower.
var_importance_rep <- 5

# Common: CPU cores for parallel runs
nb_cpu <- 8

# Leave default: global random seed (for reproducibility)
seed_val <- 1001

# ==================== 12. Ensemble models ====================

# Common: ensemble algorithm
# On a first run, EMwmean alone is usually enough.
# "EMmedian" = median ensemble
# "EMmean"   = mean ensemble
# "EMwmean"  = weighted mean ensemble (currently recommended)
# "EMca"     = committee averaging
# "EMci"     = confidence interval
# "EMcv"     = cross-validation ensemble
# For the simplest, stable result, keep EMwmean only.
ensemble_algo <- c("EMwmean")

# Common: model filtering metrics
# Only models meeting all these thresholds enter the ensemble.
# Options: "TSS", "AUCroc", "AUCprg", "BOYCE"
metric_select <- c("TSS")#, "AUCroc","AUCprg", "BOYCE"

# Common: model filtering thresholds
# Order must correspond one-to-one with metric_select

metric_select_thresh <- c(0.8)#, 0.8, 0.8, 0.6

# Advanced: ensemble grouping
# "all"     = merge all models into one ensemble (recommended)
# "algo"    = ensemble per algorithm
# "all_run" = ensemble per repeat run
em.by <- "all"

# ==================== 13. Single-model filtering ====================

# Common: model selection mode
# If unsure whether thresholds are reasonable, start with `manual` for clearer results.
# "auto"   = filter by thresholds
# "manual" = manually specify kept models
model_selection_mode <- "auto"

# Common: model names to keep under manual mode (one or more)
# e.g.: manual_selected_models <- c("MAXNET", "RF", "GLM")
manual_selected_models <- c("CTA", "GAM", "GLM", "MARS", "MAXNET", "RF")

# Advanced: max top-N models to keep under auto mode
top_n_models <- 5

# Advanced: EMwmean weight decay
# "proportional" = weight proportional to metric
# "exponential"  = exponentially boost high-score models
emwmean_decay <- "proportional"

# ==================== 14. Projection settings ====================

# Common: projection version label
# Distinguishes runs; directly affects output dir and file naming
# e.g. "v1_modern", "v2_teyf", "v3_tewt_petm"
projection_name <- "v5_for_paper_code"

# Common: paleoclimate periods to project
# NULL = scan the data dir and project all periods.
# Or specify manually:
#   1. character vector, e.g. c("tewta_0k", "tewtb_1k")
#   2. time range, e.g. list(time_range = c(0, 100))
#      = project only periods within 0-100 ka (inclusive)
# 09_paleoclimate_projection.R extracts numbers from period names for range filtering.
# For TEWT data, to test within 100 ka:
#   projection_periods <- list(time_range = c(0, 100))
projection_periods <- NULL # dynamic; NULL = auto-scan all

# Leave default: projection raster output format
# ".tif" = GeoTIFF (recommended)
# ".grd" = Raster native format
output_format <- ".tif"

# Common: threshold rule for converting continuous suitability to binary maps
# "TSS" = max TSS
# "ROC" = max ROC
# "all" = output all available rules
metric_binary <- "TSS"

# ==================== 15. Visualization ====================

# Common: map display extent
# Affects plotting only, not modeling.
map_extent <- c(-180, 180, -45, 35)  # lon min, max; lat min, max

# Advanced: suitability color ramp
suitability_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(100)

# Advanced: presence point style
presence_color <- "#336600"
presence_pch <- 19

# Advanced: pseudo-absence point style
pseudoabsence_color <- "#FF6B6B"

# Common: which ensemble metric's projection to prefer for visualization
# Options: "TSS", "AUCroc", "AUCprg", "BOYCE"
# 10_visualization.R prefers the `EMwmeanBy<metric>` tif per period;
# if that metric is missing, it falls back to other available results.
visualization_projection_metric <- "TSS"

# Advanced: title font size
title_size <- 14

# ==================== 16. Geographic parameters ====================

# Leave default: default geographic CRS
target_crs <- "EPSG:4326"  # WGS84

# Advanced: default buffer size (degrees)
buffer_degree <- 0.5

# ==================== 17. Other settings ====================

# Common: save intermediate results?
# Keep TRUE on a first run for easier debugging and traceback.
save_intermediate <- TRUE

# Common: verbose logging?
verbose <- TRUE

# Advanced: memory management
# Trigger gc() after every N periods; mainly stabilizes long batch projections.
gc_trigger <- 20  # call gc() after every N periods

# ============================================================================
# Dependencies
# ============================================================================

required_packages <- c(
  # Core modeling / spatial data
  "biomod2",
  "terra",
  "raster",
  "sf",

  # Data wrangling
  "dplyr",
  "stringr",
  "tidyr",
  "readr",
  "gtools",        # mixedsort

  # Visualization
  "ggplot2",
  "plotly",
  "tidyterra",
  "viridis",      # color maps
  "corrplot",
  "magick",

  # Statistics and auxiliary modeling
  "usdm",          # VIF
  # "spatialEco",    # pseudo.absence version mismatch
  "gamlss",

  # Parallel
  "doParallel",
  "foreach",

  # Other tools
  "reshape",
  "ggtext",
  "htmlwidgets",
  "R.utils",
  "lme4",
  "cito", # version mismatch; install manually
  "xgboost" # version mismatch; install manually
)

# Leave default: check and load dependencies
# Note: missing packages are auto-installed; on restricted servers, install manually beforehand.
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# # Advanced: example manual install for older R
# Matrixurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-1.tar.gz"
#
# # Install a specific version from source
# install.packages(Matrixurl, repos = NULL, type = "source")

# # Advanced: xgboost 1.7.7.1 archive download URL
# url_xgb <- "https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.7.1.tar.gz"
#
# # Install a specific version from source
# install.packages(url_xgb, repos = NULL, type = "source")

# ============================================================================
# Global helper functions
# ============================================================================

#' Load and standardize climate data
#' Reads raster files from a directory, sorts by natural order, renames to Bio1, Bio2...
#' @param clim_dir climate data directory
#' @param pattern filename match pattern, default "*.tif$"
#' @return terra::SpatRaster
load_standard_climate <- function(clim_dir, pattern = "*.tif$") {
  files <- list.files(clim_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop(paste("Error: No climate files found in", clim_dir))
  }

  # Read rasters
  clim <- terra::rast(files)

  # Reorder layers by natural order to avoid lexicographic 1, 10, 11, 2 issues
  clim <- clim[[gtools::mixedsort(names(clim))]]

  # Unify layer names
  names(clim) <- paste0("Bio", 1:terra::nlyr(clim))

  return(clim)
}

# ============================================================================
# Spatial thinning helper functions
# ============================================================================

#' Estimate raster resolution (km)
#' Converts lon/lat resolution of the first climate layer to km near the study-area center
#' @param clim SpatRaster
#' @param ref_points optional SpatVector points for estimating center latitude
#' @return list(x_km, y_km, mean_km)
estimate_raster_resolution_km <- function(clim, ref_points = NULL) {
  raster_res <- terra::res(clim)

  if (!is.null(ref_points) && inherits(ref_points, "SpatVector")) {
    center_lat <- mean(terra::crds(ref_points)[, 2], na.rm = TRUE)
  } else {
    center_lat <- mean(terra::ext(clim)[3:4])
  }

  lat_km <- 111.32 * abs(raster_res[2])
  lon_km <- 111.32 * cos(center_lat * pi / 180) * abs(raster_res[1])

  list(
    x_km = lon_km,
    y_km = lat_km,
    mean_km = mean(c(lon_km, lat_km), na.rm = TRUE)
  )
}

#' Guess UTM EPSG from point-set center
#' @param points SpatVector points (lon/lat)
#' @return EPSG string
guess_utm_epsg <- function(points) {
  coords <- terra::crds(points)
  center_lon <- mean(coords[, 1], na.rm = TRUE)
  center_lat <- mean(coords[, 2], na.rm = TRUE)
  zone <- floor((center_lon + 180) / 6) + 1

  if (center_lat >= 0) {
    return(sprintf("EPSG:%d", 32600 + zone))
  }

  sprintf("EPSG:%d", 32700 + zone)
}

#' Thin points by climate raster cell
#' Keep at most 1 point per raster cell
#' @param points SpatVector points
#' @param clim_ref reference climate raster
#' @param keep strategy: "random" or "first"
#' @param seed random seed
#' @return list(points, kept_index, cell_ids)
thin_points_by_cell <- function(points, clim_ref, keep = "random", seed = 123) {
  stopifnot(inherits(points, "SpatVector"))
  stopifnot(inherits(clim_ref, "SpatRaster"))

  coords <- terra::crds(points)
  cell_ids <- terra::cellFromXY(clim_ref[[1]], coords)
  valid_idx <- which(!is.na(cell_ids))

  if (length(valid_idx) == 0) {
    stop("Error: No occurrence points overlap with the predictor raster extent.")
  }

  idx_pool <- valid_idx
  if (keep == "random") {
    set.seed(seed)
    idx_pool <- sample(idx_pool)
  }

  seen_cells <- integer(0)
  keep_idx <- integer(0)

  for (idx in idx_pool) {
    cell_id <- cell_ids[idx]
    if (!(cell_id %in% seen_cells)) {
      keep_idx <- c(keep_idx, idx)
      seen_cells <- c(seen_cells, cell_id)
    }
  }

  keep_idx <- sort(keep_idx)

  list(
    points = points[keep_idx, ],
    kept_index = keep_idx,
    cell_ids = cell_ids[keep_idx]
  )
}

#' Greedy thinning by minimum distance
#' @param points SpatVector points
#' @param distance_km minimum distance threshold (km)
#' @param keep strategy: "random" or "first"
#' @param seed random seed
#' @return list(points, kept_index, crs_used)
thin_points_by_distance <- function(points, distance_km, keep = "random", seed = 123) {
  stopifnot(inherits(points, "SpatVector"))

  if (terra::nrow(points) <= 1) {
    return(list(points = points, kept_index = seq_len(terra::nrow(points)), crs_used = terra::crs(points)))
  }

  utm_epsg <- guess_utm_epsg(points)
  points_proj <- terra::project(points, utm_epsg)
  coords <- terra::crds(points_proj)
  threshold_m <- distance_km * 1000

  order_idx <- seq_len(nrow(coords))
  if (keep == "random") {
    set.seed(seed)
    order_idx <- sample(order_idx)
  }

  keep_idx <- integer(0)
  kept_coords <- matrix(numeric(0), ncol = 2)

  for (idx in order_idx) {
    candidate <- coords[idx, , drop = FALSE]

    if (length(keep_idx) == 0) {
      keep_idx <- c(keep_idx, idx)
      kept_coords <- rbind(kept_coords, candidate)
      next
    }

    distances <- sqrt((kept_coords[, 1] - candidate[1, 1])^2 + (kept_coords[, 2] - candidate[1, 2])^2)

    if (all(distances >= threshold_m)) {
      keep_idx <- c(keep_idx, idx)
      kept_coords <- rbind(kept_coords, candidate)
    }
  }

  keep_idx <- sort(keep_idx)

  list(
    points = points[keep_idx, ],
    kept_index = keep_idx,
    crs_used = utm_epsg
  )
}

#' Compute nearest-neighbor distance per point (km)
#' @param points SpatVector points
#' @return numeric vector
nearest_neighbor_distances_km <- function(points) {
  if (!inherits(points, "SpatVector") || terra::nrow(points) < 2) {
    return(NA_real_)
  }

  utm_epsg <- guess_utm_epsg(points)
  points_proj <- terra::project(points, utm_epsg)
  coords <- terra::crds(points_proj)
  dist_mat <- as.matrix(stats::dist(coords)) / 1000
  diag(dist_mat) <- Inf
  apply(dist_mat, 1, min, na.rm = TRUE)
}

#' Summarize nearest-neighbor distances (km)
#' @param points SpatVector points
#' @return named numeric vector
nearest_neighbor_summary_km <- function(points) {
  nn <- nearest_neighbor_distances_km(points)

  if (length(nn) == 1 && is.na(nn)) {
    return(c(min_km = NA_real_, median_km = NA_real_, mean_km = NA_real_))
  }

  c(
    min_km = min(nn, na.rm = TRUE),
    median_km = stats::median(nn, na.rm = TRUE),
    mean_km = mean(nn, na.rm = TRUE)
  )
}

#' Ensure points have a CRS
#' If missing, backfill with the configured target_crs
#' @param points SpatVector points
#' @param default_crs default CRS
#' @return SpatVector with CRS
ensure_points_crs <- function(points, default_crs = target_crs) {
  stopifnot(inherits(points, "SpatVector"))

  current_crs <- terra::crs(points)
  if (is.null(current_crs) || !nzchar(current_crs)) {
    terra::crs(points) <- default_crs
  }

  points
}

# ============================================================================
# Paleoclimate data loading functions
# ============================================================================

#' List all available paleoclimate periods
#' @param clim_dir paleoclimate root dir (default past_clim_dir)
#' @param subdir scan subdirectories? (default past_clim_subdir)
#' @return period name vector
list_paleo_periods <- function(clim_dir = past_clim_dir, subdir = past_clim_subdir) {
  if (subdir) {
    # Scan period subdirectories
    period_dirs <- list.dirs(clim_dir, full.names = FALSE, recursive = FALSE)
    # Strip a "_masked" suffix from dir names if present
    periods <- gsub("_masked$", "", period_dirs)
    return(sort(unique(periods))) # dedupe and sort
  } else {
    # Scan filenames directly in the directory
    clim_files <- list.files(clim_dir, pattern = "\\.(tif|nc)$", full.names = FALSE)
    periods <- unique(gsub("_bio\\d+.*$", "", clim_files))
    return(sort(periods))
  }
}

#' Load paleoclimate data for a given period
#' @param period period name (e.g. "tewta_0k", "tewtb_1000k")
#' @param clim_dir paleoclimate root dir (default past_clim_dir)
#' @param subdir is period data in a subdirectory? (default past_clim_subdir)
#' @param pattern file match pattern (default "*.nc$|*.tif$")
#' @return SpatRaster
#' Handles both `period/` and `period_masked/` directory layouts
load_paleo_climate <- function(period,
                               clim_dir = past_clim_dir,
                               subdir = past_clim_subdir,
                               pattern = "\\.(nc|tif)$") {

  cat("  Loading paleoclimate for period:", period, "\n")

  if (subdir) {
    # Try both common directory naming conventions
    path_plain <- file.path(clim_dir, period)
    path_masked <- file.path(clim_dir, paste0(period, "_masked"))

    if (dir.exists(path_plain)) {
      period_dir <- path_plain
    } else if (dir.exists(path_masked)) {
      period_dir <- path_masked
    } else {
      stop("Error: Neither '", period, "' nor '", period, "_masked' directory exists in ", clim_dir)
    }

    clim_files <- list.files(period_dir, pattern = pattern, full.names = TRUE)

    if (length(clim_files) == 0) {
      stop("Error: No climate files found in ", period_dir)
    }

    cat("    Found", length(clim_files), "files in", basename(period_dir), "\n")
    clim <- terra::rast(clim_files)

  } else {
    # Direct file-read mode
    clim_files <- list.files(clim_dir,
                             pattern = paste0("^", period, "_.*", pattern),
                             full.names = TRUE)
    if (length(clim_files) == 0) stop("Error: No files found for ", period)
    clim <- terra::rast(clim_files)
  }

  # Standardize variable names: accept bio1 / bio01 / bio_1 / bio_01 etc.
  layer_names <- names(clim)
  standard_names <- gsub(".*bio_?0*(\\d+).*", "Bio\\1", layer_names, ignore.case = TRUE)
  names(clim) <- standard_names

  cat("    Standardized layers:", paste(names(clim), collapse = ", "), "\n")

  return(clim)
}

#' Get projection output directory
#' Builds and creates the directory from projection_name
#' @param type projection type ("modern" or "paleo")
#' @param period paleoclimate period name (paleo only)
#' @return output directory path
get_projection_output_dir <- function(type = c("modern", "paleo"), period = NULL) {
  type <- match.arg(type)

  if (type == "modern") {
    out_dir <- file.path(output_dirs$projection, paste0("modern_", projection_name))
  } else {
    if (is.null(period)) {
      stop("Error: period is required for paleo projection")
    }
    out_dir <- file.path(output_dirs$projection, paste0("paleo_", period, "_", projection_name))
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  return(out_dir)
}

# ============================================================================
# Print current config summary on startup
# ============================================================================

cat("\n")
cat("=================================================================\n")
cat("  Sargentodoxa SDM - Configuration Loaded\n")
cat("=================================================================\n")
cat("  Running mode: Local (de-servered repo; climate rasters via env vars)\n")
cat("  Project:", project_dir, "\n")
cat("  Species:", species_name, "\n")
cat("  Data file:", sp_file, "\n")
cat("  Modern climate:", modern_clim_dir, "\n")
cat("  Past climate:", past_clim_dir, "\n")
cat("  Past climate source:", current_past_clim, "\n")
cat("  Past climate subdirs:", past_clim_subdir, "\n")
cat("\n")
cat("  Pseudo-absence strategy:", pa_strategy, "\n")
cat("  Pseudo-absence scheme:", pa_scheme_to_use, "\n")
cat("  Variable selection: VIF threshold =", vif_threshold, "\n")
cat("  Cross-validation:", cv_strategy, " (k =", cv_k, ", rep =", cv_nb_rep, ")\n")
cat("\n")
cat("  Models:", paste(model_list, collapse = ", "), "\n")
cat("  Model selection mode:", model_selection_mode, "\n")
if (model_selection_mode == "manual") {
  cat("  Manual selected models:", paste(manual_selected_models, collapse = ", "), "\n")
} else {
  rules <- paste(metric_select, ">=", metric_select_thresh, collapse = ", ")
  cat("  Auto selection: Top", top_n_models, "models (", rules, ")\n")
}
cat("  Ensemble:", paste(ensemble_algo, collapse = ", "), "\n")
cat("\n")
cat("  Output base:", output_base, "\n")
cat("  Projection name:", projection_name, "\n")
cat("=================================================================\n")
cat("\n")

