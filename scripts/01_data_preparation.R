# ============================================================================
# 01_data_preparation.R - 数据准备与验证
# ============================================================================
# 功能：
#   - 读取已清洗的 CVH 物种分布数据
#   - 验证数据完整性
#   - 转换为 SpatVector 格式
#   - 保存到 output 目录
# ============================================================================

# 加载配置
source("./config.R")

cat("\n")
cat("=================================================================\n")
cat("  01 Data Preparation - Starting\n")
cat("=================================================================\n")
cat("  Reading CVH data from:", sp_file_cvh, "\n")
cat("\n")

# ============================================================================
# 1. 读取已清洗的 CVH 数据
# ============================================================================

if (!file.exists(sp_file_cvh)) {
  stop("Error: CVH data file not found: ", sp_file_cvh)
}

cvh_data <- read.csv(sp_file_cvh, stringsAsFactors = FALSE)

cat("  Data loaded successfully\n")
cat("  Total records:", nrow(cvh_data), "\n")
cat("  Columns:", paste(colnames(cvh_data), collapse = ", "), "\n")

# ============================================================================
# 2. 数据验证
# ============================================================================

# 检查必要列
required_cols <- c("Longitude", "Latitude", "presence")
missing_cols <- setdiff(required_cols, colnames(cvh_data))

if (length(missing_cols) > 0) {
  stop("Error: Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("\n  Column validation: PASSED\n")

# 检查坐标范围
lon_valid <- cvh_data$Longitude >= -180 & cvh_data$Longitude <= 180
lat_valid <- cvh_data$Latitude >= -90 & cvh_data$Latitude <= 90

if (!all(lon_valid)) {
  cat("  Warning: Some Longitude values are out of range [-180, 180]\n")
  cat("  Out-of-range records:", sum(!lon_valid), "\n")
}

if (!all(lat_valid)) {
  cat("  Warning: Some Latitude values are out of range [-90, 90]\n")
  cat("  Out-of-range records:", sum(!lat_valid), "\n")
}

if (all(lon_valid) && all(lat_valid)) {
  cat("  Coordinate range validation: PASSED\n")
}

# 检查 NA 值
na_records <- apply(cvh_data, 1, function(x) any(is.na(x)))
if (any(na_records)) {
  cat("  Warning: Found", sum(na_records), "records with NA values\n")
} else {
  cat("  NA value check: PASSED\n")
}

# ============================================================================
# 3. 统计信息
# ============================================================================

cat("\n  Statistics:\n")
cat("  - Total presence points:", sum(cvh_data$presence == 1), "\n")

# 计算分布范围
lon_range <- range(cvh_data$Longitude, na.rm = TRUE)
lat_range <- range(cvh_data$Latitude, na.rm = TRUE)
cat("  - Longitude range:", lon_range[1], "to", lon_range[2], "\n")
cat("  - Latitude range:", lat_range[1], "to", lat_range[2], "\n")

# ============================================================================
# 4. 转换为 SpatVector
# ============================================================================

presence_spat <- vect(cvh_data, geom = c("Longitude", "Latitude"))
presence_spat <- ensure_points_crs(presence_spat)

cat("\n  Converted to SpatVector\n")
cat("  - Geometry type:", geomtype(presence_spat), "\n")
cat("  - CRS:", crs(presence_spat), "\n")

# ============================================================================
# 4.1 空间稀疏化（降低空间自相关）
# ============================================================================

original_n <- nrow(cvh_data)
duplicate_n <- sum(duplicated(cvh_data[, c("Longitude", "Latitude")]))
cvh_data <- cvh_data[!duplicated(cvh_data[, c("Longitude", "Latitude")]), ]
presence_spat <- vect(cvh_data, geom = c("Longitude", "Latitude"))
presence_spat <- ensure_points_crs(presence_spat)
presence_spat_before_thinning <- presence_spat

cat("\n  Spatial thinning preparation:\n")
cat("  - Removed exact duplicate coordinates:", duplicate_n, "\n")
cat("  - Records after deduplication:", nrow(cvh_data), "\n")

nn_before <- nearest_neighbor_summary_km(presence_spat)
nn_before_all <- nearest_neighbor_distances_km(presence_spat)
cell_thinning_n <- nrow(cvh_data)
distance_thinning_n <- nrow(cvh_data)
thinning_distance_used_km <- NA_real_
thinning_resolution_km <- NA_real_
thinning_crs_used <- NA_character_
thinning_applied <- FALSE

if (enable_spatial_thinning) {
  cat("\n  Spatial thinning: ENABLED\n")
  cat("  - Method:", thinning_method, "\n")
  
  current_points <- presence_spat
  
  if (thinning_use_predictor_cell && thinning_method %in% c("cell", "hybrid")) {
    modern_clim <- load_standard_climate(modern_clim_dir, modern_clim_pattern)
    thinning_resolution <- estimate_raster_resolution_km(modern_clim, current_points)
    thinning_resolution_km <- thinning_resolution$mean_km
    
    cell_result <- thin_points_by_cell(
      points = current_points,
      clim_ref = modern_clim,
      keep = thinning_keep,
      seed = thinning_seed
    )
    current_points <- cell_result$points
    cell_thinning_n <- terra::nrow(current_points)
    
    cat("  - Predictor-cell thinning applied\n")
    cat("  - Estimated raster resolution:", round(thinning_resolution_km, 2), "km\n")
    cat("  - Records after cell thinning:", cell_thinning_n, "\n")
  }
  
  if (thinning_method %in% c("distance", "hybrid")) {
    if (thinning_distance_mode == "auto") {
      if (is.na(thinning_resolution_km)) {
        modern_clim <- load_standard_climate(modern_clim_dir, modern_clim_pattern)
        thinning_resolution <- estimate_raster_resolution_km(modern_clim, current_points)
        thinning_resolution_km <- thinning_resolution$mean_km
      }
      thinning_distance_used_km <- thinning_resolution_km * thinning_distance_factor
    } else {
      thinning_distance_used_km <- thinning_distance_km
    }
    
    distance_result <- thin_points_by_distance(
      points = current_points,
      distance_km = thinning_distance_used_km,
      keep = thinning_keep,
      seed = thinning_seed
    )
    current_points <- distance_result$points
    distance_thinning_n <- terra::nrow(current_points)
    thinning_crs_used <- distance_result$crs_used
    
    cat("  - Distance thinning applied\n")
    cat("  - Minimum distance:", round(thinning_distance_used_km, 2), "km\n")
    cat("  - Distance CRS:", thinning_crs_used, "\n")
    cat("  - Records after distance thinning:", distance_thinning_n, "\n")
  } else {
    thinning_distance_used_km <- NA_real_
    thinning_crs_used <- NA_character_
  }
  
  presence_spat <- current_points
  cvh_data <- as.data.frame(presence_spat, geom = "XY")
  thinning_applied <- TRUE
} else {
  cat("\n  Spatial thinning: DISABLED\n")
}

nn_after <- nearest_neighbor_summary_km(presence_spat)
nn_after_all <- nearest_neighbor_distances_km(presence_spat)
cat("  - Final presence points:", nrow(cvh_data), "\n")
cat("  - Nearest-neighbor distance before thinning (median):", round(nn_before["median_km"], 2), "km\n")
cat("  - Nearest-neighbor distance after thinning (median):", round(nn_after["median_km"], 2), "km\n")

diagnostic_thresholds <- sort(unique(thinning_diagnostic_thresholds_km))
threshold_counts_before <- sapply(diagnostic_thresholds, function(th) sum(nn_before_all < th, na.rm = TRUE))
threshold_counts_after <- sapply(diagnostic_thresholds, function(th) sum(nn_after_all < th, na.rm = TRUE))

cat("  - Nearby pairs by threshold before thinning:\n")
for (i in seq_along(diagnostic_thresholds)) {
  cat("    <", diagnostic_thresholds[i], "km:", threshold_counts_before[i], "points\n")
}
cat("  - Nearby pairs by threshold after thinning:\n")
for (i in seq_along(diagnostic_thresholds)) {
  cat("    <", diagnostic_thresholds[i], "km:", threshold_counts_after[i], "points\n")
}

# ============================================================================
# 5. 保存数据
# ============================================================================

# 确保输出目录存在
dir.create(output_dirs$data_prep, showWarnings = FALSE, recursive = TRUE)

# 保存 CSV
output_csv <- file.path(output_dirs$data_prep, "sargentodoxa_occ_cvh_clean.csv")
write.csv(cvh_data, output_csv, row.names = FALSE)
cat("\n  Saved CSV:", output_csv, "\n")

# 保存 SpatVector
output_rds <- file.path(output_dirs$data_prep, "sargentodoxa_presence_spat.rds")
saveRDS(presence_spat, output_rds)
cat("  Saved SpatVector:", output_rds, "\n")

# 保存稀疏化统计
thinning_stats <- data.frame(
  stage = c("original", "deduplicated", "cell_thinned", "distance_thinned", "final"),
  n_points = c(
    original_n,
    original_n - duplicate_n,
    cell_thinning_n,
    distance_thinning_n,
    nrow(cvh_data)
  )
)

thinning_stats_file <- file.path(output_dirs$data_prep, "presence_thinning_summary.csv")
write.csv(thinning_stats, thinning_stats_file, row.names = FALSE)
cat("  Saved thinning summary:", thinning_stats_file, "\n")

thinning_diag_file <- file.path(output_dirs$data_prep, "presence_thinning_diagnostics.csv")
thinning_diag <- data.frame(
  threshold_km = diagnostic_thresholds,
  n_points_before = threshold_counts_before,
  n_points_after = threshold_counts_after,
  pct_points_before = round(threshold_counts_before / original_n * 100, 2),
  pct_points_after = round(threshold_counts_after / nrow(cvh_data) * 100, 2)
)
write.csv(thinning_diag, thinning_diag_file, row.names = FALSE)
cat("  Saved thinning diagnostics:", thinning_diag_file, "\n")

# 保存稀疏化前后对比图
thinning_plot_file <- file.path(output_dirs$data_prep, "presence_thinning_check.pdf")
pdf(thinning_plot_file, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(
  presence_spat_before_thinning,
  pch = 16,
  col = "#1f77b4",
  main = "Before thinning",
  axes = TRUE
)
plot(
  presence_spat,
  pch = 16,
  col = "#d62728",
  main = "After thinning",
  axes = TRUE
)
dev.off()
cat("  Saved thinning plot:", thinning_plot_file, "\n")

thinning_hist_file <- file.path(output_dirs$data_prep, "presence_thinning_nn_histograms.pdf")
zoom_xlim <- thinning_histogram_xlim_km
nn_before_zoom <- nn_before_all[nn_before_all <= zoom_xlim]
nn_after_zoom <- nn_after_all[nn_after_all <= zoom_xlim]
breaks_zoom <- seq(0, zoom_xlim, by = 5)
hist_before_counts <- hist(nn_before_zoom, breaks = breaks_zoom, plot = FALSE)$counts
hist_after_counts <- hist(nn_after_zoom, breaks = breaks_zoom, plot = FALSE)$counts
ymax_zoom <- max(c(hist_before_counts, hist_after_counts), na.rm = TRUE)

pdf(thinning_hist_file, width = 10, height = 10)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(
  nn_before_zoom,
  breaks = breaks_zoom,
  col = "#1f77b4",
  main = paste0("NN distance before thinning (0-", zoom_xlim, " km)"),
  xlab = "Nearest-neighbor distance (km)",
  xlim = c(0, zoom_xlim),
  ylim = c(0, ymax_zoom)
)
abline(v = diagnostic_thresholds, col = "firebrick", lty = 2)
mtext(paste0(">", zoom_xlim, " km: ", sum(nn_before_all > zoom_xlim, na.rm = TRUE), " points"), side = 3, line = 0.2, cex = 0.8)

hist(
  nn_after_zoom,
  breaks = breaks_zoom,
  col = "#d62728",
  main = paste0("NN distance after thinning (0-", zoom_xlim, " km)"),
  xlab = "Nearest-neighbor distance (km)",
  xlim = c(0, zoom_xlim),
  ylim = c(0, ymax_zoom)
)
abline(v = diagnostic_thresholds, col = "firebrick", lty = 2)
mtext(paste0(">", zoom_xlim, " km: ", sum(nn_after_all > zoom_xlim, na.rm = TRUE), " points"), side = 3, line = 0.2, cex = 0.8)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(
  ecdf(nn_before_all),
  verticals = TRUE,
  do.points = FALSE,
  col = "#1f77b4",
  lwd = 2,
  main = "NN distance cumulative distribution",
  xlab = "Nearest-neighbor distance (km)",
  ylab = "Cumulative proportion",
  xlim = c(0, zoom_xlim)
)
lines(
  ecdf(nn_after_all),
  verticals = TRUE,
  do.points = FALSE,
  col = "#d62728",
  lwd = 2
)
abline(v = diagnostic_thresholds, col = "firebrick", lty = 2)
legend(
  "bottomright",
  legend = c("Before thinning", "After thinning"),
  col = c("#1f77b4", "#d62728"),
  lwd = 2,
  bty = "n"
)
dev.off()
cat("  Saved NN histogram:", thinning_hist_file, "\n")

# ============================================================================
# 6. 保存元数据
# ============================================================================

metadata <- list(
  input_file = sp_file_cvh,
  total_records = nrow(cvh_data),
  original_records = original_n,
  duplicate_records_removed = duplicate_n,
  presence_records = sum(cvh_data$presence == 1),
  lon_range = lon_range,
  lat_range = lat_range,
  crs = crs(presence_spat),
  thinning = list(
    enabled = enable_spatial_thinning,
    applied = thinning_applied,
    method = thinning_method,
    keep = thinning_keep,
    predictor_cell = thinning_use_predictor_cell,
    distance_mode = thinning_distance_mode,
    distance_km = thinning_distance_used_km,
    predictor_resolution_km = thinning_resolution_km,
    distance_crs = thinning_crs_used,
    nearest_neighbor_before_km = nn_before,
    nearest_neighbor_after_km = nn_after,
    diagnostic_thresholds_km = diagnostic_thresholds,
    histogram_xlim_km = zoom_xlim,
    diagnostic_file = thinning_diag_file,
    nearest_neighbor_histogram = thinning_hist_file,
    summary_file = thinning_stats_file,
    check_plot = thinning_plot_file
  )
)

metadata_file <- file.path(output_dirs$data_prep, "data_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  01 Data Preparation - Completed\n")
cat("=================================================================\n")
cat("\n")
