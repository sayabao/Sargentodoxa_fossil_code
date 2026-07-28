# ============================================================================
# 10_visualization.R - 结果可视化
# ============================================================================
# 功能：
#   - 可视化古气候投影结果
#   - 生成静态图和动画
#   - 对比化石记录与模型预测
# ============================================================================

# 加载配置
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  10 Visualization - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. 设置输出目录
# ============================================================================

viz_dir <- output_dirs$visualization
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(viz_dir, "static"), showWarnings = FALSE)
dir.create(file.path(viz_dir, "animation"), showWarnings = FALSE)
dir.create(file.path(viz_dir, "comparison"), showWarnings = FALSE)

cat("  Output directory:", viz_dir, "\n")

# ============================================================================
# 2. 加载古气候投影结果
# ============================================================================

paleo_proj_dir <- output_dirs$paleoclimate

# 检查输出目录
if (!dir.exists(paleo_proj_dir)) {
  stop("Error: Paleoclimate projection output not found. Run 09 first.")
}

# 获取所有 tif 文件
tif_files <- list.files(
  file.path(paleo_proj_dir, "projections"),
  pattern = "\\.tif$",
  full.names = TRUE
)

if (length(tif_files) == 0) {
  stop("Error: No projection TIF files found.")
}

# 一个时期通常会输出多张集成图（例如 EMwmeanByTSS / AUCroc / AUCprg / BOYCE）。
# 可视化阶段只保留每个时期的一张代表性结果，避免同一时期名称重复导致栅格图层重名。
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

# --- 动态提取时期名称 (不再硬编码 "_Sargentodoxa") ---
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

# 加载投影结果为栅格栈
proj_stack <- rast(tif_files)
proj_stack <- normalize_prediction_raster(proj_stack)

# --- 关键修复：在这里手动为 proj_stack 的图层命名 ---
names(proj_stack) <- period_names
# --- 修复结束 ---

cat("  Loaded projection stack:", nrow(proj_stack), "x", ncol(proj_stack), 
    "x", nlyr(proj_stack), "\n")

# (可选) 添加一个验证步骤，确保名称已正确设置
cat("  Verifying layer names. First 5 names in stack are:\n")
print(head(names(proj_stack)))


# 加载统计摘要
stats_dir <- file.path(paleo_proj_dir, "stats")
if (dir.exists(stats_dir)) {
  stats_files <- list.files(stats_dir, pattern = "summary\\.csv$", full.names = TRUE)
  stats_all <- lapply(stats_files, read.csv)
  names(stats_all) <- gsub("_summary\\.csv$", "", basename(stats_files))
  cat("  Loaded", length(stats_files), "stat files\n")
}

# ============================================================================
# 3. 加载现代投影结果（用于对比）
# ============================================================================

modern_proj_file <- file.path(output_dirs$projection, "projection_modern_map.pdf")
if (file.exists(modern_proj_file)) {
  cat("  Modern projection found:", modern_proj_file, "\n")
}

# ============================================================================
# 4. 加载化石记录（用于验证）
# ============================================================================

# 检查化石记录文件
fossil_data <- NULL
if (file.exists(fossil_record_file)) {
  fossil_data <- read.csv(fossil_record_file)
  # 检查必要的列：paleoLong, paleoLat, age
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
# 5. 可视化函数 (全面支持真实 Ma 年代)
# ============================================================================

# 5.1 单时期分布图
plot_single_period <- function(r, period_name, output_file, 
                               period_label = period_name, # <-- 新增：接收真实的 Ma 标签
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
      title = paste(species_name, "-", period_label), # <-- 关键：标题使用真实 Ma 标签
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

# 5.2 多时期对比图 (内存极致优化版)
plot_multi_period_comparison <- function(proj_stack, periods, output_file, 
                                         period_labels, ncol = 4, 
                                         downsample_factor = 3) { # <-- 新增降采样参数
  plot_data <- list()
  
  cat("    Extracting and downsampling raster layers (factor=", downsample_factor, ")...\n", sep="")
  
  for (i in seq_along(periods)) {
    if (periods[i] %in% names(proj_stack)) {
      r <- proj_stack[[periods[i]]]
      
      # --- 💡核心内存优化：降低空间分辨率 ---
      # 如果 factor = 3，数据量直接变成原来的 1/9！
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
  
  # 根据面板数量动态调整图片高度 (一行大概占 2.5 inch)
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
  
  # --- 💡主动释放内存 ---
  rm(combined_df, plot_data)
  gc()
  # ----------------------
  
  return(p)
}

# 5.3 栖息地面积变化图 (数值型真实年代 X 轴)
plot_area_trend <- function(stats_all, period_order, period_ma_numeric, output_file) { # <-- 新增：数值字典
  
  valid_periods <- intersect(period_order, names(stats_all))
  ordered_stats <- stats_all[valid_periods]
  
  area_data <- data.frame(
    period = valid_periods,
    mean = sapply(ordered_stats, function(x) x$mean),
    sd = sapply(ordered_stats, function(x) x$sd)
  )
  
  # --- 关键：将 X 轴映射为真实的数值年代 ---
  area_data$Ma <- period_ma_numeric[valid_periods]
  
  p <- ggplot(area_data, aes(x = Ma, y = mean)) +  # <-- 使用数值 Ma 作为 X 轴
    geom_line(color = "steelblue", linewidth = 1) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, fill = "steelblue") +
    geom_point(color = "darkblue", size = 2) +
    scale_x_reverse() + # --- 关键：X轴反转，97.2 在左边，0 在右边 ---
    labs(
      title = paste0(species_name, " - Habitat Suitability Trend"),
      subtitle = "Trend over geological time",
      x = "Age (Ma)", # <-- 坐标轴改名为 Age (Ma)
      y = "Mean Suitability Score"
    ) +
    theme_minimal()
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  return(p)
}

# 5.4 单时期化石记录对比图
plot_fossil_comparison <- function(proj_stack, period_name, fossil_data, 
                                   period_ma, output_file, threshold = 0.5) {
  # 获取对应时期的栅格
  if (!period_name %in% names(proj_stack)) {
    cat("    Warning: period", period_name, "not found in proj_stack\n")
    return(NULL)
  }
  r <- proj_stack[[period_name]]
  
  # 筛选该时期的化石记录
  fossil_subset <- fossil_data[fossil_data$age == period_ma, ]
  if (nrow(fossil_subset) == 0) {
    cat("    Warning: no fossil records for period", period_name, "at", period_ma, "Ma\n")
    return(NULL)
  }
  
  # 转换为 SpatVector (使用 paleoLong, paleoLat)
  fossil_vect <- vect(fossil_subset, geom = c("paleoLong", "paleoLat"), crs = "EPSG:4326")
  
  # 转换为数据框用于 ggplot
  r_df <- as.data.frame(r, xy = TRUE)
  names(r_df)[3] <- "value"
  
  # 获取时期标签
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

# 5.5 多时期化石对比拼图 (终极修复版 - 化石点精准匹配)
plot_fossil_multi_comparison <- function(proj_stack, period_names, fossil_data,
                                         period_ma_numeric, output_file, 
                                         ncol = 5, downsample_factor = 4) {
  
  # --- 1. 数据准备阶段 ---
  
  # 准备两个列表，一个装栅格数据，一个装化石数据
  raster_data_list <- list()
  fossil_data_list <- list()
  
  # 只保留有化石记录的时期
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
  
  # 循环处理每个有效时期
  for (i in seq_along(valid_periods)) {
    p_name <- valid_periods[i]
    ma_val <- period_ma_numeric[p_name]
    label <- period_labels[p_name]
    
    # 准备栅格数据
    r <- proj_stack[[p_name]]
    if (downsample_factor > 1) {
      r <- aggregate(r, fact = downsample_factor)
    }
    r_df <- as.data.frame(r, xy = TRUE)
    names(r_df)[3] <- "value"
    r_df$period_label <- factor(label, levels = period_labels[valid_periods])
    raster_data_list[[i]] <- r_df
    
    # --- 关键：准备对应时期的化石数据 ---
    fossil_subset <- fossil_data[fossil_data$age == ma_val, ]
    if(nrow(fossil_subset) > 0) {
      fossil_subset$period_label <- factor(label, levels = period_labels[valid_periods])
      fossil_data_list[[i]] <- fossil_subset
    }
  }
  
  # 合并成两个大数据框
  combined_raster_df <- do.call(rbind, raster_data_list)
  combined_fossil_df <- do.call(rbind, fossil_data_list)
  
  # --- 2. 绘图阶段 ---
  
  p <- ggplot() +
    # 第一层：画降采样后的栅格底图
    geom_tile(data = combined_raster_df, aes(x = x, y = y, fill = value)) +
    
    # 第二层：画精准匹配的化石点！
    geom_point(data = combined_fossil_df, aes(x = paleoLong, y = paleoLat), 
               color = "red", size = 2, shape = 1, stroke = 1.2) +
    
    # 第三层：切分面板
    facet_wrap(~ period_label, ncol = ncol) +
    
    # 第四层：美化
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
  
  # 根据面板数量动态调整图片高度
  dynamic_height <- max(8, 2.5 * ceiling(length(valid_periods) / ncol))
  
  ggsave(output_file, p, width = 16, height = dynamic_height, dpi = 300)
  cat("    Created:", basename(output_file), "\n")
  return(p)
}

# 5.6 化石点 suitability 柱状图
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

# 5.7 化石点 vs 背景 suitability 分布对比
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
# 6. 生成可视化
# ============================================================================

cat("\n  Generating visualizations...\n")

# 6.1 单时期分布图
cat("  Generating high-resolution static maps for all periods...\n")
for (period in period_names) {
  output_file <- file.path(viz_dir, "static", paste0(period, "_distribution.png"))
  if (!file.exists(output_file)) {
    if (period %in% names(proj_stack)) {
      r <- proj_stack[[period]]
      # --- 传参：将对应时期的年代标签传进去 ---
      plot_single_period(r, period, output_file, period_label = period_labels[period])
      cat("    Created:", basename(output_file), "\n")
    }
  }
}
cat("  All static maps are up to date.\n")

# 6.2 多时期对比图 (全时期展示)
if (length(period_names) >= 4) {
  output_file <- file.path(viz_dir, "comparison", "multi_period_comparison.png")
  
  cat("  Generating multi-period comparison map...\n")
  
  # downsample_factor = 3 或 4 可以极大地拯救内存
  plot_multi_period_comparison(proj_stack, period_names[1:length(period_names)], 
                               output_file, ncol = 4, 
                               period_labels = period_labels,
                               downsample_factor = 3) # <-- 调整此参数来平衡速度与清晰度
  
  cat("    Created: multi_period_comparison.png\n")
}

# 6.3 栖息地面积变化趋势
if (exists("stats_all")) {
  output_file <- file.path(viz_dir, "static", "area_trend.png")
  # --- 传参：把提取好的年代数值字典传进去 ---
  plot_area_trend(stats_all, period_names, period_ma_numeric, output_file)
  cat("    Created: area_trend.png\n")
}

# 6.4 与化石记录对比
if (!is.null(fossil_data) && "paleoLong" %in% names(fossil_data) && "paleoLat" %in% names(fossil_data)) {
  cat("\n  Generating fossil comparison plots...\n")
  
  # 6.4.1 单时期化石对比图
  for (period in period_names) {
    ma_val <- period_ma_numeric[period]
    if (is.na(ma_val)) next
    
    # 检查该时期是否有化石记录
    fossil_subset <- fossil_data[fossil_data$age == ma_val, ]
    if (nrow(fossil_subset) == 0) next
    
    output_file <- file.path(viz_dir, "comparison", paste0("fossil_comparison_", period, ".png"))
    if (!file.exists(output_file)) {
      plot_fossil_comparison(proj_stack, period, fossil_data, ma_val, output_file)
    }
  }
  
  # 6.4.2 多时期化石对比拼图
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
# 7. 生成 GIF 动画 (从高分辨率文件降维)
# ============================================================================
cat("\n  Generating animation from high-resolution files...\n")

if (require(magick)) {
  
  # 1. 按正确的顺序获取所有高分辨率PNG文件的路径
  high_res_pngs <- file.path(viz_dir, "static", paste0(period_names, "_distribution.png"))
  high_res_pngs <- high_res_pngs[file.exists(high_res_pngs)]
  
  if (length(high_res_pngs) > 0) {
    
    cat("  Found", length(high_res_pngs), "images to process.\n")
    cat("  Resizing images for animation (this may take a moment)...\n")
    
    # 2. 读取每个高分辨率文件，将其缩小，然后存入一个列表
    #    image_scale 的 "800x" 语法会按比例将宽度缩放到800像素
    images_resized <- image_read(high_res_pngs) %>%
      image_scale("600x")
    
    # 打印新尺寸以供验证
    cat("  Images resized to:\n")
    print(image_info(images_resized)[, 1:3])
    
    # 3. 将缩小后的图像列表合并成动画
    animation <- image_animate(images_resized, fps = 2) # fps=2 表示每秒2帧
    
    # 4. 保存GIF
    output_file <- file.path(viz_dir, "animation", "paleoclimate_timeline.gif")
    image_write(animation, output_file)
    cat("    Successfully created animation:", basename(output_file), "\n")
    
    # 5. 清理内存 (重要！)
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
# 8. 保存元数据
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
