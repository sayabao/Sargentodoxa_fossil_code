# ============================================================================
# 09_paleoclimate_projection.R - 古气候投影
# ============================================================================
# 功能：
#   - 将模型投影到古气候数据
#   - 支持分时期运行（通过 config 中的 projection_periods 控制）
#   - 输出扁平化目录结构
# ============================================================================

# 加载配置
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  09 Paleoclimate Projection - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. 加载模型和数据
# ============================================================================

# 加载单模型训练结果
model_results_file <- file.path(output_dirs$models, "model_results.RData")
if (!file.exists(model_results_file)) {
  stop("Error: Single model results not found. Run 05_model_training.R first.")
}
load(model_results_file)
cat("  Loaded single model results\n")

# 加载集成模型
ensemble_file <- file.path(output_dirs$ensemble, "ensemble_model.RData")
if (!file.exists(ensemble_file)) {
  stop("Error: Ensemble model not found. Run 07_ensemble_modeling.R first.")
}
load(ensemble_file)
cat("  Loaded ensemble model:", ensemble_file, "\n")

# 加载选中的变量
variables_file <- file.path(output_dirs$variable_select, "variables_selected.csv")
if (file.exists(variables_file)) {
  selected_vars <- read.csv(variables_file)$variable
  cat("  Selected variables:", paste(selected_vars, collapse = ", "), "\n")
} else {
  selected_vars <- NULL
  cat("  Variables not found, using all available\n")
}

# 加载格式化的数据（用于获取坐标等）
data_file_temp <- file.path(output_dirs$formatted, "biomod2_formatted_data.RData")
if (!file.exists(data_file_temp)) {
  stop("Error: Formatted data not found.")
}
load(data_file_temp)
cat("  Loaded formatted data\n")

# ============================================================================
# 2. 获取要投影的时期列表
# ============================================================================

cat("\n  Scanning available paleoclimate periods...\n")

# 使用 config 中的函数获取时期列表
all_periods <- list_paleo_periods()
cat("  Found", length(all_periods), "periods in", current_past_clim, "\n")

# 根据 config 中的 projection_periods 过滤
if (is.null(projection_periods)) {
  # 默认：投影所有时期
  periods_to_project <- all_periods
  cat("  Projection periods: ALL (", length(periods_to_project), ")\n")
} else if (is.character(projection_periods)) {
  # 指定时期列表
  periods_to_project <- intersect(projection_periods, all_periods)
  cat("  Projection periods: CUSTOM (", length(periods_to_project), ")\n")
} else if (is.list(projection_periods) && !is.null(projection_periods$time_range)) {
  # 按时间范围过滤（需要从时期名称中提取时间）
  time_range <- projection_periods$time_range
  # 提取时期名称中的数字（单位：ka）
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
# 3. 创建输出目录
# ============================================================================

# 主输出目录
output_paleo_base <- "./output/09_paleoclimate"
dir.create(output_paleo_base, showWarnings = FALSE, recursive = TRUE)

# 投影结果目录
output_proj_dir <- file.path(output_paleo_base, "projections")
dir.create(output_proj_dir, showWarnings = FALSE, recursive = TRUE)

# 统计结果目录
output_stats_dir <- file.path(output_paleo_base, "stats")
dir.create(output_stats_dir, showWarnings = FALSE, recursive = TRUE)

cat("  Output base:", output_paleo_base, "\n")

# ============================================================================
# 4. 获取单模型列表（用于投影）
# ============================================================================

selected_file <- file.path(output_dirs$evaluation, "selected_models.rds")
selected_models <- readRDS(selected_file)
all_built_models <- get_built_models(myBiomodModel)

models_to_project <- unlist(sapply(selected_models, function(algo) {
  pattern <- paste0(algo, "$")
  grep(pattern, all_built_models, value = TRUE)
}))
names(models_to_project) <- NULL

cat("  Models to project:", length(models_to_project), "\n")

# ============================================================================
# 5. 逐时期投影
# ============================================================================

cat("\n")
cat("=================================================================\n")
cat("  Starting projection for", length(periods_to_project), "periods\n")
cat("=================================================================\n")

# 存储统计结果
all_stats <- data.frame()

# 设置投影名称前缀
proj_prefix <- paste0("paleo_", projection_name)

# 进度跟踪
n_periods <- length(periods_to_project)

for (i in seq_along(periods_to_project)) {
  period <- periods_to_project[i]
  
  cat("\n  [", i, "/", n_periods, "] Processing:", period, "\n")
  
  # ----- 5.1 加载该时期的气候数据 -----
  load_error <- FALSE
  
  tryCatch({
    paleo_clim <- load_paleo_climate(period)
    
    # 只保留选中的变量
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
  
  # 跳过加载失败的时期
  if (load_error) {
    cat("    Skipping this period due to error\n")
    next
  }
  
  # ----- 5.2 单模型投影 -----
  proj_name <- paste0(proj_prefix, "_", period)
  
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
  
  # ----- 5.3 集成模型预测 -----
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
  
  # ----- 5.4 提取并保存统计信息 -----
  cat("    Extracting statistics...")
  
  sp_name <- myBiomodModel@sp.name
  proj_dir <- file.path(getwd(), "Models", sp_name, paste0("proj_", proj_name), "individual_projections")
  
  if (dir.exists(proj_dir)) {
    proj_files <- list.files(proj_dir, pattern = "\\.tif$", full.names = TRUE)
    
    # 提取 EMwmean 的结果（推荐）
    em_files <- proj_files[grepl("EMwmean", proj_files) & !grepl("_TSSbin", proj_files)]
    
    if (length(em_files) > 0) {
      # 复制到输出目录
      for (f in em_files) {
        new_name <- paste0(period, "_", basename(f))
        file.copy(f, file.path(output_proj_dir, new_name), overwrite = TRUE)
      }
      
      # 计算统计
      r <- rast(em_files[1])
      stats <- data.frame(
        period = period,
        min = global(r, "min", na.rm = TRUE)[1, 1],
        max = global(r, "max", na.rm = TRUE)[1, 1],
        mean = global(r, "mean", na.rm = TRUE)[1, 1],
        sd = global(r, "sd", na.rm = TRUE)[1, 1],
        na_count = global(is.na(r), "sum", na.rm = TRUE)[1, 1]
      )
      
      # 保存单时期统计
      stats_file <- file.path(output_stats_dir, paste0(period, "_summary.csv"))
      write.csv(stats, stats_file, row.names = FALSE)
      
      # 添加到汇总
      all_stats <- rbind(all_stats, stats)
    }
  }
  
  cat(" done\n")
  
  # ----- 5.5 内存清理 -----
  if (i %% gc_trigger == 0) {
    gc()
    cat("    [Garbage collected]\n")
  }
  
  # 清理当前时期的模型对象，释放内存
  rm(myBiomodProj, myBiomodEF, paleo_clim)
}

# ============================================================================
# 6. 保存汇总统计
# ============================================================================

cat("\n  Collecting all available statistics files...\n")

# 1. 获取 stats 目录下所有的汇总 CSV 文件
stat_files <- list.files(output_stats_dir, pattern = "_summary.csv$", full.names = TRUE)

if (length(stat_files) > 0) {
  cat("  Found", length(stat_files), "individual period stats. Merging...\n")
  
  # 2. 读取并合并所有文件
  all_stats <- do.call(rbind, lapply(stat_files, function(f) {
    read.csv(f, stringsAsFactors = FALSE)
  }))
  
  # 3. 按时期名称排序（可选，根据你的时期命名规则）
  # all_stats <- all_stats[order(all_stats$period), ] 
  
  # 4. 保存最终汇总文件
  summary_file <- file.path(output_paleo_base, "summary_all_periods.csv")
  write.csv(all_stats, summary_file, row.names = FALSE)
  cat("  Saved final summary:", summary_file, "\n")
  
  cat("\n  Projection summary:\n")
  print(head(all_stats, 10))
} else {
  cat("  Warning: No statistics files found in", output_stats_dir, "\n")
}

# ============================================================================
# 7. 元数据
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

# 清理
gc()