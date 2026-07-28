# ============================================================================
# 06_model_evaluation.R - 模型评估与选择
# ============================================================================
# 功能：
#   - 基于评估指标筛选优质模型
#   - 支持手动指定模型或自动筛选
#   - 保存选中模型列表
# ============================================================================

# 加载配置
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  06 Model Evaluation - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. 加载数据
# ============================================================================

# 加载模型训练结果
model_results_file <- file.path(output_dirs$models, "model_results.RData")
if (!file.exists(model_results_file)) {
  stop("Error: Model results not found. Run 05_model_training.R first.")
}

load(model_results_file)
cat("  Loaded model results:", model_results_file, "\n")

training_metadata_file <- file.path(output_dirs$models, "model_training_metadata.rds")
training_metadata <- NULL
if (file.exists(training_metadata_file)) {
  training_metadata <- readRDS(training_metadata_file)
}

# 加载评估摘要
eval_summary_file <- file.path(output_dirs$models, "model_evaluation_summary.csv")
if (!file.exists(eval_summary_file)) {
  stop("Error: Evaluation summary not found.")
}

eval_summary <- read.csv(eval_summary_file)
cat("  Loaded evaluation summary:", eval_summary_file, "\n")

# 显示所有模型评估结果
cat("\n  All models evaluation:\n")
print(eval_summary)

# ============================================================================
# 2. 模型筛选
# ============================================================================

cat("\n  Selecting models...\n")
cat("  Selection mode:", model_selection_mode, "\n")

if (model_selection_mode == "manual") {
  # 手动模式：使用配置的模型列表
  selected_models <- manual_selected_models
  
  # 检查选中的模型是否在评估结果中
  valid_models <- selected_models[selected_models %in% eval_summary$algo]
  invalid_models <- selected_models[!selected_models %in% eval_summary$algo]
  
  if (length(invalid_models) > 0) {
    cat("  Warning: Models not found in results:", paste(invalid_models, collapse = ", "), "\n")
  }
  
  cat("  Manual selected models:", paste(valid_models, collapse = ", "), "\n")
  
} else {
  # 自动模式：根据阈值筛选
  if (length(metric_select) != length(metric_select_thresh)) {
    stop("Error: metric_select and metric_select_thresh must have the same length.")
  }
  
  cat("  Thresholds:\n")
  for (i in seq_along(metric_select)) {
    cat("   -", metric_select[i], ">=", metric_select_thresh[i], "\n")
  }
  
  qualified <- eval_summary
  
  for (i in seq_along(metric_select)) {
    metric_name <- metric_select[i]
    metric_col <- paste0("mean_", metric_name)
    metric_thresh <- metric_select_thresh[i]
    
    if (!metric_col %in% names(qualified)) {
      stop("Error: Required metric column not found in evaluation summary: ", metric_col)
    }
    
    qualified <- qualified[!is.na(qualified[[metric_col]]) & qualified[[metric_col]] >= metric_thresh, ]
  }
  
  # 按 TSS 降序排列
  qualified <- qualified[order(-qualified$mean_TSS), ]
  
  # 选择 Top N
  if (nrow(qualified) >= top_n_models) {
    selected_models <- qualified$algo[1:top_n_models]
  } else {
    selected_models <- qualified$algo
    cat("  Warning: Only", nrow(qualified), "models meet threshold\n")
  }
  
  cat("  Auto selected (Top", length(selected_models), "):", paste(selected_models, collapse = ", "), "\n")
}

# 获取选中模型的详细信息
selected_eval <- eval_summary[eval_summary$algo %in% selected_models, ]
selected_eval <- selected_eval[order(-selected_eval$mean_TSS), ]

all_built_models <- get_built_models(myBiomodModel)
selected_response_models <- unlist(lapply(selected_models, function(algo) {
  run1_models <- grep(paste0("_RUN1_", algo, "$"), all_built_models, value = TRUE)
  if (length(run1_models) > 0) {
    return(run1_models)
  }
  grep(paste0("_allRun_", algo, "$"), all_built_models, value = TRUE)
}))
selected_response_models <- unique(selected_response_models)
response_curve_variables <- if (!is.null(training_metadata) && !is.null(training_metadata$variables)) {
  training_metadata$variables
} else {
  names(get_formal_data(myBiomodModel, "expl.var"))
}

cat("\n  Selected models details:\n")
print(selected_eval)

# ============================================================================
# 3. 保存结果
# ============================================================================

cat("\n  Saving results...\n")

# 确保输出目录存在
dir.create(output_dirs$evaluation, showWarnings = FALSE, recursive = TRUE)

# 保存选中模型列表
selected_file <- file.path(output_dirs$evaluation, "selected_models.rds")
saveRDS(selected_models, selected_file)
cat("  Saved selected models:", selected_file, "\n")

# 保存模型排名
ranking_file <- file.path(output_dirs$evaluation, "model_ranking.csv")
ranking <- eval_summary[order(-eval_summary$mean_TSS), ]
ranking$rank <- 1:nrow(ranking)
write.csv(ranking, ranking_file, row.names = FALSE)
cat("  Saved model ranking:", ranking_file, "\n")

# 保存选中模型的评估详情
selected_eval_file <- file.path(output_dirs$evaluation, "selected_models_evaluation.csv")
write.csv(selected_eval, selected_eval_file, row.names = FALSE)
cat("  Saved selected evaluation:", selected_eval_file, "\n")

# ============================================================================
# 4. 可视化
# ============================================================================

cat("\n  Generating visualization...\n")

pdf_file <- file.path(output_dirs$evaluation, "model_evaluation_summary.pdf")
pdf(pdf_file, width = 14, height = 10)

# 为 metric_select 中的每个指标单独绘制全模型比较页
metric_bar_colors <- c("forestgreen", "steelblue", "darkorange", "orchid4", "brown3")
for (i in seq_along(metric_select)) {
  metric_name <- metric_select[i]
  metric_col <- paste0("mean_", metric_name)
  metric_thresh <- metric_select_thresh[i]
  bar_color <- metric_bar_colors[((i - 1) %% length(metric_bar_colors)) + 1]
  
  if (!metric_col %in% names(eval_summary)) {
    next
  }
  
  par(mar = c(8, 4, 4, 2))
  barplot(
    eval_summary[[metric_col]],
    names.arg = eval_summary$algo,
    las = 2,
    col = ifelse(eval_summary$algo %in% selected_models, bar_color, "gray"),
    main = paste("Model", metric_name, "Comparison"),
    ylab = metric_name,
    ylim = c(0, 1)
  )
  abline(h = metric_thresh, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("Selected", "Not Selected"), fill = c(bar_color, "gray"))
}

# 选中模型的指标对比（使用主要排序指标 TSS）
selected_data <- eval_summary[eval_summary$algo %in% selected_models, ]
selected_data <- selected_data[order(-selected_data$mean_TSS), ]

x <- barplot(
  selected_data$mean_TSS,
  names.arg = selected_data$algo,
  las = 2,
  col = "steelblue",
  main = "Selected Models - TSS",
  ylab = "TSS",
  ylim = c(0, 1)
)

# 添加误差棒
segments(x, selected_data$mean_TSS - selected_data$sd_TSS,
         x, selected_data$mean_TSS + selected_data$sd_TSS,
         col = "black")
segments(x - 0.2, selected_data$mean_TSS - selected_data$sd_TSS,
         x + 0.2, selected_data$mean_TSS - selected_data$sd_TSS,
         col = "black")
segments(x - 0.2, selected_data$mean_TSS + selected_data$sd_TSS,
         x + 0.2, selected_data$mean_TSS + selected_data$sd_TSS,
         col = "black")

# 最后一页：模型选择总结
plot.new()
text(0.5, 0.9, "Model Selection Summary", cex = 1.5, font = 2)
text(0.5, 0.75, paste("Selection mode:", model_selection_mode), cex = 1.2)
text(0.5, 0.65, paste("Number of models selected:", length(selected_models)), cex = 1.2)
text(0.5, 0.5, "Selected models:", cex = 1.1)
text(0.5, 0.4, paste(selected_models, collapse = ", "), cex = 1)
text(0.5, 0.28, paste("Best TSS:", round(max(selected_data$mean_TSS), 4),
                      "(", selected_data$algo[which.max(selected_data$mean_TSS)], ")"), cex = 1)
text(0.5, 0.18, paste("Best AUC:", round(max(selected_data$mean_AUCroc), 4),
                      "(", selected_data$algo[which.max(selected_data$mean_AUCroc)], ")"), cex = 1)
if ("mean_BOYCE" %in% names(selected_data)) {
  text(0.5, 0.08, paste("Best BOYCE:", round(max(selected_data$mean_BOYCE), 4),
                        "(", selected_data$algo[which.max(selected_data$mean_BOYCE)], ")"), cex = 1)
}

dev.off()

cat("  Saved visualization:", pdf_file, "\n")

# 单独输出筛选后模型的响应曲线，避免和 summary PDF 混用设备时丢图
response_pdf_file <- file.path(output_dirs$evaluation, "selected_model_response_curves.pdf")
if (length(selected_response_models) > 0) {
  pdf(response_pdf_file, width = 14, height = 10)
  tryCatch({
    bm_PlotResponseCurves(
      bm.out = myBiomodModel,
      models.chosen = selected_response_models,
      show_variables = response_curve_variables,
      fixed.var = "median",
      do.bivariate = FALSE,
      do.plot = TRUE,
      main = paste0(
        "Response Curves: Selected Models (",
        paste(selected_models, collapse = ", "),
        ")"
      )
    )
  }, error = function(e) {
    plot.new()
    text(
      0.5,
      0.5,
      paste("Selected-model response curves could not be generated:", conditionMessage(e)),
      cex = 0.9
    )
  })
  dev.off()
  cat("  Saved response curves:", response_pdf_file, "\n")
}

# ============================================================================
# 5. 保存元数据
# ============================================================================

metadata <- list(
  selection_mode = model_selection_mode,
  manual_selected = if (model_selection_mode == "manual") manual_selected_models else NULL,
  auto_metrics = if (model_selection_mode == "auto") metric_select else NULL,
  auto_thresholds = if (model_selection_mode == "auto") metric_select_thresh else NULL,
  top_n = if (model_selection_mode == "auto") top_n_models else NULL,
  selected_models = selected_models,
  selected_response_models = selected_response_models,
  response_curve_variables = response_curve_variables,
  response_curve_file = if (length(selected_response_models) > 0) response_pdf_file else NULL,
  n_selected = length(selected_models),
  all_models = eval_summary$algo,
  ranking = ranking
)

metadata_file <- file.path(output_dirs$evaluation, "model_evaluation_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  06 Model Evaluation - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Summary:\n")
cat("  - Selection mode:", model_selection_mode, "\n")
cat("  - Models selected:", length(selected_models), "\n")
cat("  - Selected:", paste(selected_models, collapse = ", "), "\n")
cat("\n")
# 清理内存
gc()
