# ============================================================================
# 05_model_training.R - 单模型训练
# ============================================================================
# 功能：
#   - 使用 biomod2 训练多个物种分布模型
#   - 进行 k-fold 交叉验证
#   - 评估模型性能
# ============================================================================

# 加载配置
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  05 Model Training - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. 加载数据
# ============================================================================

# 加载格式化数据
formatted_file <- file.path(output_dirs$formatted, "biomod2_formatted_data.RData")
if (!file.exists(formatted_file)) {
  stop("Error: Formatted data not found. Run 04_data_formatting.R first.")
}

load(formatted_file)
cat("  Loaded formatted data:", formatted_file, "\n")

# 读取元数据获取变量信息
metadata_file <- file.path(output_dirs$formatted, "data_formatting_metadata.rds")
metadata <- readRDS(metadata_file)

cat("  Training with", metadata$n_variables, "variables\n")
cat("  Variables:", paste(metadata$variables, collapse = ", "), "\n")

cv_table <- NULL
cv_info <- NULL
if (identical(cv_strategy, "user.defined")) {
  cv_table_file <- file.path(output_dirs$formatted, "spatial_cv_table.csv")
  cv_info_file <- file.path(output_dirs$formatted, "spatial_cv_info.rds")
  
  if (!file.exists(cv_table_file)) {
    stop("Error: Spatial CV table not found. Run 04_data_formatting.R first.")
  }
  
  cv_table <- read.csv(cv_table_file, check.names = FALSE)
  cv_table[] <- lapply(cv_table, as.logical)
  
  if (file.exists(cv_info_file)) {
    cv_info <- readRDS(cv_info_file)
  }
  
  cat("  Loaded custom spatial CV table:", cv_table_file, "\n")
  cat("  Custom CV columns:", paste(colnames(cv_table), collapse = ", "), "\n")
  if (!is.null(cv_info)) {
    cat("  Spatial CV method:", cv_info$method, "\n")
    cat("  Spatial block size:", cv_info$block_size_km, "km\n")
  }
}
# formatted_data
# plot(formatted_data)
# ============================================================================
# 2. 配置模型参数
# ============================================================================

cat("\n  Configuring model parameters...\n")

# 检查哪些模型可用
available_models <- model_list[sapply(model_list, function(m) {
  tryCatch({
    # --- 新增：专门为 MAXENT (Java 版) 做的检查 ---
    if (m == "MAXENT") {
      # 检查1: 工作目录下是否有maxent.jar
      if (!file.exists(paste0(project_dir,"/maxent.jar"))) {
        cat("   - MAXENT check: maxent.jar not found at path:", project_dir, ". Skipping.\n")
        return(FALSE)
      }
      # 检查2 (可选但推荐): 检查 Java 是否安装
      java_check <- system("/usr/bin/java --version", intern = TRUE, ignore.stderr = TRUE)
      if (length(java_check) == 0) {
        cat("   - MAXENT check: Java command not found. Skipping.\n")
        return(FALSE)
      }
      # 如果以上三项都通过，才认为 MAXENT 可用
      return(TRUE)
    } 
    # --- MAXENT 检查结束 ---
    else if (m == "MAXNET") {
      require(maxnet)
    } else if (m == "MARS") {
      require(earth)
    } else if (m == "ANN") {
      require(nnet)
    } else if (m == "CTA") {
      require(rpart)
    } else if (m == "FDA") {
      require(mda)
    } else if (m == "RF") {
      require(randomForest)
    } else if (m == "RFd") {
      require(randomForest)
    } else if (m == "GAM") {
      require(gam)
    } else if (m == "GBM") {
      require(gbm)
    } else if (m == "GLM") {
      TRUE  # GLM 是基础 R 的一部分，永远为 TRUE
    } else if (m == "SRE") {
      TRUE # SRE 是 biomod2 内置的，也永远为 TRUE
    }
    else {
      FALSE # 对于其他未知的模型，返回 FALSE
    }
  }, error = function(e) {
    # 如果 require(package) 失败，会抛出错误，我们在这里捕获它
    cat("   - Check for", m, "failed (package likely not installed).\n")
    return(FALSE)
  })
})]

# 模型选项
bm.options <- bm_ModelingOptions(
  data.type = "binary",          # 数据类型：物种分布（出现/伪缺席）属于二元分类 "binary"
  models = available_models,     # 传入刚刚检查好的可用模型列表
  strategy = opt_strategy,           # 策略："default" 代表使用默认参数（相当于旧版的空号）
  bm.format = formatted_data
  # 你也可以改用 "bigboss"（biomod2 团队推荐的预设优选参数）或 "tuned"（自动调参，但很慢）
)

cat("  Available models:", paste(available_models, collapse = ", "), "\n")
cat("  CV strategy:", cv_strategy, "\n")
if (identical(cv_strategy, "user.defined")) {
  custom_fold_cols <- setdiff(colnames(cv_table), "_allData_allRun")
  cat("  Custom CV folds:", length(custom_fold_cols), "\n")
} else {
  cat("  CV repeats:", cv_nb_rep, "\n")
  cat("  CV folds:", cv_k, "\n")
}
cat("  Evaluation metrics:", paste(metric_eval, collapse = ", "), "\n")

# ============================================================================
# 3. 训练模型
# ============================================================================

cat("\n  Training models...\n")

# 确保输出目录存在
models_dir <- "./Models"
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

# 设置 modeling.id
modeling.id <- paste0(species_name, "_", format(Sys.time(), "%Y%m%d%H%M"))

modeling_args <- list(
  bm.format = formatted_data,
  modeling.id = modeling.id,
  models = available_models,
  CV.strategy = cv_strategy,
  OPT.strategy = "user.defined",
  OPT.user = bm.options,
  metric.eval = metric_eval,
  var.import = var_importance_rep,
  nb.cpu = nb_cpu,
  seed.val = seed_val,
  do.progress = TRUE
)

if (identical(cv_strategy, "user.defined")) {
  modeling_args$CV.user.table <- cv_table
} else {
  modeling_args$CV.nb.rep <- cv_nb_rep
  modeling_args$CV.perc <- cv_perc
  modeling_args$CV.k <- cv_k
}

myBiomodModel <- do.call(BIOMOD_Modeling, modeling_args)

cat("  Model training completed\n")

# ============================================================================
# 4. 提取评估结果
# ============================================================================

cat("\n  Extracting evaluation results...\n")

# 提取评估指标
model_eval <- get_evaluations(myBiomodModel)

# 提取变量重要性
model_varimp <- get_variables_importance(myBiomodModel)

# 方案 B 没有独立 evaluation 集，因此优先使用 validation 列
score_dataset <- if ("evaluation" %in% names(model_eval) &&
                     any(is.finite(model_eval$evaluation), na.rm = TRUE)) {
  "evaluation"
} else {
  "validation"
}

cat("  Aggregating model scores from dataset:", score_dataset, "\n")

# 计算各模型平均 TSS 和 AUC
eval_summary <- model_eval %>%
  group_by(algo, metric.eval) %>%
  summarise(
    mean = mean(.data[[score_dataset]], na.rm = TRUE),
    sd = sd(.data[[score_dataset]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(id_cols = "algo", names_from = "metric.eval", values_from = c("mean", "sd"))


cat("  Models evaluated:", length(available_models), "\n")

# 找出最佳模型
if ("mean_TSS" %in% names(eval_summary)) {
  best_tss_idx <- which.max(eval_summary$mean_TSS)
  best_tss_model <- eval_summary$algo[best_tss_idx]
  best_tss_value <- eval_summary$mean_TSS[best_tss_idx]
} else {
  best_tss_model <- NA
  best_tss_value <- NA
}

if ("mean_AUCroc" %in% names(eval_summary)) {
  best_auc_idx <- which.max(eval_summary$mean_AUCroc)
  best_auc_model <- eval_summary$algo[best_auc_idx]
  best_auc_value <- eval_summary$mean_AUCroc[best_auc_idx]
} else {
  best_auc_model <- NA
  best_auc_value <- NA
}

response_curve_algos <- available_models
if (length(response_curve_algos) == 0) {
  response_curve_algos <- unique(na.omit(c(best_tss_model, best_auc_model)))
}

all_built_models <- get_built_models(myBiomodModel)
response_curve_models <- unlist(lapply(response_curve_algos, function(algo) {
  grep(paste0("_allRun_", algo, "$"), all_built_models, value = TRUE)
}))
response_curve_models <- unique(response_curve_models)

# ============================================================================
# 5. 保存结果
# ============================================================================

cat("\n  Saving results...\n")

# 确保输出目录存在
dir.create(output_dirs$models, showWarnings = FALSE, recursive = TRUE)

# 保存完整模型结果
model_results_file <- file.path(output_dirs$models, "model_results.RData")
save(myBiomodModel, file = model_results_file)
cat("  Saved model results:", model_results_file, "\n")

# 保存评估指标
eval_file <- file.path(output_dirs$models, "model_evaluations.csv")
write.csv(model_eval, eval_file, row.names = FALSE)
cat("  Saved evaluations:", eval_file, "\n")

# 保存变量重要性
varimp_file <- file.path(output_dirs$models, "variable_importance.csv")
write.csv(model_varimp, varimp_file, row.names = FALSE)
cat("  Saved variable importance:", varimp_file, "\n")

# 保存评估摘要
eval_summary_file <- file.path(output_dirs$models, "model_evaluation_summary.csv")
write.csv(eval_summary, eval_summary_file, row.names = FALSE)
cat("  Saved evaluation summary:", eval_summary_file, "\n")

# ============================================================================
# 6. 可视化
# ============================================================================

cat("\n  Generating visualization...\n")

pdf_file <- file.path(output_dirs$models, "model_training_summary.pdf")
pdf(pdf_file, width = 14, height = 10)

# 第1页：模型平均评估指标 2D图 (TSS vs AUCroc)
invisible(bm_PlotEvalMean(
  bm.out = myBiomodModel,
  metric.eval = c("TSS", "AUCroc"),  # <--- 必须同时传 2 个指标
  dataset = score_dataset,
  group.by = "algo"
))

# 第2页：模型平均评估指标 2D图 (AUCprg vs BOYCE)
available_eval_metrics <- intersect(metric_eval, c("AUCprg", "BOYCE"))
if (length(available_eval_metrics) == 2) {
  invisible(bm_PlotEvalMean(
    bm.out = myBiomodModel,
    metric.eval = available_eval_metrics,
    dataset = score_dataset,
    group.by = "algo"
  ))
}

# 第3页：模型评估箱线图
invisible(bm_PlotEvalBoxplot(
  bm.out = myBiomodModel,
  dataset = score_dataset,
  group.by = c("algo", "run")
))

# 第4页：变量重要性
invisible(bm_PlotVarImpBoxplot(
  bm.out = myBiomodModel,
  group.by = c("expl.var", "run", "algo"),
  main = "Variable Importance"
))

# 第5页起：重点模型响应曲线
if (length(response_curve_models) > 0) {
  tryCatch({
    invisible(bm_PlotResponseCurves(
      bm.out = myBiomodModel,
      models.chosen = response_curve_models,
      show_variables = metadata$variables,
      fixed.var = "median",
      do.bivariate = FALSE,
      do.plot = TRUE,
      main = paste0(
        "Response Curves: ",
        paste(response_curve_algos, collapse = ", ")
      )
    ))
  }, error = function(e) {
    plot.new()
    text(
      0.5,
      0.5,
      paste("Response curves could not be generated:", conditionMessage(e)),
      cex = 0.9
    )
  })
}

dev.off()

cat("  Saved visualization:", pdf_file, "\n")

# ============================================================================
# 7. 保存元数据
# ============================================================================

metadata_output <- list(
  modeling_id = modeling.id,
  species_name = species_name,
  models_trained = available_models,
  n_models = length(available_models),
  cv_strategy = cv_strategy,
  cv_nb_rep = cv_nb_rep,
  cv_k = cv_k,
  cv_perc = cv_perc,
  score_dataset = score_dataset,
  metric_eval = metric_eval,
  n_variables = metadata$n_variables,
  variables = metadata$variables,
  total_presence = metadata$n_presence,
  total_pseudo = metadata$n_pseudo,
  seed = seed_val,
  best_model_tss = best_tss_model,
  best_tss_value = best_tss_value,
  best_model_auc = best_auc_model,
  best_auc_value = best_auc_value,
  response_curve_algos = response_curve_algos,
  response_curve_models = response_curve_models
)

metadata_file <- file.path(output_dirs$models, "model_training_metadata.rds")
saveRDS(metadata_output, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  05 Model Training - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Summary:\n")
cat("  - Models trained:", paste(available_models, collapse = ", "), "\n")
if (!is.na(best_tss_model)) {
  cat("  - Best model (TSS):", best_tss_model, "(", round(best_tss_value, 3), ")\n")
}
if (!is.na(best_auc_model)) {
  cat("  - Best model (AUC):", best_auc_model, "(", round(best_auc_value, 3), ")\n")
}
cat("\n")
# 清理内存
gc()
