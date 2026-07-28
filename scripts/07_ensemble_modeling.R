# ============================================================================
# 07_ensemble_modeling.R - Ensemble Modeling
# ============================================================================
# Functionality:
#   - Combine multiple quality models into an ensemble
#   - Support multiple ensemble algorithms
#   - Evaluate ensemble performance
#   - Extract and visualize ensemble variable importance and response curves
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  07 Ensemble Modeling - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Load data
# ============================================================================

# Load model training results
model_results_file <- file.path(output_dirs$models, "model_results.RData")
if (!file.exists(model_results_file)) {
  stop("Error: Model results not found. Run 05_model_training.R first.")
}

load(model_results_file)
cat("  Loaded model results:", model_results_file, "\n")

# Load selected model list
selected_file <- file.path(output_dirs$evaluation, "selected_models.rds")
if (!file.exists(selected_file)) {
  stop("Error: Selected models not found. Run 06_model_evaluation.R first.")
}

selected_models <- readRDS(selected_file)
cat("  Loaded selected models:", paste(selected_models, collapse = ", "), "\n")
cat("  Number of selected models:", length(selected_models), "\n")

training_metadata_file <- file.path(output_dirs$models, "model_training_metadata.rds")
training_metadata <- NULL
if (file.exists(training_metadata_file)) {
  training_metadata <- readRDS(training_metadata_file)
}

# ============================================================================
# 2. Ensemble modeling
# ============================================================================

cat("\n  Building ensemble models...\n")
cat("  Ensemble algorithms:", paste(ensemble_algo, collapse = ", "), "\n")
cat("  Selection metric:", paste(metric_select, collapse = ", "), "\n")
cat("  Selection thresholds:", paste(metric_select_thresh, collapse = ", "), "\n")

# 1. Get full names of all built models from myBiomodModel
all_built_models <- get_built_models(myBiomodModel)

# 2. Match selected_models algo short-names to full model names
#    Append "$" to anchor the algo name at string end (distinguishes "RF" from "RFd")
models_to_use <- unlist(sapply(selected_models, function(algo) {
  # Build an end-anchored regex, e.g. "RF$"
  pattern <- paste0(algo, "$")
  grep(pattern, all_built_models, value = TRUE)
}))

# Remove names from the result (sapply may name the vector, which can cause issues)
names(models_to_use) <- NULL

# Check whether models were found
if (length(models_to_use) == 0) {
  stop("Error: No full model names were found matching the selected algorithms.")
}

cat("  Found", length(models_to_use), "full model runs to use for ensemble.\n")
cat("  Selected algorithms for ensemble:", paste(ensemble_algo, collapse = ", "), "\n")
cat("  Metric for model selection:", paste(metric_select, collapse = ", "), "\n")

# Ensure output directory exists
models_dir <- "./Models"
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

# Set modeling.id (same as during training)
modeling.id <- myBiomodModel@modeling.id

# Ensemble modeling
ensemble_model <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModel,
  models.chosen = models_to_use,
  em.by = em.by,
  em.algo = ensemble_algo,
  metric.select = metric_select,
  metric.select.thresh = metric_select_thresh,
  metric.eval = metric_eval,
  EMwmean.decay = emwmean_decay,
  var.import = var_importance_rep,
  nb.cpu = nb_cpu,
  seed.val = seed_val
)

cat("  Ensemble modeling completed\n")

# ============================================================================
# 3. Extract evaluation results
# ============================================================================

cat("\n Extracting ensemble evaluation results...\n")

# Extract ensemble evaluation metrics
ensemble_eval <- get_evaluations(ensemble_model)

pick_score_dataset <- function(eval_df, preferred = NULL) {
  candidate_order <- unique(c(
    preferred,
    "evaluation",
    "validation",
    "calibration"
  ))
  candidate_order <- candidate_order[!is.na(candidate_order) & nzchar(candidate_order)]
  
  for (dataset_name in candidate_order) {
    if (dataset_name %in% names(eval_df) && any(is.finite(eval_df[[dataset_name]]), na.rm = TRUE)) {
      return(dataset_name)
    }
  }
  
  stop("Error: No finite scores found in ensemble evaluation table.")
}

preferred_dataset <- if (!is.null(training_metadata) && !is.null(training_metadata$score_dataset)) {
  training_metadata$score_dataset
} else {
  NULL
}

score_dataset <- pick_score_dataset(ensemble_eval, preferred = preferred_dataset)

cat("  Aggregating ensemble scores from dataset:", score_dataset, "\n")

# Compute mean performance per ensemble algorithm
ensemble_summary <- ensemble_eval %>%
  group_by(algo, metric.eval) %>%
  summarise(
    mean = mean(.data[[score_dataset]], na.rm = TRUE),
    sd = sd(.data[[score_dataset]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(id_cols = "algo", names_from = "metric.eval", values_from = c("mean", "sd"))

cat("  Ensemble algorithms evaluated:", length(ensemble_algo), "\n")

# Show ensemble evaluation results
cat("\n  Ensemble evaluation summary:\n")
print(ensemble_summary)

# ============================================================================
# 4. Save results
# ============================================================================

cat("\n  Saving results...\n")

# Ensure output directory exists
dir.create(output_dirs$ensemble, showWarnings = FALSE, recursive = TRUE)

# Save ensemble model
ensemble_file <- file.path(output_dirs$ensemble, "ensemble_model.RData")
save(ensemble_model, file = ensemble_file)
cat("  Saved ensemble model:", ensemble_file, "\n")

# Save ensemble evaluation metrics
ensemble_eval_file <- file.path(output_dirs$ensemble, "ensemble_evaluations.csv")
write.csv(ensemble_eval, ensemble_eval_file, row.names = FALSE)
cat("  Saved ensemble evaluations:", ensemble_eval_file, "\n")

# Save ensemble evaluation summary
ensemble_summary_file <- file.path(output_dirs$ensemble, "ensemble_summary.csv")
write.csv(ensemble_summary, ensemble_summary_file, row.names = FALSE)
cat("  Saved ensemble summary:", ensemble_summary_file, "\n")

# ============================================================================
# 5. Visualization
# ============================================================================

cat("\n  Generating visualization...\n")

# --- Extract ensemble variable importance ---
ensemble_var_imp <- NULL
tryCatch({
  ensemble_var_imp <- get_variables_importance(ensemble_model)
  if (!is.null(ensemble_var_imp) && nrow(ensemble_var_imp) > 0) {
    ensemble_var_imp_file <- file.path(output_dirs$ensemble, "ensemble_variable_importance.csv")
    write.csv(ensemble_var_imp, ensemble_var_imp_file, row.names = FALSE)
    cat("  Saved ensemble variable importance:", ensemble_var_imp_file, "\n")
  }
}, error = function(e) {
  cat("  Warning: Could not extract ensemble variable importance:", conditionMessage(e), "\n")
})

# --- Model names needed for ensemble response curves ---
ensemble_response_models <- NULL
tryCatch({
  ensemble_built <- get_built_models(ensemble_model)
  if (length(ensemble_built) > 0) {
    ensemble_response_models <- ensemble_built
    cat("  Ensemble models available for response curves:", length(ensemble_response_models), "\n")
  }
}, error = function(e) {
  cat("  Warning: Could not get ensemble built models:", conditionMessage(e), "\n")
})

# Get variable list (from training metadata or raw data)
response_curve_variables <- if (!is.null(training_metadata) && !is.null(training_metadata$variables)) {
  training_metadata$variables
} else {
  tryCatch(names(get_formal_data(myBiomodModel, "expl.var")), error = function(e) NULL)
}

pdf_file <- file.path(output_dirs$ensemble, "ensemble_summary.pdf")
pdf(pdf_file, width = 12, height = 8)

# --- Page 1: ensemble performance bar chart (TSS & AUC) ---
# Data preparation
plot_data_summary <- ensemble_summary %>%
  tidyr::pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    names_prefix = "mean_",
    values_to = "mean_value"
  ) %>%
  dplyr::filter(metric %in% c("TSS", "AUCroc"))

p1 <- ggplot(plot_data_summary, aes(x = algo, y = mean_value, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(mean_value, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 4) +
  scale_fill_manual(values = c("TSS" = "skyblue", "AUCroc" = "orange")) +
  labs(title = "Ensemble Model Performance Comparison",
       x = "Ensemble Algorithm",
       y = "Mean Score") +
  theme_bw() +
  ylim(0, 1.1)

print(p1)

# --- Page 2: ensemble evaluation boxplots (TSS & AUC) ---
# Data preparation
plot_data_boxplot <- ensemble_eval %>%
  dplyr::filter(metric.eval %in% c("TSS", "AUCroc")) %>%
  dplyr::mutate(score_value = .data[[score_dataset]]) %>%
  dplyr::filter(is.finite(score_value))

if (nrow(plot_data_boxplot) > 0) {
  p2 <- ggplot(plot_data_boxplot, aes(x = algo, y = score_value, fill = algo)) +
    geom_boxplot() +
    facet_wrap(~metric.eval, scales = "free_y") +
    labs(title = "Ensemble Model Evaluation Boxplots",
         x = "Ensemble Algorithm",
         y = paste0(score_dataset, " score")) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(p2)
} else {
  plot.new()
  text(0.5, 0.5, paste("No finite", score_dataset, "scores available for ensemble boxplots"), cex = 1.1)
}


# --- Page 3: single-model vs ensemble comparison (boxplots) ---
# (kept from the previous working version)
# Get single-model evaluation results
single_model_eval <- get_evaluations(myBiomodModel)
single_model_eval <- single_model_eval[single_model_eval$algo %in% selected_models, ]
single_score_dataset <- pick_score_dataset(single_model_eval, preferred = score_dataset)
single_score_column <- single_model_eval[[single_score_dataset]]
ensemble_score_column <- ensemble_eval[[score_dataset]]

# Prepare plotting data
plot_data_tss <- rbind(
  data.frame(Type = "Selected Single Models", Score = single_score_column[single_model_eval$metric.eval == "TSS"]),
  data.frame(Type = "Ensemble Models", Score = ensemble_score_column[ensemble_eval$metric.eval == "TSS"])
)

plot_data_auc <- rbind(
  data.frame(Type = "Selected Single Models", Score = single_score_column[single_model_eval$metric.eval == "AUCroc"]),
  data.frame(Type = "Ensemble Models", Score = ensemble_score_column[ensemble_eval$metric.eval == "AUCroc"])
)

plot_data_tss <- plot_data_tss[is.finite(plot_data_tss$Score), , drop = FALSE]
plot_data_auc <- plot_data_auc[is.finite(plot_data_auc$Score), , drop = FALSE]

# Set plot layout to 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Draw TSS comparison boxplot
if (nrow(plot_data_tss) > 0) {
  boxplot(Score ~ Type, data = plot_data_tss,
          main = "TSS: Single vs Ensemble",
          ylab = "TSS Score",
          col = c("orange", "steelblue"),
          ylim = c(0, 1))
} else {
  plot.new()
  text(0.5, 0.5, "No finite TSS scores available", cex = 1.1)
}

# Draw AUC comparison boxplot
if (nrow(plot_data_auc) > 0) {
  boxplot(Score ~ Type, data = plot_data_auc,
          main = "AUC: Single vs Ensemble",
          ylab = "AUC Score",
          col = c("orange", "steelblue"),
          ylim = c(0, 1))
} else {
  plot.new()
  text(0.5, 0.5, "No finite AUC scores available", cex = 1.1)
}

# Restore default plot settings
par(mfrow = c(1, 1))

# --- Page 4: ensemble variable importance ---
if (!is.null(ensemble_var_imp) && nrow(ensemble_var_imp) > 0) {
  tryCatch({
    invisible(bm_PlotVarImpBoxplot(
      bm.out = ensemble_model,
      group.by = c("expl.var", "algo"),
      main = "Ensemble Variable Importance"
    ))
  }, error = function(e) {
    # If bm_PlotVarImpBoxplot fails, draw manually with ggplot2
    tryCatch({
      var_imp_plot <- ensemble_var_imp %>%
        dplyr::group_by(expl.var, algo) %>%
        dplyr::summarise(
          mean_imp = mean(var.imp, na.rm = TRUE),
          sd_imp = sd(var.imp, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(mean_imp))

      p_varimp <- ggplot(var_imp_plot,
                          aes(x = reorder(expl.var, mean_imp),
                              y = mean_imp, fill = algo)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = mean_imp - sd_imp,
                          ymax = mean_imp + sd_imp),
                      position = position_dodge(width = 0.9),
                      width = 0.25) +
        coord_flip() +
        labs(title = "Ensemble Variable Importance",
             x = "Variable",
             y = "Importance") +
        theme_bw() +
        theme(legend.position = "bottom")
      print(p_varimp)
    }, error = function(e2) {
      plot.new()
      text(0.5, 0.5,
           paste("Ensemble variable importance plot failed:",
                 conditionMessage(e2)),
           cex = 0.9)
    })
  })
} else {
  plot.new()
  text(0.5, 0.5,
       "Ensemble variable importance not available",
       cex = 1.1)
}

# --- Page 5: ensemble response curves ---
if (!is.null(ensemble_response_models) &&
    length(ensemble_response_models) > 0 &&
    !is.null(response_curve_variables)) {
  tryCatch({
    invisible(bm_PlotResponseCurves(
      bm.out = ensemble_model,
      models.chosen = ensemble_response_models,
      show_variables = response_curve_variables,
      fixed.var = "median",
      do.bivariate = FALSE,
      do.plot = TRUE,
      main = "Ensemble Model Response Curves"
    ))
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5,
         paste("Ensemble response curves could not be generated:",
               conditionMessage(e)),
         cex = 0.9)
  })
} else {
  plot.new()
  text(0.5, 0.5,
       "Ensemble response curves not available",
       cex = 1.1)
}

dev.off()

cat("  Saved visualization:", pdf_file, "\n")

# Ensemble response curves in a separate PDF (avoids device conflict with summary PDF)
if (!is.null(ensemble_response_models) &&
    length(ensemble_response_models) > 0 &&
    !is.null(response_curve_variables)) {
  response_pdf <- file.path(output_dirs$ensemble,
                            "ensemble_response_curves.pdf")
  tryCatch({
    pdf(response_pdf, width = 14, height = 10)
    invisible(bm_PlotResponseCurves(
      bm.out = ensemble_model,
      models.chosen = ensemble_response_models,
      show_variables = response_curve_variables,
      fixed.var = "median",
      do.bivariate = FALSE,
      do.plot = TRUE,
      main = "Ensemble Model Response Curves"
    ))
    dev.off()
    cat("  Saved ensemble response curves:", response_pdf, "\n")
  }, error = function(e) {
    tryCatch(dev.off(), error = function(x) NULL)
    cat("  Warning: Ensemble response curves PDF failed:",
        conditionMessage(e), "\n")
  })
}

# ============================================================================
# 6. Save metadata
# ============================================================================

metadata <- list(
  modeling_id = modeling.id,
  selected_models = selected_models,
  n_selected = length(selected_models),
  ensemble_algo = ensemble_algo,
  metric_select = metric_select,
  metric_select_thresh = metric_select_thresh,
  score_dataset = score_dataset,
  emwmean_decay = emwmean_decay,
  ensemble_summary = ensemble_summary,
  best_ensemble = ensemble_summary$algo[which.max(ensemble_summary$mean_TSS)]
)

metadata_file <- file.path(output_dirs$ensemble, "ensemble_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  07 Ensemble Modeling - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Summary:\n")
cat("  - Single models:", length(selected_models), "\n")
cat("  - Ensemble algorithms:", paste(ensemble_algo, collapse = ", "), "\n")
cat("  - Best ensemble:", ensemble_summary$algo[which.max(ensemble_summary$mean_TSS)], 
    "(TSS =", round(max(ensemble_summary$mean_TSS, na.rm = TRUE), 4), ")\n")
cat("\n")
# Free memory
gc()
