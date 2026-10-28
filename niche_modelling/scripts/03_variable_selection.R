# ============================================================================
# 03_variable_selection.R - VIF variable selection
# ============================================================================
# Functionality:
#   - Select climate variables via Pearson correlation + VIF
#   - Drop highly collinear variables to reduce multicollinearity
# ============================================================================

# Load config
source("config.R")

cat("\n")
cat("=================================================================\n")
cat("  03 Variable Selection - Starting\n")
cat("=================================================================\n")

# ============================================================================
# 1. Read data
# ============================================================================

# Read presence SpatVector
presence_file <- file.path(output_dirs$data_prep, "sargentodoxa_presence_spat.rds")
if (!file.exists(presence_file)) {
  stop("Error: Presence data not found. Run 01_data_preparation.R first.")
}

presence_spat <- readRDS(presence_file)
n_presence <- nrow(presence_spat)

cat("  Loading presence data:", presence_file, "\n")
cat("  Presence points:", n_presence, "\n")

# Read modern climate data
cat("\n  Loading modern climate data from:", modern_clim_dir, "\n")
clim_files <- list.files(modern_clim_dir, pattern = modern_clim_pattern, full.names = TRUE)

if (length(clim_files) == 0) {
  stop("Error: No climate files found in ", modern_clim_dir)
}

cat("  Found", length(clim_files), "climate files\n")

# Read climate rasters and order them
modern_clim <- rast(clim_files)
modern_clim <- modern_clim[[mixedsort(names(modern_clim))]]

names(modern_clim) <- paste0("Bio", 1:nlyr(modern_clim)) 

cat("  Climate variables:", length(modern_clim), "\n")
cat("  Variable names:", paste(names(modern_clim), collapse = ", "), "\n")

# ============================================================================
# 2. Extract climate values at presence points
# ============================================================================

cat("\n  Extracting climate values at presence points...\n")

# Extract values with terra::extract
clim_values <- terra::extract(modern_clim, presence_spat)

# Fix: convert to df directly and drop the first ID column
clim_df <- as.data.frame(clim_values)

if ("ID" %in% colnames(clim_df)) {
  clim_df <- clim_df[, !colnames(clim_df) %in% "ID"]
}

colnames(clim_df) <- names(modern_clim)

cat("  Extracted values matrix:", nrow(clim_df), "x", ncol(clim_df), "\n")

# Check for NA values
na_count <- sum(is.na(clim_df))
if (na_count > 0) {
  cat("  Warning: Found", na_count, "NA values in extracted data\n")
  cat("  Removing rows with NA values...\n")
  valid_rows <- complete.cases(clim_df)
  clim_df <- clim_df[valid_rows, ]
  cat("  Remaining rows:", nrow(clim_df), "\n")
}

# ============================================================================
# 3. Correlation analysis
# ============================================================================

cat("\n  Computing Pearson correlation matrix...\n")

# Compute the correlation matrix
corr_matrix <- cor(clim_df, method = "pearson")

# Identify highly correlated pairs (|r| > corr_threshold); upper triangle only to avoid A-B / B-A duplication
high_corr <- which(abs(corr_matrix) > corr_threshold & upper.tri(corr_matrix), arr.ind = TRUE)

# Build a data frame of highly correlated pairs
high_corr_pairs <- data.frame(
  var1 = rownames(corr_matrix)[high_corr[, 1]],
  var2 = colnames(corr_matrix)[high_corr[, 2]],
  correlation = round(corr_matrix[high_corr], 4)
)

if (nrow(high_corr_pairs) > 0) {
  high_corr_pairs$abs_correlation <- abs(high_corr_pairs$correlation)
  high_corr_pairs <- high_corr_pairs[order(high_corr_pairs$abs_correlation, decreasing = TRUE), ]
}

cat("  Found", nrow(high_corr_pairs), "highly correlated pairs (|r| > ", corr_threshold, ")\n")

# ============================================================================
# 4. Manual variable selection (based on ecological judgement)
# ============================================================================

cat("\n  Manual variable selection (based on ecological knowledge)...\n")

# Preset variables to keep (reference: main.r variable choice)
# Full Bio1-19 list
all_vars <- names(modern_clim)

if (!is.null(manual_var_select)) {
  selected_vars <- intersect(manual_var_select, all_vars)
} else {
  selected_vars <- all_vars
}

# Ensure preset variables exist in the climate data
selected_vars <- intersect(selected_vars, all_vars)

cat("  Selected variables (manual):", length(selected_vars), "\n")
cat("  ", paste(selected_vars, collapse = ", "), "\n")

# Update climate data (keep selected variables only)
clim_selected <- modern_clim[[selected_vars]]
clim_df_selected <- clim_df[, selected_vars]

# ============================================================================
# 5. VIF (variance inflation factor) analysis
# ============================================================================

if (use_vif_step) {
  cat("\n  Computing VIF (Variance Inflation Factor)...\n")
  cat("  VIF threshold:", vif_threshold, "\n")
  
  # Fix 1: ensure the data frame has only numeric columns and drop the ID column
  # The first column of terra::extract output is usually an ID
  if ("ID" %in% colnames(clim_df_selected)) {
    clim_df_selected <- clim_df_selected[, !colnames(clim_df_selected) %in% "ID"]
  }
  
  # Check for constant columns (zero variance); VIF cannot be computed for them
  variances <- apply(clim_df_selected, 2, var, na.rm = TRUE)
  if (any(variances == 0)) {
    const_vars <- names(variances[variances == 0])
    cat("  Warning: Removing constant variables:", paste(const_vars, collapse = ", "), "\n")
    clim_df_selected <- clim_df_selected[, variances > 0]
  }
  
  # Record initial VIF of all candidates to gauge collinearity before stepwise VIF
  initial_vif_results <- usdm::vif(clim_df_selected)
  initial_vif_df <- as.data.frame(initial_vif_results)
  colnames(initial_vif_df) <- c("variable", "initial_vif")
  initial_vif_df$variable <- trimws(as.character(initial_vif_df$variable))
  initial_vif_df$initial_vif <- round(initial_vif_df$initial_vif, 4)
  
  # Compute VIF with usdm::vifstep
  # train = TRUE means a data frame is passed in
  vif_results <- usdm::vifstep(clim_df_selected, th = vif_threshold)
  
  # Fix 2: correctly retrieve the final selected variable names
  # In the usdm VIF object, @variables holds the retained variable names
  final_vars <- trimws(as.character(vif_results@results$Variables))
  
  if (length(final_vars) == 0) {
    stop("Error: VIF filtering removed all variables. Check your input data for extreme collinearity.")
  }
  
  cat("  VIF step completed\n")
  cat("  Remaining variables after VIF filtering:", length(final_vars), "\n")
  cat("  ", paste(final_vars, collapse = ", "), "\n")
  
  # Fix 3: correctly extract VIF results for display
  # vif_results@results is a data.frame with Variables and VIF columns
  vif_df <- vif_results@results
  colnames(vif_df) <- c("variable", "final_vif")
  vif_df$variable <- trimws(as.character(vif_df$variable))
  vif_df$final_vif <- round(vif_df$final_vif, 4)
  
} else {
  cat("\n  Skipping VIF analysis (use_vif_step = FALSE)\n")
  final_vars <- selected_vars
  vif_df <- NULL
  initial_vif_df <- NULL
}

# Update the final climate data
clim_final <- modern_clim[[final_vars]]

# Build an advice table (does not change selection): list each variable collinearity risk and keep/drop status
bio_labels <- c(
  Bio1 = "Annual mean temperature",
  Bio2 = "Mean diurnal range",
  Bio3 = "Isothermality",
  Bio4 = "Temperature seasonality",
  Bio5 = "Max temperature of warmest month",
  Bio6 = "Min temperature of coldest month",
  Bio7 = "Temperature annual range",
  Bio8 = "Mean temperature of wettest quarter",
  Bio9 = "Mean temperature of driest quarter",
  Bio10 = "Mean temperature of warmest quarter",
  Bio11 = "Mean temperature of coldest quarter",
  Bio12 = "Annual precipitation",
  Bio13 = "Precipitation of wettest month",
  Bio14 = "Precipitation of driest month",
  Bio15 = "Precipitation seasonality",
  Bio16 = "Precipitation of wettest quarter",
  Bio17 = "Precipitation of driest quarter",
  Bio18 = "Precipitation of warmest quarter",
  Bio19 = "Precipitation of coldest quarter"
)

get_high_corr_partners <- function(var_name, pairs_df) {
  if (nrow(pairs_df) == 0) {
    return("")
  }
  pair_rows <- pairs_df[pairs_df$var1 == var_name | pairs_df$var2 == var_name, ]
  if (nrow(pair_rows) == 0) {
    return("")
  }
  partners <- ifelse(pair_rows$var1 == var_name, pair_rows$var2, pair_rows$var1)
  paste(
    paste0(partners, " (r=", sprintf("%.2f", pair_rows$correlation), ")"),
    collapse = "; "
  )
}

max_abs_corr <- sapply(all_vars, function(v) {
  vals <- abs(corr_matrix[v, setdiff(all_vars, v)])
  if (length(vals) == 0) {
    return(NA_real_)
  }
  max(vals, na.rm = TRUE)
})

n_high_corr_links <- sapply(all_vars, function(v) {
  if (nrow(high_corr_pairs) == 0) {
    return(0)
  }
  sum(high_corr_pairs$var1 == v | high_corr_pairs$var2 == v)
})

final_corr_matrix <- corr_matrix[final_vars, final_vars, drop = FALSE]
final_high_corr_pairs <- high_corr_pairs[
  high_corr_pairs$var1 %in% final_vars & high_corr_pairs$var2 %in% final_vars,
]

final_max_abs_corr <- sapply(all_vars, function(v) {
  if (!(v %in% final_vars) || length(final_vars) < 2) {
    return(NA_real_)
  }
  vals <- abs(corr_matrix[v, setdiff(final_vars, v)])
  if (length(vals) == 0) {
    return(NA_real_)
  }
  max(vals, na.rm = TRUE)
})

final_n_high_corr_links <- sapply(all_vars, function(v) {
  if (!(v %in% final_vars) || nrow(final_high_corr_pairs) == 0) {
    return(0)
  }
  sum(final_high_corr_pairs$var1 == v | final_high_corr_pairs$var2 == v)
})

advice_df <- data.frame(
  variable = all_vars,
  description = unname(bio_labels[all_vars]),
  in_manual_candidates = all_vars %in% selected_vars,
  selected_final = all_vars %in% final_vars,
  max_abs_corr = round(max_abs_corr, 4),
  n_high_corr_links = as.integer(n_high_corr_links),
  final_max_abs_corr = round(final_max_abs_corr, 4),
  final_n_high_corr_links = as.integer(final_n_high_corr_links),
  initial_vif = NA_real_,
  final_vif = NA_real_,
  vif_status = ifelse(all_vars %in% final_vars, "selected", "removed_or_not_tested"),
  high_corr_partners = sapply(all_vars, get_high_corr_partners, pairs_df = high_corr_pairs),
  suggestion = "",
  stringsAsFactors = FALSE
)

if (!is.null(initial_vif_df)) {
  advice_df$initial_vif[match(initial_vif_df$variable, advice_df$variable)] <- initial_vif_df$initial_vif
}

if (!is.null(vif_df)) {
  advice_df$final_vif[match(vif_df$variable, advice_df$variable)] <- vif_df$final_vif
  advice_df$vif_status[advice_df$variable %in% vif_df$variable] <- "selected_after_vif"
  advice_df$vif_status[advice_df$in_manual_candidates & !(advice_df$variable %in% vif_df$variable)] <- "removed_by_vif"
}

advice_df$suggestion <- ifelse(
  advice_df$selected_final & advice_df$final_n_high_corr_links == 0 &
    (is.na(advice_df$final_vif) | advice_df$final_vif <= vif_threshold),
  "Keep: low residual collinearity after VIF.",
  ifelse(
    advice_df$selected_final & advice_df$final_n_high_corr_links > 0,
    "Keep with caution: final set still has high-correlation partners.",
    ifelse(
      !advice_df$in_manual_candidates,
      "Not in manual candidate list: only reconsider if ecologically important.",
      ifelse(
        !advice_df$selected_final & !is.na(advice_df$initial_vif) & advice_df$initial_vif > vif_threshold,
        "Removed by VIF: high initial VIF among candidate variables.",
        "Removed by VIF/manual step: reconsider only with a clear ecological reason."
      )
    )
  )
)

# ============================================================================
# 6. Save results
# ============================================================================

cat("\n  Saving results...\n")

# Ensure output directory exists
dir.create(output_dirs$variable_select, showWarnings = FALSE, recursive = TRUE)

# Save the correlation matrix
corr_file <- file.path(output_dirs$variable_select, "correlation_matrix.csv")
write.csv(corr_matrix, corr_file)
cat("  Saved correlation matrix:", corr_file, "\n")

# Save highly correlated pairs
high_corr_file <- file.path(output_dirs$variable_select, "high_correlation_pairs.csv")
write.csv(high_corr_pairs, high_corr_file, row.names = FALSE)
cat("  Saved high correlation pairs:", high_corr_file, "\n")

# Save the variable-selection advice
advice_file <- file.path(output_dirs$variable_select, "variable_selection_advice.csv")
write.csv(advice_df, advice_file, row.names = FALSE)
cat("  Saved variable selection advice:", advice_file, "\n")

# Save variables remaining after correlation filtering
corr_vars_file <- file.path(output_dirs$variable_select, "variables_after_corr.csv")
write.csv(data.frame(variable = selected_vars), corr_vars_file, row.names = FALSE)
cat("  Saved variables after correlation:", corr_vars_file, "\n")

# Save VIF results
if (!is.null(initial_vif_df)) {
  initial_vif_file <- file.path(output_dirs$variable_select, "vif_initial_all_candidates.csv")
  write.csv(initial_vif_df, initial_vif_file, row.names = FALSE)
  cat("  Saved initial VIF for all candidates:", initial_vif_file, "\n")
}

if (!is.null(vif_df)) {
  vif_file <- file.path(output_dirs$variable_select, "vif_results.csv")
  write.csv(vif_df, vif_file, row.names = FALSE)
  cat("  Saved VIF results:", vif_file, "\n")
}

# Save the final selected variables
final_vars_file <- file.path(output_dirs$variable_select, "variables_selected.csv")
write.csv(data.frame(variable = final_vars), final_vars_file, row.names = FALSE)
cat("  Saved final selected variables:", final_vars_file, "\n")

# ============================================================================
# 7. Visualization
# ============================================================================

cat("\n  Generating visualization...\n")

pdf_file <- file.path(output_dirs$variable_select, "variable_selection.pdf")
if (isTRUE(capabilities("cairo"))) {
  grDevices::cairo_pdf(pdf_file, width = 14, height = 12, onefile = TRUE, family = "sans")
} else {
  warning("Cairo graphics is not available; Chinese text in PDF may not render correctly.")
  grDevices::pdf(pdf_file, width = 14, height = 12, onefile = TRUE)
}
pdf_device <- grDevices::dev.cur()
cat("  PDF device opened:", names(pdf_device), "(", as.integer(pdf_device), ")\n")

bio_order_key <- function(x) as.integer(sub("^Bio", "", x))

# Page 1: BioClim variable descriptions
plot.new()
par(mar = c(4, 4, 4, 2))
bio_desc_df <- data.frame(
  variable = names(bio_labels),
  description = unname(bio_labels),
  stringsAsFactors = FALSE
)
bio_desc_df <- bio_desc_df[order(bio_order_key(bio_desc_df$variable)), ]
text(0.02, 0.95, "BioClim Variable Reference", adj = 0, cex = 1.4, font = 2)
text(0.02, 0.89,
     "Use this page to interpret the correlation and VIF diagnostics.",
     adj = 0, cex = 0.9)
y_desc <- seq(0.82, 0.08, length.out = nrow(bio_desc_df))
for (i in seq_len(nrow(bio_desc_df))) {
  text(0.05, y_desc[i], bio_desc_df$variable[i], adj = 0, cex = 0.78, font = 2)
  text(0.18, y_desc[i], bio_desc_df$description[i], adj = 0, cex = 0.78)
}

# Page 2: correlation heatmap
corrplot(corr_matrix, type = "upper", order = "original",
         tl.cex = 0.8, tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         title = "Correlation Matrix (Pearson)",
         mar = c(0, 0, 2, 0))

# Page 3: highly correlated pairs
if (nrow(high_corr_pairs) > 0) {
  par(mar = c(7, 9, 4, 2))
  high_corr_plot <- high_corr_pairs[order(bio_order_key(high_corr_pairs$var1),
                                          bio_order_key(high_corr_pairs$var2)), ]
  # With horiz = TRUE the first item is at the bottom; reverse so left labels read Bio1-Bio19 top to bottom.
  high_corr_plot <- high_corr_plot[rev(seq_len(nrow(high_corr_plot))), ]
  pair_cols <- ifelse(high_corr_plot$correlation > 0, "firebrick3", "steelblue3")
  bar_mid <- barplot(high_corr_plot$correlation,
                     names.arg = paste(high_corr_plot$var1, "-", high_corr_plot$var2),
                     horiz = TRUE, las = 2, cex.names = 0.7,
                     main = paste("High Correlation Pairs (|r| > ", corr_threshold, ")"),
                     col = pair_cols, border = NA,
                     xlim = c(-1, 1), axes = FALSE,
                     xlab = "Pearson correlation (r)")
  axis(1, at = seq(-1, 1, by = 0.25), las = 1)
  abline(v = 0, col = "grey25", lwd = 1.3)
  abline(v = corr_threshold, col = "firebrick3", lty = 2, lwd = 1.2)
  abline(v = -corr_threshold, col = "steelblue3", lty = 2, lwd = 1.2)
  value_x <- ifelse(high_corr_plot$correlation >= 0,
                    pmin(high_corr_plot$correlation + 0.03, 0.98),
                    pmax(high_corr_plot$correlation - 0.03, -0.98))
  value_adj <- ifelse(high_corr_plot$correlation >= 0, 0, 1)
  text(value_x, bar_mid, sprintf("%.2f", high_corr_plot$correlation),
       adj = value_adj, cex = 0.62)
  n_positive_pairs <- sum(high_corr_plot$correlation > 0)
  n_negative_pairs <- sum(high_corr_plot$correlation < 0)
  legend("topright",
         legend = c(paste0("Positive: ", n_positive_pairs),
                    paste0("Negative: ", n_negative_pairs)),
         fill = c("firebrick3", "steelblue3"), cex = 0.8, bg = "white")
  if (n_negative_pairs == 0) {
    mtext("No negative high-correlation pairs exceed the threshold in this run.",
          side = 1, line = 4.5, col = "steelblue4", cex = 0.8)
  }
} else {
  plot.new()
  text(0.5, 0.5, "No high correlation pairs found", cex = 1.5)
}

# Page 4: VIF bar chart
if (!is.null(vif_df)) {
  par(mar = c(5, 8, 4, 2))
  vif_df_sorted <- vif_df[order(bio_order_key(vif_df$variable)), ]
  # With horiz = TRUE the first item is at the bottom; reverse so left labels read in Bio order top to bottom.
  vif_df_sorted <- vif_df_sorted[rev(seq_len(nrow(vif_df_sorted))), ]
  barplot(vif_df_sorted$final_vif, names.arg = vif_df_sorted$variable,
          horiz = TRUE, las = 2, cex.names = 0.8,
          main = paste("Final VIF Values After Stepwise Filtering (threshold =", vif_threshold, ")"),
          col = ifelse(vif_df_sorted$final_vif > vif_threshold, "red", "green"),
          xlim = c(0, max(vif_df_sorted$final_vif, vif_threshold, na.rm = TRUE) * 1.2),
          xlab = "Final VIF after vifstep")
  abline(v = vif_threshold, col = "red", lty = 2)
  legend("bottomright", legend = c("Above threshold", "Below threshold"),
         fill = c("red", "green"), cex = 0.8)
} else {
  plot.new()
  text(0.5, 0.5, "VIF analysis not performed", cex = 1.5)
}

# Page 5: variable-selection workflow summary
plot.new()
text(0.5, 0.9, "Variable Selection Summary", cex = 1.5, font = 2)
text(0.5, 0.75, paste("Initial variables:", length(all_vars)), cex = 1.2)
text(0.5, 0.65, paste("After correlation filtering (|r| >", corr_threshold, "):", length(selected_vars)), cex = 1.2)
if (use_vif_step) {
  text(0.5, 0.55, paste("After VIF filtering (threshold =", vif_threshold, "):", length(final_vars)), cex = 1.2)
} else {
  text(0.5, 0.55, paste("VIF filtering: skipped"), cex = 1.2)
}
text(0.5, 0.35, "Final selected variables:", cex = 1.1)
text(0.5, 0.25, paste(final_vars, collapse = ", "), cex = 0.8)

# Page 6: residual correlation check of final variables
if (length(final_vars) >= 2) {
  final_vars_plot <- final_vars[order(bio_order_key(final_vars))]
  final_corr_matrix_plot <- final_corr_matrix[final_vars_plot, final_vars_plot, drop = FALSE]
  corrplot(final_corr_matrix_plot, type = "upper", order = "original",
           tl.cex = 0.9, tl.col = "black",
           col = colorRampPalette(c("blue", "white", "red"))(100),
           title = "Residual Correlation Among Final Variables",
           mar = c(0, 0, 2, 0))
  if (nrow(final_high_corr_pairs) > 0) {
    mtext("Warning: final variables still contain high-correlation pairs.", side = 1, line = 1, col = "red")
  } else {
    mtext("No final variable pairs exceed the correlation threshold.", side = 1, line = 1, col = "darkgreen")
  }
} else {
  plot.new()
  text(0.5, 0.5, "Only one final variable selected; residual correlation check is not applicable.", cex = 1.2)
}

# Page 7: variable-selection advice
plot.new()
par(mar = c(4, 4, 4, 2))
text(0.02, 0.95, "Variable Selection Advice", adj = 0, cex = 1.4, font = 2)
text(0.02, 0.89,
     "All columns are calculated before final filtering; Final columns are recalculated among final selected variables only.",
     adj = 0, cex = 0.85)

advice_plot <- advice_df[order(bio_order_key(advice_df$variable)), ]

y_pos <- seq(0.77, 0.17, length.out = nrow(advice_plot))
headers <- c("Var", "Final", "All |r|", "All links", "Final links", "Init VIF", "Final VIF", "Action")
x_pos <- c(0.02, 0.10, 0.18, 0.28, 0.39, 0.51, 0.61, 0.72)
text(x_pos, rep(0.83, length(x_pos)), headers, adj = 0, cex = 0.68, font = 2)
abline(h = 0.805, col = "grey70")

for (i in seq_len(nrow(advice_plot))) {
  row_col <- ifelse(advice_plot$selected_final[i], "darkgreen", "grey30")
  text(x_pos[1], y_pos[i], advice_plot$variable[i], adj = 0, cex = 0.66, col = row_col)
  text(x_pos[2], y_pos[i], ifelse(advice_plot$selected_final[i], "Yes", "No"), adj = 0, cex = 0.66, col = row_col)
  text(x_pos[3], y_pos[i], sprintf("%.2f", advice_plot$max_abs_corr[i]), adj = 0, cex = 0.66, col = row_col)
  text(x_pos[4], y_pos[i], advice_plot$n_high_corr_links[i], adj = 0, cex = 0.66, col = row_col)
  text(x_pos[5], y_pos[i], advice_plot$final_n_high_corr_links[i], adj = 0, cex = 0.66, col = row_col)
  text(x_pos[6], y_pos[i], ifelse(is.na(advice_plot$initial_vif[i]), "-", sprintf("%.2f", advice_plot$initial_vif[i])),
       adj = 0, cex = 0.66, col = row_col)
  text(x_pos[7], y_pos[i], ifelse(is.na(advice_plot$final_vif[i]), "-", sprintf("%.2f", advice_plot$final_vif[i])),
       adj = 0, cex = 0.66, col = row_col)
  text(x_pos[8], y_pos[i], advice_plot$suggestion[i], adj = 0, cex = 0.43, col = row_col)
}

text(0.02, 0.05,
     paste("Full table saved to:", advice_file),
     adj = 0, cex = 0.75)

if (grDevices::dev.cur() == pdf_device) {
  grDevices::dev.off()
} else if (as.integer(pdf_device) %in% grDevices::dev.list()) {
  grDevices::dev.set(pdf_device)
  grDevices::dev.off()
} else {
  warning("PDF device was already closed before visualization finished.")
}

cat("  Saved visualization:", pdf_file, "\n")
cat("  Expected PDF pages: 7\n")

# ============================================================================
# 8. Save metadata
# ============================================================================

metadata <- list(
  n_presence = n_presence,
  initial_vars = all_vars,
  n_initial_vars = length(all_vars),
  correlation_threshold = corr_threshold,
  n_high_corr_pairs = nrow(high_corr_pairs),
  selected_vars_after_corr = selected_vars,
  n_vars_after_corr = length(selected_vars),
  use_vif_step = use_vif_step,
  vif_threshold = vif_threshold,
  final_vars = final_vars,
  n_final_vars = length(final_vars),
  n_final_high_corr_pairs = nrow(final_high_corr_pairs),
  initial_vif_values = if (!is.null(initial_vif_df)) setNames(initial_vif_df$initial_vif, initial_vif_df$variable) else NULL,
  final_vif_values = if (!is.null(vif_df)) setNames(vif_df$final_vif, vif_df$variable) else NULL
)

metadata_file <- file.path(output_dirs$variable_select, "variable_selection_metadata.rds")
saveRDS(metadata, metadata_file)
cat("  Saved metadata:", metadata_file, "\n")

cat("\n")
cat("=================================================================\n")
cat("  03 Variable Selection - Completed\n")
cat("=================================================================\n")
cat("\n")
cat("  Summary:\n")
cat("  - Initial variables:", length(all_vars), "\n")
cat("  - After correlation filtering:", length(selected_vars), "\n")
if (use_vif_step) {
  cat("  - After VIF filtering:", length(final_vars), "\n")
}
cat("  - Final selected variables:", paste(final_vars, collapse = ", "), "\n")
cat("\n")
# Free memory
gc()
