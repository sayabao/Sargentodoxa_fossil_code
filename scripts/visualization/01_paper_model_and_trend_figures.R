# Paper figures: model diagnostics and post-fossil-start trend
#
# Outputs:
#   M2 variable correlation and final VIF
#   M3 model evaluation metrics
#   M4 ensemble variable importance
#   R3 suitability area trend from the earliest fossil-matched period onward
#   R3b latitudinal suitability-area-by-class TABLE (input to R3c)
#   (the >= 0.3 latitudinal HEATMAP is drawn by R3c_latitudinal_area_heatmap.R)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(patchwork)
library(terra)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "visualization") {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
} else if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

# Output root: defaults to "output" (frozen shipped outputs). Set the
# SARGENTO_OUTPUT_ROOT env var to regenerate figures from an alternative
# output tree (e.g. a fresh server copy under output_new/) without touching
# the frozen outputs. Only the "output" prefix is redirected; "data" stays.
output_root <- Sys.getenv("SARGENTO_OUTPUT_ROOT", "output")
path_here <- function(...) {
  parts <- list(...)
  if (length(parts) >= 1L && identical(parts[[1]], "output")) parts[[1]] <- output_root
  do.call(file.path, c(list(project_root), parts))
}

# Paleoclimate projection dir: current pipeline writes "09_paleoclimate";
# older frozen outputs used "09_paleoclimate_teyf". Accept either.
paleo_proj_dir <- function() {
  for (nm in c("09_paleoclimate_teyf", "09_paleoclimate")) {
    d <- path_here("output", nm, "projections")
    if (dir.exists(d)) return(d)
  }
  path_here("output", "09_paleoclimate_teyf", "projections")
}

dir.create(path_here("output", "figures_paper", "model_diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_here("output", "figures_paper", "trends"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_here("output", "figures_paper", "tables"), recursive = TRUE, showWarnings = FALSE)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      axis.title = element_text(color = "grey15"),
      axis.text = element_text(color = "grey20"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}

save_figure <- function(plot, stem, width, height, dpi = 600) {
  pdf_file <- paste0(stem, ".pdf")
  png_file <- paste0(stem, ".png")
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
}

metric_labels <- c(
  mean_TSS = "TSS",
  mean_AUCroc = "AUCroc",
  mean_AUCprg = "AUCprg",
  mean_BOYCE = "Boyce"
)

metric_sd_cols <- c(
  mean_TSS = "sd_TSS",
  mean_AUCroc = "sd_AUCroc",
  mean_AUCprg = "sd_AUCprg",
  mean_BOYCE = "sd_BOYCE"
)

bio_labels <- c(
  Bio1 = "Bio1 Annual mean temp.",
  Bio2 = "Bio2 Mean diurnal range",
  Bio3 = "Bio3 Isothermality",
  Bio4 = "Bio4 Temp. seasonality",
  Bio5 = "Bio5 Max temp. warmest month",
  Bio6 = "Bio6 Min temp. coldest month",
  Bio7 = "Bio7 Temp. annual range",
  Bio8 = "Bio8 Mean temp. wettest quarter",
  Bio9 = "Bio9 Mean temp. driest quarter",
  Bio10 = "Bio10 Mean temp. warmest quarter",
  Bio11 = "Bio11 Mean temp. coldest quarter",
  Bio12 = "Bio12 Annual precipitation",
  Bio13 = "Bio13 Precip. wettest month",
  Bio14 = "Bio14 Precip. driest month",
  Bio15 = "Bio15 Precip. seasonality",
  Bio16 = "Bio16 Precip. wettest quarter",
  Bio17 = "Bio17 Precip. driest quarter",
  Bio18 = "Bio18 Precip. warmest quarter",
  Bio19 = "Bio19 Precip. coldest quarter"
)

label_bio <- function(x) {
  ifelse(x %in% names(bio_labels), bio_labels[x], x)
}

# -------------------------------------------------------------------------
# M2. Variable filtering and final VIF
# -------------------------------------------------------------------------

correlation_matrix <- read_csv(path_here("output", "03_variable_select", "correlation_matrix.csv"), show_col_types = FALSE)
selected_variables <- read_csv(path_here("output", "03_variable_select", "variables_selected.csv"), show_col_types = FALSE) %>%
  pull(variable)
vif_results <- read_csv(path_here("output", "03_variable_select", "vif_results.csv"), show_col_types = FALSE)

correlation_long <- correlation_matrix %>%
  rename(var1 = 1) %>%
  pivot_longer(
    cols = -var1,
    names_to = "var2",
    values_to = "correlation"
  ) %>%
  mutate(
    var1 = factor(var1, levels = rev(names(bio_labels))),
    var2 = factor(var2, levels = names(bio_labels)),
    selected_pair = as.character(var1) %in% selected_variables & as.character(var2) %in% selected_variables
  )

selected_diag <- tibble::tibble(
  var1 = factor(selected_variables, levels = rev(names(bio_labels))),
  var2 = factor(selected_variables, levels = names(bio_labels))
)

vif_plot_data <- vif_results %>%
  mutate(
    variable_label = variable,
    variable_label = factor(variable_label, levels = variable_label[order(final_vif)])
  ) %>%
  arrange(desc(final_vif))

write_csv(
  correlation_long,
  path_here("output", "figures_paper", "tables", "M2_variable_correlation_long.csv")
)
write_csv(
  vif_plot_data,
  path_here("output", "figures_paper", "tables", "M2_final_vif_values.csv")
)

p_m2_corr <- ggplot(correlation_long, aes(x = var2, y = var1, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.12) +
  geom_tile(
    data = selected_diag,
    aes(x = var2, y = var1),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.35
  ) +
  coord_fixed() +
  scale_fill_gradient2(
    low = "#3C65A6",
    mid = "white",
    high = "#B33A3A",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  labs(
    title = "Candidate variable correlation",
    subtitle = "Black diagonal boxes mark final selected variables",
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 7) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6.6),
    panel.grid = element_blank()
  )

p_m2_vif <- ggplot(vif_plot_data, aes(x = variable_label, y = final_vif)) +
  geom_col(fill = "#5B8E7D", color = "grey20", linewidth = 0.15, width = 0.72) +
  geom_hline(yintercept = 10, linewidth = 0.35, linetype = "dashed", color = "#B33A3A") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Final variable VIF",
    subtitle = "Dashed line marks VIF = 10 reference",
    x = NULL,
    y = "VIF"
  ) +
  theme_paper(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.2)
  )

p_m2 <- p_m2_corr + p_m2_vif +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A")

save_figure(
  p_m2,
  path_here("output", "figures_paper", "model_diagnostics", "M2_variable_correlation_vif"),
  width = 9,
  height = 5.2
)

# -------------------------------------------------------------------------
# M3. Model evaluation
# -------------------------------------------------------------------------

model_ranking <- read_csv(path_here("output", "06_evaluation", "model_ranking.csv"), show_col_types = FALSE)
selected_eval <- read_csv(path_here("output", "06_evaluation", "selected_models_evaluation.csv"), show_col_types = FALSE) %>%
  arrange(desc(mean_TSS)) %>%
  mutate(rank = row_number())
ensemble_summary <- read_csv(path_here("output", "07_ensemble", "ensemble_summary.csv"), show_col_types = FALSE)

selected_algorithms <- selected_eval$algo

model_mean_long <- model_ranking %>%
  select(algo, rank, all_of(names(metric_labels))) %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric_col",
    values_to = "mean_value"
  )

model_sd_long <- model_ranking %>%
  select(algo, rank, all_of(unname(metric_sd_cols))) %>%
  pivot_longer(
    cols = all_of(unname(metric_sd_cols)),
    names_to = "sd_col",
    values_to = "sd_value"
  ) %>%
  mutate(metric_col = names(metric_sd_cols)[match(sd_col, unname(metric_sd_cols))]) %>%
  select(algo, rank, metric_col, sd_value)

model_long <- model_mean_long %>%
  left_join(model_sd_long, by = c("algo", "rank", "metric_col")) %>%
  mutate(
    metric = factor(unname(metric_labels[metric_col]), levels = unname(metric_labels)),
    model_type = if_else(algo %in% selected_algorithms, "Selected top 5 models", "Not selected models"),
    model_type = factor(model_type, levels = c("Selected top 5 models", "Not selected models")),
    algo = factor(algo, levels = model_ranking$algo[order(model_ranking$rank)])
  )

ensemble_long <- ensemble_summary %>%
  select(algo, all_of(names(metric_labels))) %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric_col",
    values_to = "mean_value"
  ) %>%
  mutate(
    sd_value = NA_real_,
    metric = factor(unname(metric_labels[metric_col]), levels = unname(metric_labels)),
    model_type = factor("Ensemble", levels = c("Selected top 5 models", "Not selected models", "Ensemble")),
    algo = factor(algo, levels = c(levels(model_long$algo), ensemble_summary$algo))
  )

eval_plot_order <- c(ensemble_summary$algo, model_ranking$algo[order(model_ranking$rank)])

eval_plot_data <- bind_rows(model_long, ensemble_long) %>%
  mutate(algo_label = factor(as.character(algo), levels = rev(eval_plot_order)))

write_csv(eval_plot_data, path_here("output", "figures_paper", "tables", "M3_model_evaluation_metrics.csv"))
write_csv(selected_eval, path_here("output", "figures_paper", "tables", "M3_selected_top5_models.csv"))

p_m3 <- ggplot(eval_plot_data, aes(x = algo_label, y = mean_value)) +
  geom_col(aes(fill = model_type), width = 0.72, color = "grey20", linewidth = 0.15) +
  geom_errorbar(
    aes(ymin = pmax(mean_value - sd_value, 0), ymax = pmin(mean_value + sd_value, 1)),
    width = 0.22,
    linewidth = 0.25,
    na.rm = TRUE
  ) +
  facet_wrap(~ metric, nrow = 1) +
  coord_flip(ylim = c(0, 1)) +
  scale_fill_manual(values = c(
    "Selected top 5 models" = "#4C78A8",
    "Not selected models" = "grey70",
    "Ensemble" = "#D95F02"
  )) +
  labs(
    title = "Evaluation metrics of candidate models",
    subtitle = "Models ranked by mean TSS; grey bars indicate models not selected for the ensemble",
    x = NULL,
    y = "Metric value",
    fill = NULL
  ) +
  theme_paper()

save_figure(
  p_m3,
  path_here("output", "figures_paper", "model_diagnostics", "M3_model_evaluation_metrics"),
  width = 9,
  height = 4.2
)

# -------------------------------------------------------------------------
# M4. Ensemble variable importance
# -------------------------------------------------------------------------

ensemble_varimp <- read_csv(path_here("output", "07_ensemble", "ensemble_variable_importance.csv"), show_col_types = FALSE)

varimp_summary <- ensemble_varimp %>%
  group_by(expl.var) %>%
  summarise(
    mean_importance = mean(var.imp, na.rm = TRUE),
    sd_importance = sd(var.imp, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance)) %>%
  mutate(
    variable_label = expl.var,
    variable_label = factor(variable_label, levels = rev(variable_label))
  )

write_csv(varimp_summary, path_here("output", "figures_paper", "tables", "M4_ensemble_variable_importance.csv"))

p_m4 <- ggplot(varimp_summary, aes(x = variable_label, y = mean_importance)) +
  geom_col(fill = "#5B8E7D", color = "grey20", linewidth = 0.15, width = 0.72) +
  geom_errorbar(
    aes(ymin = pmax(mean_importance - sd_importance, 0), ymax = mean_importance + sd_importance),
    width = 0.22,
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Ensemble variable importance",
    subtitle = "Mean permutation importance across repeated ensemble importance runs",
    x = NULL,
    y = "Variable importance"
  ) +
  theme_paper()

save_figure(
  p_m4,
  path_here("output", "figures_paper", "model_diagnostics", "M4_variable_importance"),
  width = 6.8,
  height = 4.8
)

# -------------------------------------------------------------------------
# R3. Suitability area trend from earliest fossil period onward
# -------------------------------------------------------------------------

teyf_simulations <- read_csv(path_here("data", "teyf_simulations.csv"), show_col_types = FALSE)
fossil_records <- read_csv(path_here("data", "fossil_record", "FossilRecord-paleocor-v2.1.csv"), show_col_types = FALSE)

earliest_fossil_ma <- max(fossil_records$age, na.rm = TRUE)

trend_periods <- teyf_simulations %>%
  transmute(period = `Simulation Name`, Ma = as.numeric(Ma)) %>%
  filter(is.finite(Ma), Ma <= earliest_fossil_ma) %>%
  arrange(desc(Ma)) %>%
  mutate(period_order = row_number(), ma_label = as.character(Ma))

fossil_counts <- fossil_records %>%
  count(age, name = "fossil_count")

suitability_classes <- tibble::tibble(
  suitability_class = factor(
    c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable"),
    levels = c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable")
  ),
  lower = c(0, 0.1, 0.3, 0.7),
  upper = c(0.1, 0.3, 0.7, 1),
  class_label = c("0-0.1", "0.1-0.3", "0.3-0.7", "0.7-1")
)

class_suitability <- function(x) {
  cut(
    x,
    breaks = c(0, 0.1, 0.3, 0.7, 1),
    labels = levels(suitability_classes$suitability_class),
    include.lowest = TRUE,
    right = FALSE
  )
}

normalize_prediction_raster <- function(r) {
  xmax <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (is.finite(xmax) && xmax > 1) {
    r <- r / 1000
  }
  r
}

projection_dir <- paleo_proj_dir()
projection_files <- list.files(projection_dir, pattern = "\\.tif$", full.names = TRUE)

find_period_file <- function(period) {
  hit <- projection_files[grepl(paste0("^", period, "_"), basename(projection_files))]
  if (length(hit) == 0) {
    stop("No projection raster found for period: ", period)
  }
  hit[1]
}

calculate_area_by_class <- function(period, ma) {
  raster_file <- find_period_file(period)
  r <- normalize_prediction_raster(terra::rast(raster_file))
  area_r <- terra::cellSize(r, unit = "km")
  names(r) <- "suitability"
  names(area_r) <- "area_km2"
  
  values_df <- as.data.frame(c(r, area_r), na.rm = TRUE)
  values_df %>%
    mutate(
      period = period,
      Ma = ma,
      suitability_class = class_suitability(suitability)
    ) %>%
    filter(!is.na(suitability_class)) %>%
    group_by(period, Ma, suitability_class) %>%
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
}

calculate_latitudinal_area_by_class <- function(period, ma, lat_bin_width = 10) {
  raster_file <- find_period_file(period)
  r <- normalize_prediction_raster(terra::rast(raster_file))
  area_r <- terra::cellSize(r, unit = "km")
  lat_r <- terra::init(r, "y")
  names(r) <- "suitability"
  names(area_r) <- "area_km2"
  names(lat_r) <- "latitude"
  
  values_df <- as.data.frame(c(r, area_r, lat_r), na.rm = TRUE)
  values_df %>%
    mutate(
      period = period,
      Ma = ma,
      lat_bin = floor(latitude / lat_bin_width) * lat_bin_width,
      lat_bin_mid = lat_bin + lat_bin_width / 2,
      lat_band = paste0(lat_bin, " to ", lat_bin + lat_bin_width),
      suitability_class = class_suitability(suitability)
    ) %>%
    filter(!is.na(suitability_class)) %>%
    group_by(period, Ma, lat_bin, lat_bin_mid, lat_band, suitability_class) %>%
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
}

calculate_latitudinal_continuous_suitability <- function(period, ma, lat_bin_width = 10) {
  raster_file <- find_period_file(period)
  r <- normalize_prediction_raster(terra::rast(raster_file))
  area_r <- terra::cellSize(r, unit = "km")
  lat_r <- terra::init(r, "y")
  names(r) <- "suitability"
  names(area_r) <- "area_km2"
  names(lat_r) <- "latitude"
  
  as.data.frame(c(r, area_r, lat_r), na.rm = TRUE) %>%
    mutate(
      period = period,
      Ma = ma,
      lat_bin = floor(latitude / lat_bin_width) * lat_bin_width,
      lat_bin_mid = lat_bin + lat_bin_width / 2,
      lat_band = paste0(lat_bin, " to ", lat_bin + lat_bin_width),
      suitability_weighted_area_km2 = area_km2 * pmin(pmax(suitability, 0), 1)
    ) %>%
    group_by(period, Ma, lat_bin, lat_bin_mid, lat_band) %>%
    summarise(
      weighted_area_km2 = sum(suitability_weighted_area_km2, na.rm = TRUE),
      .groups = "drop"
    )
}

message("Calculating suitability area classes for ", nrow(trend_periods), " TEYF periods...")
area_by_class <- dplyr::bind_rows(
  Map(calculate_area_by_class, trend_periods$period, trend_periods$Ma)
) %>%
  left_join(trend_periods %>% select(period, Ma, period_order, ma_label), by = c("period", "Ma")) %>%
  right_join(
    tidyr::crossing(
      trend_periods %>% select(period, Ma, period_order, ma_label),
      suitability_class = suitability_classes$suitability_class
    ),
    by = c("period", "Ma", "period_order", "ma_label", "suitability_class")
  ) %>%
  mutate(area_km2 = replace_na(area_km2, 0)) %>%
  left_join(suitability_classes %>% select(suitability_class, class_label), by = "suitability_class") %>%
  arrange(desc(Ma), suitability_class)

area_wide <- area_by_class %>%
  select(period, Ma, suitability_class, area_km2) %>%
  tidyr::pivot_wider(names_from = suitability_class, values_from = area_km2)

trend_data <- trend_periods %>%
  left_join(fossil_counts, by = c("Ma" = "age")) %>%
  mutate(fossil_count = replace_na(fossil_count, 0L))

extract_fossil_suitability <- function(row_index) {
  fossil_row <- fossil_records[row_index, , drop = FALSE]
  raster_file <- find_period_file(fossil_row$simulation_name)
  r <- normalize_prediction_raster(terra::rast(raster_file))
  point <- terra::vect(
    data.frame(x = fossil_row$paleoLong, y = fossil_row$paleoLat),
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  suitability_value <- terra::extract(r, point)[1, 2]
  fossil_row %>%
    mutate(
      extracted_suitability = as.numeric(suitability_value),
      suitability_class = class_suitability(extracted_suitability)
    )
}

message("Extracting fossil-point suitability values...")
fossil_suitability <- bind_rows(lapply(seq_len(nrow(fossil_records)), extract_fossil_suitability))

write_csv(area_by_class, path_here("output", "figures_paper", "tables", "R3_suitability_area_by_class_from_earliest_fossil.csv"))
write_csv(area_wide, path_here("output", "figures_paper", "tables", "R3_suitability_area_by_class_wide_from_earliest_fossil.csv"))
write_csv(fossil_suitability, path_here("output", "figures_paper", "tables", "R3_fossil_point_suitability_values.csv"))

p_r3_area <- ggplot(area_by_class, aes(x = period_order, y = area_km2 / 1e6, fill = suitability_class)) +
  geom_area(position = "stack", alpha = 0.95, color = "white", linewidth = 0.15) +
  scale_x_continuous(
    breaks = trend_periods$period_order,
    labels = trend_periods$ma_label,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.04)),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_fill_manual(
    values = c(
      "Least suitable" = "#E8E8E8",
      "Slightly suitable" = "#BFD8A9",
      "Moderately suitable" = "#66A65B",
      "Highly suitable" = "#1F6B3A"
    ),
    labels = c(
      "Least suitable (0-0.1)",
      "Slightly suitable (0.1-0.3)",
      "Moderately suitable (0.3-0.7)",
      "Highly suitable (0.7-1)"
    )
  ) +
  labs(
    title = "Paleoclimate suitability area trend",
    subtitle = paste0("Equally spaced TEYF periods from the earliest fossil-matched period (", earliest_fossil_ma, " Ma) to present"),
    x = NULL,
    y = expression("Area ("*10^6~km^2*")"),
    fill = NULL
  ) +
  theme_paper() +
  theme(
    legend.position = "top",
    legend.key.width = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_r3 <- p_r3_area +
  labs(x = "Age (Ma)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(
  p_r3,
  path_here("output", "figures_paper", "trends", "R3_suitability_area_trend_from_earliest_fossil"),
  width = 7.4,
  height = 4.6
)

# -------------------------------------------------------------------------
# R3b. Latitudinal suitability-area gradient from earliest fossil period onward
# -------------------------------------------------------------------------

lat_bin_width <- 10

message("Calculating latitudinal suitability area gradients for ", nrow(trend_periods), " TEYF periods...")
lat_area_by_class <- dplyr::bind_rows(
  Map(
    function(period, ma) calculate_latitudinal_area_by_class(period, ma, lat_bin_width = lat_bin_width),
    trend_periods$period,
    trend_periods$Ma
  )
) %>%
  left_join(trend_periods %>% select(period, Ma, period_order, ma_label), by = c("period", "Ma")) %>%
  right_join(
    tidyr::crossing(
      trend_periods %>% select(period, Ma, period_order, ma_label),
      tibble::tibble(
        lat_bin = seq(-90, 80, by = lat_bin_width),
        lat_bin_mid = lat_bin + lat_bin_width / 2,
        lat_band = paste0(lat_bin, " to ", lat_bin + lat_bin_width)
      ),
      suitability_class = suitability_classes$suitability_class
    ),
    by = c("period", "Ma", "period_order", "ma_label", "lat_bin", "lat_bin_mid", "lat_band", "suitability_class")
  ) %>%
  mutate(area_km2 = replace_na(area_km2, 0)) %>%
  left_join(suitability_classes %>% select(suitability_class, class_label), by = "suitability_class") %>%
  arrange(desc(Ma), lat_bin, suitability_class)

write_csv(
  lat_area_by_class,
  path_here("output", "figures_paper", "tables", "R3b_latitudinal_suitability_area_by_class.csv")
)

# The manuscript uses the >= 0.3 (moderately + highly suitable) latitudinal
# heatmap, produced separately by R3c_latitudinal_area_heatmap.R from the table
# above. The former R3b (0.7-1 only) heatmap figure is intentionally not drawn.

message("All requested first-pass paper figures were generated.")
