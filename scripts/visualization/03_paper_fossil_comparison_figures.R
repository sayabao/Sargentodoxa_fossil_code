# Paper figures: fossil-model comparison
#
# Outputs:
#   R4 fossil records over matched paleoclimate projections
#   R5 fossil-point suitability extraction

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(terra)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "visualization") {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
} else if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

path_here <- function(...) file.path(project_root, ...)

dir.create(path_here("output", "figures_paper", "fossil_comparison"), recursive = TRUE, showWarnings = FALSE)
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

theme_map <- function(base_size = 8) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "grey35", margin = margin(b = 6)),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.text = element_text(size = base_size - 1),
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.35, "cm"),
      plot.margin = margin(5, 5, 5, 5)
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

suitability_levels <- c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable")

class_suitability <- function(x) {
  cut(
    x,
    breaks = c(0, 0.1, 0.3, 0.7, 1),
    labels = suitability_levels,
    include.lowest = TRUE,
    right = FALSE
  )
}

suitability_colors <- c(
  "Least suitable" = "#E8E8E8",
  "Slightly suitable" = "#BFD8A9",
  "Moderately suitable" = "#66A65B",
  "Highly suitable" = "#1F6B3A"
)

suitability_scale <- scale_fill_manual(
  values = suitability_colors,
  labels = c(
    "Least (0-0.1)",
    "Slight (0.1-0.3)",
    "Moderate (0.3-0.7)",
    "High (0.7-1)"
  ),
  name = NULL,
  guide = guide_legend(nrow = 2, byrow = TRUE),
  drop = FALSE,
  na.value = "transparent"
)

organ_shapes <- c("Seed" = 21, "Pollen" = 24)

normalize_prediction_raster <- function(r) {
  xmax <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (is.finite(xmax) && xmax > 1) {
    r <- r / 1000
  }
  r
}

find_teyf_file <- function(period) {
  projection_dir <- path_here("output", "09_paleoclimate_teyf", "projections")
  hit <- list.files(
    projection_dir,
    pattern = paste0("^", period, "_.*\\.tif$"),
    full.names = TRUE
  )
  if (length(hit) == 0) {
    stop("No TEYF projection raster found for period: ", period)
  }
  hit[1]
}

raster_to_plot_df <- function(raster_file, aggregate_factor = 5, plot_extent = terra::ext(-130, 150, 15, 60)) {
  r <- normalize_prediction_raster(terra::rast(raster_file))
  r <- terra::crop(r, plot_extent)
  r <- terra::aggregate(r, fact = aggregate_factor, fun = "mean", na.rm = TRUE)
  names(r) <- "suitability"
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    mutate(
      suitability = pmin(pmax(suitability, 0), 1),
      suitability_class = factor(class_suitability(suitability), levels = suitability_levels)
    )
}

extract_fossil_suitability <- function(fossil_row) {
  raster_file <- find_teyf_file(fossil_row$simulation_name)
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
      suitability_class = factor(class_suitability(extracted_suitability), levels = suitability_levels)
    )
}

teyf_simulations <- read_csv(path_here("data", "teyf_simulations.csv"), show_col_types = FALSE)
fossil_records <- read_csv(path_here("data", "fossil_record", "FossilRecord-paleocor-v2.1.csv"), show_col_types = FALSE) %>%
  mutate(
    number = as.integer(number),
    organ = factor(organ),
    simulation_name = as.character(simulation_name)
  )

fossil_suitability <- bind_rows(lapply(seq_len(nrow(fossil_records)), function(i) {
  extract_fossil_suitability(fossil_records[i, , drop = FALSE])
})) %>%
  left_join(
    teyf_simulations %>%
      transmute(
        simulation_name = `Simulation Name`,
        Ma = as.numeric(Ma)
      ),
    by = "simulation_name"
  ) %>%
  mutate(
    Ma = ifelse(is.na(Ma), age, Ma),
    panel_label = paste0(Ma, " Ma"),
    fossil_label = as.character(number),
    age_label = paste0(Age_stage, " / ", Age, " Ma")
  ) %>%
  arrange(desc(Ma), number)

period_order <- fossil_suitability %>%
  distinct(simulation_name, Ma, panel_label) %>%
  arrange(desc(Ma))

fossil_suitability <- fossil_suitability %>%
  mutate(panel_label = factor(panel_label, levels = period_order$panel_label))

write_csv(
  fossil_suitability,
  path_here("output", "figures_paper", "tables", "R5_fossil_point_suitability.csv")
)

fossil_summary <- fossil_suitability %>%
  group_by(simulation_name, Ma, panel_label, suitability_class) %>%
  summarise(
    fossil_count = dplyr::n(),
    mean_suitability = mean(extracted_suitability, na.rm = TRUE),
    min_suitability = min(extracted_suitability, na.rm = TRUE),
    max_suitability = max(extracted_suitability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Ma), suitability_class)

write_csv(
  fossil_summary,
  path_here("output", "figures_paper", "tables", "R5_fossil_point_suitability_summary.csv")
)

# -------------------------------------------------------------------------
# R4. Fossil records over matched paleoclimate projections
# -------------------------------------------------------------------------

panel_map_data <- bind_rows(lapply(seq_len(nrow(period_order)), function(i) {
  raster_to_plot_df(find_teyf_file(period_order$simulation_name[i])) %>%
    mutate(
      simulation_name = period_order$simulation_name[i],
      Ma = period_order$Ma[i],
      panel_label = period_order$panel_label[i]
    )
}))

panel_map_data <- panel_map_data %>%
  mutate(panel_label = factor(panel_label, levels = period_order$panel_label))

p_r4 <- ggplot(panel_map_data, aes(x = x, y = y, fill = suitability_class)) +
  geom_raster() +
  geom_point(
    data = fossil_suitability,
    aes(x = paleoLong, y = paleoLat, shape = organ),
    inherit.aes = FALSE,
    size = 2.2,
    stroke = 0.35,
    color = "black",
    fill = "#D95F02"
  ) +
  facet_wrap(~ panel_label, ncol = 2) +
  suitability_scale +
  scale_shape_manual(values = organ_shapes, name = "Fossil organ") +
  coord_quickmap(xlim = c(-130, 150), ylim = c(15, 60), expand = FALSE) +
  labs(
    title = "Fossil records over matched paleoclimate projections",
    subtitle = "Paleo-coordinates are plotted on their corresponding TEYF projection period"
  ) +
  theme_map(base_size = 8) +
  theme(panel.spacing = unit(0.6, "lines"))

save_figure(
  p_r4,
  path_here("output", "figures_paper", "fossil_comparison", "R4_fossil_projection_panels"),
  width = 7.2,
  height = 5.8
)

# -------------------------------------------------------------------------
# R5. Fossil-point suitability values
# -------------------------------------------------------------------------

fossil_plot_data <- fossil_suitability %>%
  mutate(
    fossil_label = factor(fossil_label, levels = rev(unique(fossil_label))),
    suitability_class = factor(suitability_class, levels = suitability_levels)
  )

threshold_bands <- tibble::tibble(
  ymin = c(0, 0.1, 0.3, 0.7),
  ymax = c(0.1, 0.3, 0.7, 1),
  suitability_class = factor(suitability_levels, levels = suitability_levels)
)

p_r5 <- ggplot() +
  geom_rect(
    data = threshold_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = suitability_class),
    alpha = 0.18,
    color = NA
  ) +
  geom_hline(yintercept = c(0.1, 0.3, 0.7), color = "grey40", linewidth = 0.25, linetype = "dashed") +
  geom_point(
    data = fossil_plot_data,
    aes(x = fossil_label, y = extracted_suitability, fill = suitability_class, shape = organ),
    size = 2.8,
    stroke = 0.35,
    color = "black"
  ) +
  coord_flip(ylim = c(0, 1)) +
  scale_fill_manual(
    values = suitability_colors,
    labels = c(
      "Least (0-0.1)",
      "Slight (0.1-0.3)",
      "Moderate (0.3-0.7)",
      "High (0.7-1)"
    ),
    name = "Suitability class",
    guide = guide_legend(override.aes = list(shape = 22, size = 5, color = NA, alpha = 1)),
    drop = FALSE
  ) +
  scale_shape_manual(values = organ_shapes, name = "Fossil organ") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Fossil-point suitability values",
    subtitle = "Suitability extracted from the matched-period projection at each fossil paleo-coordinate",
    x = NULL,
    y = "Extracted suitability",
    caption = "Fossil record order follows the source table."
  ) +
  theme_paper(base_size = 8) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(color = "grey35", hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.2)
  )

save_figure(
  p_r5,
  path_here("output", "figures_paper", "fossil_comparison", "R5_fossil_point_suitability"),
  width = 7.2,
  height = 5.8
)

message("R4 and R5 fossil comparison figures were generated.")
