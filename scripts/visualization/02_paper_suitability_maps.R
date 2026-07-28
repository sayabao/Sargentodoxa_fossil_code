# Paper figures: suitability maps
#
# Outputs:
#   R1 present-day suitability map
#   R2 key TEYF paleoclimate suitability panels
#   R2b all TEYF paleoclimate suitability panels from 69 Ma to present
#
# Note:
#   The modern biomod2 projection in output/08_projection is currently stored as
#   a biomod2 RData object with server-side raster links. For a reproducible
#   local figure, R1 uses the TEYF 0 Ma projection (teyfa), which is also used as
#   the present-day panel in R2.

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

dir.create(path_here("output", "figures_paper", "maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_here("output", "figures_paper", "tables"), recursive = TRUE, showWarnings = FALSE)

theme_map <- function(base_size = 9) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "grey35", margin = margin(b = 6)),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
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

raster_to_plot_df <- function(raster_file, aggregate_factor = 8) {
  r <- terra::rast(raster_file)
  r <- normalize_prediction_raster(r)
  r <- terra::aggregate(r, fact = aggregate_factor, fun = "mean", na.rm = TRUE)
  names(r) <- "suitability"
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    mutate(
      suitability = pmin(pmax(suitability, 0), 1),
      suitability_class = cut(
        suitability,
        breaks = c(0, 0.1, 0.3, 0.7, 1),
        labels = c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable"),
        include.lowest = TRUE,
        right = FALSE
      )
    )
}

class_suitability <- function(x) {
  cut(
    x,
    breaks = c(0, 0.1, 0.3, 0.7, 1),
    labels = c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable"),
    include.lowest = TRUE,
    right = FALSE
  )
}

calculate_area_by_class <- function(raster_file, period, ma) {
  r <- normalize_prediction_raster(terra::rast(raster_file))
  area_r <- terra::cellSize(r, unit = "km")
  names(r) <- "suitability"
  names(area_r) <- "area_km2"
  
  as.data.frame(c(r, area_r), na.rm = TRUE) %>%
    mutate(
      period = period,
      Ma = ma,
      suitability_class = class_suitability(pmin(pmax(suitability, 0), 1))
    ) %>%
    filter(!is.na(suitability_class)) %>%
    group_by(period, Ma, suitability_class) %>%
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
}

suitability_scale <- scale_fill_manual(
  values = c(
    "Least suitable" = "#E8E8E8",
    "Slightly suitable" = "#BFD8A9",
    "Moderately suitable" = "#66A65B",
    "Highly suitable" = "#1F6B3A"
  ),
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

plot_global_suitability <- function(plot_data) {
  ggplot(plot_data, aes(x = x, y = y, fill = suitability_class)) +
    geom_raster() +
    suitability_scale +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    theme_map()
}

teyf_simulations <- read_csv(path_here("documents", "teyf_simulations.csv"), show_col_types = FALSE)

key_periods <- tibble::tibble(
  period = c("teyfo", "teyfk", "teyfj", "teyfe", "teyfd", "teyfb", "teyfa"),
  Ma = c(69, 52, 44, 19, 14, 4, 0)
) %>%
  mutate(
    panel_label = paste0(Ma, " Ma"),
    raster_file = vapply(period, find_teyf_file, character(1))
  )

all_periods <- teyf_simulations %>%
  transmute(
    period = `Simulation Name`,
    Ma = as.numeric(Ma)
  ) %>%
  filter(is.finite(Ma), Ma <= 69) %>%
  arrange(desc(Ma)) %>%
  mutate(
    panel_label = paste0(Ma, " Ma"),
    raster_file = vapply(period, find_teyf_file, character(1))
  )

write_csv(all_periods, path_here("output", "figures_paper", "tables", "R2b_all_period_map_inputs.csv"))

write_csv(key_periods, path_here("output", "figures_paper", "tables", "R2_key_period_map_inputs.csv"))

key_period_area <- bind_rows(lapply(seq_len(nrow(key_periods)), function(i) {
  calculate_area_by_class(key_periods$raster_file[i], key_periods$period[i], key_periods$Ma[i])
}))

write_csv(
  key_period_area,
  path_here("output", "figures_paper", "tables", "R2_key_period_suitability_area_by_class.csv")
)

all_period_area <- bind_rows(lapply(seq_len(nrow(all_periods)), function(i) {
  calculate_area_by_class(all_periods$raster_file[i], all_periods$period[i], all_periods$Ma[i])
}))

write_csv(
  all_period_area,
  path_here("output", "figures_paper", "tables", "R2b_all_period_suitability_area_by_class.csv")
)

# -------------------------------------------------------------------------
# R1. Present-day suitability map
# -------------------------------------------------------------------------

present_file <- find_teyf_file("teyfa")
present_map_data <- raster_to_plot_df(present_file, aggregate_factor = 6)

p_r1 <- plot_global_suitability(present_map_data) +
  labs(
    title = "Present-day habitat suitability",
    subtitle = "TEYF 0 Ma ensemble projection, classified by shared suitability thresholds"
  )

save_figure(
  p_r1,
  path_here("output", "figures_paper", "maps", "R1_modern_suitability_map"),
  width = 7.8,
  height = 3.9
)

# -------------------------------------------------------------------------
# R2. Key paleoclimate suitability panels
# -------------------------------------------------------------------------

panel_data <- bind_rows(lapply(seq_len(nrow(key_periods)), function(i) {
  raster_to_plot_df(key_periods$raster_file[i], aggregate_factor = 8) %>%
    mutate(
      period = key_periods$period[i],
      Ma = key_periods$Ma[i],
      panel_label = key_periods$panel_label[i]
    )
}))

panel_data <- panel_data %>%
  mutate(panel_label = factor(panel_label, levels = key_periods$panel_label))

p_r2 <- ggplot(panel_data, aes(x = x, y = y, fill = suitability_class)) +
  geom_raster() +
  facet_wrap(~ panel_label, ncol = 2) +
  suitability_scale +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = "Paleoclimate habitat suitability through key periods",
    subtitle = "Selected TEYF periods from the earliest fossil-matched period to present, classified using shared suitability thresholds"
  ) +
  theme_map(base_size = 8) +
  theme(
    panel.spacing = unit(0.55, "lines"),
    strip.text = element_text(margin = margin(b = 2))
  )

save_figure(
  p_r2,
  path_here("output", "figures_paper", "maps", "R2_paleoclimate_key_period_panels"),
  width = 7.2,
  height = 8.9
)

# -------------------------------------------------------------------------
# R2b. All paleoclimate suitability panels from 69 Ma to present
# -------------------------------------------------------------------------

all_panel_data <- bind_rows(lapply(seq_len(nrow(all_periods)), function(i) {
  raster_to_plot_df(all_periods$raster_file[i], aggregate_factor = 9) %>%
    mutate(
      period = all_periods$period[i],
      Ma = all_periods$Ma[i],
      panel_label = all_periods$panel_label[i]
    )
}))

all_panel_data <- all_panel_data %>%
  mutate(panel_label = factor(panel_label, levels = all_periods$panel_label))

p_r2b <- ggplot(all_panel_data, aes(x = x, y = y, fill = suitability_class)) +
  geom_raster() +
  facet_wrap(~ panel_label, ncol = 5) +
  suitability_scale +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = "Paleoclimate habitat suitability across all modeled periods",
    subtitle = "All TEYF periods from 69 Ma to present, classified using shared suitability thresholds"
  ) +
  theme_map(base_size = 7) +
  theme(
    panel.spacing = unit(0.35, "lines"),
    strip.text = element_text(margin = margin(b = 1.5), size = 7.2),
    legend.position = "bottom",
    plot.margin = margin(3, 5, 3, 5)
  )

save_figure(
  p_r2b,
  path_here("output", "figures_paper", "maps", "R2b_all_period_suitability_panels"),
  width = 11,
  height = 5.6
)

message("R1, R2, and R2b suitability maps were generated.")
