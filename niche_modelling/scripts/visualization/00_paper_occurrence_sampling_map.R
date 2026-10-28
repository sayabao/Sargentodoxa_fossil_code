# Paper figure M1: occurrence and pseudo-absence sampling map
#
# Outputs:
#   M1 occurrence and pseudo-absence sampling map
#   M1 sampling point table and summary

library(dplyr)
library(readr)
library(ggplot2)
library(maps)

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

dir.create(path_here("output", "figures_paper", "model_diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_here("output", "figures_paper", "tables"), recursive = TRUE, showWarnings = FALSE)

theme_map <- function(base_size = 9) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "grey35", margin = margin(b = 6)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(0.65, "cm"),
      plot.margin = margin(6, 6, 6, 6)
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

read_sampling_points <- function() {
  final_points_file <- path_here("output", "04_formatted", "all_points.csv")
  if (file.exists(final_points_file)) {
    points <- read_csv(final_points_file, show_col_types = FALSE) %>%
      transmute(
        presence = as.integer(presence),
        cv_fold = if ("cv_fold" %in% names(.)) as.integer(cv_fold) else NA_integer_,
        x = as.numeric(x),
        y = as.numeric(y),
        source = "final_model_points"
      )
    return(points)
  }

  occ <- read_csv(path_here("output", "01_data_prep", "sargentodoxa_occ_cvh_clean.csv"), show_col_types = FALSE) %>%
    transmute(
      presence = 1L,
      cv_fold = NA_integer_,
      x = as.numeric(x),
      y = as.numeric(y),
      source = "presence_clean"
    )

  # Use the pseudo-absence scheme that 04 actually formatted (read from its
  # metadata); default to scheme_1to1 if the metadata is unavailable.
  pa_scheme <- "scheme_1to1"
  fmt_meta <- path_here("output", "04_formatted", "data_formatting_metadata.rds")
  tryCatch({
    m <- readRDS(fmt_meta)
    if (!is.null(m$pa_scheme)) pa_scheme <- m$pa_scheme
  }, error = function(e) {})
  pa_suffix <- sub("^scheme_", "", pa_scheme)
  pa_candidates <- c(
    path_here("output", "02_pseudo_absence", paste0("sargentodoxa_pseudoabsence_", pa_suffix, ".csv")),
    path_here("output", "02_pseudo_absence", "sargentodoxa_pseudoabsence_1to1.csv"),
    path_here("output", "02_pseudo_absence", "sargentodoxa_pseudoabsence_1000.csv")
  )
  pa_file <- pa_candidates[file.exists(pa_candidates)][1]
  pa <- read_csv(pa_file, show_col_types = FALSE) %>%
    transmute(
      presence = 0L,
      cv_fold = NA_integer_,
      x = as.numeric(x),
      y = as.numeric(y),
      source = "pseudo_absence"
    )

  bind_rows(occ, pa)
}

sampling_points <- read_sampling_points() %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(
    point_type = ifelse(presence == 1, "Occurrence", "Pseudo-absence"),
    point_type = factor(point_type, levels = c("Pseudo-absence", "Occurrence"))
  )

sampling_summary <- sampling_points %>%
  count(point_type, name = "n_points") %>%
  mutate(
    x_min = min(sampling_points$x),
    x_max = max(sampling_points$x),
    y_min = min(sampling_points$y),
    y_max = max(sampling_points$y)
  )

write_csv(
  sampling_points,
  path_here("output", "figures_paper", "tables", "M1_occurrence_pseudoabsence_points.csv")
)
write_csv(
  sampling_summary,
  path_here("output", "figures_paper", "tables", "M1_occurrence_pseudoabsence_summary.csv")
)

plot_xlim <- range(sampling_points$x, na.rm = TRUE) + c(-5, 5)
plot_ylim <- range(sampling_points$y, na.rm = TRUE) + c(-5, 5)
plot_xlim <- c(max(-180, plot_xlim[1]), min(180, plot_xlim[2]))
plot_ylim <- c(max(-90, plot_ylim[1]), min(90, plot_ylim[2]))

world_map <- map_data("world") %>%
  filter(
    long >= plot_xlim[1] - 5,
    long <= plot_xlim[2] + 5,
    lat >= plot_ylim[1] - 5,
    lat <= plot_ylim[2] + 5
  )

p_m1 <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey94",
    color = NA
  ) +
  geom_point(
    data = sampling_points %>% filter(point_type == "Pseudo-absence"),
    aes(x = x, y = y, color = point_type, shape = point_type),
    size = 1.2,
    alpha = 0.5,
    stroke = 0
  ) +
  geom_point(
    data = sampling_points %>% filter(point_type == "Occurrence"),
    aes(x = x, y = y, color = point_type, shape = point_type),
    size = 1.25,
    alpha = 0.82,
    stroke = 0.25
  ) +
  scale_color_manual(
    values = c("Pseudo-absence" = "#6C8EBF", "Occurrence" = "#D95F02"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Pseudo-absence" = 16, "Occurrence" = 17),
    name = NULL
  ) +
  coord_quickmap(xlim = plot_xlim, ylim = plot_ylim, expand = FALSE) +
  labs(
    title = "Occurrence and pseudo-absence sampling",
    subtitle = paste0(
      "Final modeling points: ",
      sum(sampling_points$point_type == "Occurrence"), " occurrences and ",
      sum(sampling_points$point_type == "Pseudo-absence"), " pseudo-absences"
    )
  ) +
  theme_map(base_size = 9)

save_figure(
  p_m1,
  path_here("output", "figures_paper", "model_diagnostics", "M1_occurrence_pseudoabsence_map"),
  width = 6.8,
  height = 4.8
)

message("Done.")
