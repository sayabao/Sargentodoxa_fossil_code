# R3c. Latitudinal area-gradient heatmap for suitability >= 0.3
#
# Draws the manuscript's latitudinal heatmap of suitable habitat for the
# moderately + highly suitable classes (suitability >= 0.3), from the earliest
# fossil period to present.
#
# It reads the already-computed latitudinal_suitability_area_by_class.csv
# (which contains all four classes) and does NOT recompute rasters.
# It honours SARGENTO_OUTPUT_ROOT to point at the desired output tree
# (e.g. output_new, output_new_2/output).
#
# Usage:
#   export SARGENTO_OUTPUT_ROOT=output_new_2/output
#   Rscript scripts/visualization/R3c_latitudinal_area_heatmap.R

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "visualization") {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
} else if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

output_root <- Sys.getenv("SARGENTO_OUTPUT_ROOT", "output")
path_here <- function(...) {
  parts <- list(...)
  if (length(parts) >= 1L && identical(parts[[1]], "output")) parts[[1]] <- output_root
  do.call(file.path, c(list(project_root), parts))
}

dir.create(path_here("output", "figures_paper", "trends"), recursive = TRUE, showWarnings = FALSE)

theme_paper <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

save_figure <- function(plot, stem, width, height, dpi = 600) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Wrote: ", stem, ".pdf / .png")
}

lat_area <- read_csv(
  path_here("output", "figures_paper", "tables", "latitudinal_suitability_area_by_class.csv"),
  show_col_types = FALSE
)

periods <- lat_area %>%
  distinct(period, Ma, period_order, ma_label) %>%
  arrange(desc(Ma))

# Suitability >= 0.3 (moderately + highly suitable)
panel_ge03 <- lat_area %>%
  filter(suitability_class %in% c("Moderately suitable", "Highly suitable")) %>%
  group_by(period, Ma, period_order, ma_label, lat_bin, lat_bin_mid, lat_band) %>%
  summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(panel = ">= 0.3 (moderate + high)")

# ---- R3c: latitudinal area gradient for suitability >= 0.3 ----
p_r3c <- ggplot(panel_ge03, aes(x = period_order, y = lat_bin_mid, fill = area_km2 / 1e6)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_x_continuous(
    breaks = periods$period_order,
    labels = periods$ma_label,
    minor_breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(-80, 80, by = 20),
    limits = c(-90, 90),
    expand = c(0, 0)
  ) +
  scale_fill_gradientn(
    colours = c("#F7FCF5", "#C7E9C0", "#74C476", "#238B45", "#00441B"),
    name = expression("Suitable area ("*10^6~km^2*")"),
    trans = "sqrt"
  ) +
  labs(
    title = "Latitudinal area gradient of suitable habitat (suitability >= 0.3)",
    subtitle = paste0("Moderately + highly suitable (0.3-1); 10-degree latitude bands from earliest fossil period to present"),
    x = "Age (Ma)",
    y = "Latitude band midpoint"
  ) +
  theme_paper()

save_figure(
  p_r3c,
  path_here("output", "figures_paper", "trends", "R3c_latitudinal_area_heatmap_ge03"),
  width = 7.4,
  height = 4.8
)

message("Done. Standalone R3c (0.3-1) written.")