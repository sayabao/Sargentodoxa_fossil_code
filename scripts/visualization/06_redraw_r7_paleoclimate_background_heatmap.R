# Redraw R7 heatmap from the synchronized server summary table
#
# Purpose:
#   - Order variables from top to bottom by M4 variable importance.
#   - Use a taller output height matching M4 (4.8 in) so the heatmap is not flat.
#
# Inputs:
#   output/figures_paper/tables/R7_important_variable_paleoclimate_background.csv
#
# Outputs:
#   output/figures_paper/trends/R7_important_variable_paleoclimate_background_heatmap.pdf
#   output/figures_paper/trends/R7_important_variable_paleoclimate_background_heatmap.png

library(dplyr)
library(readr)
library(ggplot2)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "visualization") {
  project_root <- normalizePath(file.path(project_root, "..", ".."), mustWork = TRUE)
} else if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

path_here <- function(...) file.path(project_root, ...)

dir.create(path_here("output", "figures_paper", "trends"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_here("output", "figures_paper", "tables"), recursive = TRUE, showWarnings = FALSE)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      axis.title = element_text(color = "grey15"),
      axis.text = element_text(color = "grey20"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(6, 7, 6, 7)
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

r7_table <- read_csv(
  path_here("output", "figures_paper", "tables", "R7_important_variable_paleoclimate_background.csv"),
  show_col_types = FALSE
)

variable_order <- r7_table %>%
  distinct(variable, mean_importance) %>%
  arrange(desc(mean_importance)) %>%
  pull(variable)

period_order <- r7_table %>%
  distinct(Ma) %>%
  arrange(desc(Ma)) %>%
  mutate(period_label = paste0(Ma, " Ma")) %>%
  pull(period_label)

plot_data <- r7_table %>%
  filter(!is.na(mean_high)) %>%
  group_by(variable) %>%
  mutate(
    mean_high_z = as.numeric(scale(mean_high)),
    period_label = paste0(Ma, " Ma")
  ) %>%
  ungroup() %>%
  mutate(
    variable = factor(variable, levels = rev(variable_order)),
    period_label = factor(period_label, levels = period_order)
  )

r7_heatmap_order <- tibble::tibble(
  variable = variable_order,
  importance_rank = seq_along(variable_order)
)

write_csv(
  r7_heatmap_order,
  path_here("output", "figures_paper", "tables", "R7_heatmap_variable_order.csv")
)

p_heat <- ggplot(plot_data, aes(x = period_label, y = variable, fill = mean_high_z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3C65A6",
    mid = "white",
    high = "#B33A3A",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    title = "High-suitability paleoclimate background",
    subtitle = "Variables are ordered from top to bottom by ensemble importance; cells use suitability >= 0.7",
    x = "Time slice",
    y = NULL
  ) +
  theme_paper(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.key.width = unit(1.4, "cm"),
    legend.key.height = unit(0.35, "cm")
  )

save_figure(
  p_heat,
  path_here("output", "figures_paper", "trends", "R7_important_variable_paleoclimate_background_heatmap"),
  width = 6.8,
  height = 4.8
)

message("Variable order, top to bottom: ", paste(variable_order, collapse = " > "))
message("Done.")
