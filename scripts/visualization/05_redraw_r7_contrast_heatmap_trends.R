# Redraw R7 as high-suitability vs background contrast
#
# This figure emphasizes how strongly the high-suitability cells differ from
# the full valid climate background for each important variable and period.
# It does not replace M4 variable importance; it complements M4 by showing
# temporal environmental differentiation for the important variables.
#
# Inputs:
#   output/figures_paper/tables/R7_important_variable_paleoclimate_background.csv
#
# Outputs:
#   output/figures_paper/trends/R7_important_variable_background_contrast_heatmap.pdf
#   output/figures_paper/trends/R7_important_variable_background_contrast_heatmap.png
#   output/figures_paper/tables/R7_important_variable_background_contrast.csv

library(dplyr)
library(readr)
library(ggplot2)

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

contrast_data <- r7_table %>%
  filter(!is.na(mean_high), !is.na(mean_all)) %>%
  group_by(variable) %>%
  mutate(
    pooled_mean = mean(c(mean_high, mean_all), na.rm = TRUE),
    pooled_sd = sd(c(mean_high, mean_all), na.rm = TRUE),
    high_z = (mean_high - pooled_mean) / pooled_sd,
    background_z = (mean_all - pooled_mean) / pooled_sd,
    signed_contrast_z = high_z - background_z,
    abs_contrast_z = abs(signed_contrast_z),
    period_label = paste0(Ma, " Ma")
  ) %>%
  ungroup() %>%
  mutate(
    variable = factor(variable, levels = rev(variable_order)),
    variable_line = factor(variable, levels = variable_order),
    period_label = factor(period_label, levels = period_order)
  )

write_csv(
  contrast_data,
  path_here("output", "figures_paper", "tables", "R7_important_variable_background_contrast.csv")
)

p_heat <- ggplot(contrast_data, aes(x = period_label, y = variable, fill = abs_contrast_z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FBF7F2", "#E7C8BF", "#C4776D", "#8B2C2B"),
    values = c(0, 0.35, 0.58, 0.78, 1),
    name = "|High - background|\nstandardized"
  ) +
  labs(
    title = "High-suitability vs background climate contrast",
    subtitle = "Heatmap shows absolute standardized difference; variables ordered by ensemble importance",
    x = NULL,
    y = NULL
  ) +
  coord_fixed(ratio = 1) +
  theme_paper(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.key.width = unit(1.8, "cm"),
    legend.key.height = unit(0.35, "cm")
  )

save_figure(
  p_heat,
  path_here("output", "figures_paper", "trends", "R7_important_variable_background_contrast_heatmap"),
  width = 6.8,
  height = 3.2
)

message("Variable order, top to bottom in heatmap: ", paste(variable_order, collapse = " > "))
message("Done.")
