#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(png)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

fig_dir <- file.path(root, "output", "figures_paper")
plate_dir <- file.path(fig_dir, "plates")
dir.create(plate_dir, recursive = TRUE, showWarnings = FALSE)

image_path <- function(...) {
  path <- file.path(fig_dir, ...)
  if (!file.exists(path)) {
    stop("Missing figure input: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

crop_top <- function(img, pixels = 0) {
  if (pixels <= 0) {
    return(img)
  }
  keep_from <- min(dim(img)[1], pixels + 1)
  img[keep_from:dim(img)[1], , , drop = FALSE]
}

panel_grob <- function(path, tag = NULL, tag_size = 16, crop_top_px = 0) {
  img <- png::readPNG(path)
  img <- crop_top(img, crop_top_px)
  grobs <- list(
    rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"),
               interpolate = TRUE)
  )

  if (!is.null(tag)) {
    grobs <- c(
      grobs,
      list(
        textGrob(
          tag,
          x = unit(0.012, "npc"),
          y = unit(0.985, "npc"),
          just = c("left", "top"),
          gp = gpar(fontsize = tag_size, fontface = "bold", col = "black")
        )
      )
    )
  }

  grobTree(children = do.call(gList, grobs))
}

save_plate <- function(grob, basename, width, height, dpi = 300) {
  pdf_path <- file.path(plate_dir, paste0(basename, ".pdf"))
  png_path <- file.path(plate_dir, paste0(basename, ".png"))

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  grid.newpage()
  grid.draw(grob)
  dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in",
                 res = dpi, type = "cairo-png")
  grid.newpage()
  grid.draw(grob)
  dev.off()

  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

single_plate <- function(src_png, basename, width, height, crop_top_px = 0) {
  grob <- panel_grob(src_png, crop_top_px = crop_top_px)
  save_plate(grob, basename, width, height)
}

# Main Figure 1: M4 + R7
fig1 <- arrangeGrob(
  panel_grob(image_path("model_diagnostics", "M4_variable_importance.png"), "A", crop_top_px = 265),
  panel_grob(image_path("trends", "R7_important_variable_background_contrast_heatmap.png"), "B", crop_top_px = 310),
  ncol = 1,
  heights = unit.c(unit(4.8, "null"), unit(3.2, "null"))
)
save_plate(fig1, "Figure_1_M4_R7_variable_importance_and_background", 6.8, 8.0)

# Main Figure 2: R3 + R3b
fig2 <- arrangeGrob(
  panel_grob(image_path("trends", "R3_suitability_area_trend_from_earliest_fossil.png"), "A", crop_top_px = 275),
  panel_grob(image_path("trends", "R3b_latitudinal_suitability_area_gradient.png"), "B", crop_top_px = 205),
  ncol = 1,
  heights = unit.c(unit(4.2, "null"), unit(3.8, "null"))
)
save_plate(fig2, "Figure_2_R3_R3b_temporal_latitudinal_dynamics", 7.2, 8.0)

# Main Figure 3: R2b all-period panels. This is already a complete plate.
single_plate(
  image_path("maps", "R2b_all_period_suitability_panels.png"),
  "Figure_3_R2b_all_period_suitability_panels",
  11.0,
  5.2,
  crop_top_px = 350
)

# Main Figure 4: R4 + R5
fig4 <- arrangeGrob(
  panel_grob(image_path("fossil_comparison", "R4_fossil_projection_panels.png"), "A", crop_top_px = 390),
  panel_grob(image_path("fossil_comparison", "R5_fossil_point_suitability.png"), "B", crop_top_px = 230),
  ncol = 1,
  heights = unit.c(unit(5.2, "null"), unit(3.2, "null"))
)
save_plate(fig4, "Figure_4_R4_R5_fossil_model_consistency", 7.4, 8.4)

# Appendix Figure S1: M1 + M2 + M3
figs1 <- arrangeGrob(
  panel_grob(image_path("model_diagnostics", "M1_occurrence_pseudoabsence_map.png"), "A", crop_top_px = 270),
  arrangeGrob(
    panel_grob(image_path("model_diagnostics", "M2_variable_correlation_vif.png"), "B", tag_size = 14, crop_top_px = 310),
    panel_grob(image_path("model_diagnostics", "M3_model_evaluation_metrics.png"), "C", tag_size = 14, crop_top_px = 250),
    ncol = 2,
    widths = unit.c(unit(1.15, "null"), unit(1, "null"))
  ),
  ncol = 1,
  heights = unit.c(unit(4.4, "null"), unit(4.2, "null"))
)
save_plate(figs1, "Appendix_Figure_S1_M1_M2_M3_model_diagnostics", 9.2, 8.6)

message("Done. Combined plates are in: ", plate_dir)
