# Server-side paper figure: paleoclimate background of important variables
#
# This script needs the original TEYF bioclimatic rasters/nc files, which are
# not stored in the local project checkout. Run it on the server that contains:
#   /home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min/
#
# RStudio use:
#   1. Open this file on the server.
#   2. Edit the "RStudio/direct-run settings" block below if needed.
#   3. Click Source.
#
# Command-line example:
# Rscript scripts/visualization/04_server_important_variable_paleoclimate_background.R \
#   --teyf-climate-dir "/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min" \
#   --project-root "/home/cedrus/path/to/proj_sargentodoxa_niche_sim" \
#   --top-n 4 \
#   --threshold 0.7

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(terra)

# ============================================================================
# RStudio/direct-run settings
# ============================================================================

# Leave project_root blank to auto-detect it from the current working directory
# or from the folder containing this script. Existing scripts/config.R uses:
#   project_dir <- "/home/staff/shz/proj_sargentodoxa"
# If auto-detection fails, set this to the absolute server project path.
rstudio_project_root <- ""

# Server directory containing the original TEYF Bio nc/tif files.
rstudio_teyf_climate_dir <- "/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min"

# Existing scripts/config.R uses this simulation table path.
rstudio_teyf_simulations_file <- "/home/cedrus/data/shz/data/from linlin/teyf_simulations.csv"

# Leave blank to auto-search common projection folders under project_root.
# If needed, set this to the directory that contains files like:
#   teyfo_Sargentodoxa_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo.tif
rstudio_projection_dir <- ""

# Default variable selection: top N variables from ensemble variable importance.
rstudio_top_n <- 4L

# High-suitability threshold, matching current R1/R2/R3/R4 plotting logic.
rstudio_threshold <- 0.7

# Optional manual variables. Example: c("Bio11", "Bio9", "Bio15", "Bio18")
# Leave character(0) to use top N by variable importance.
rstudio_variables <- character(0)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      i <- i + 1
      next
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  markers <- c(
    file.path("data", "teyf_simulations.csv"),
    file.path("output", "07_ensemble", "ensemble_variable_importance.csv")
  )

  repeat {
    if (any(file.exists(file.path(current, markers)))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  script_name <- "04_server_important_variable_paleoclimate_background.R"
  script_hit <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (!is.na(script_hit) && basename(script_hit) == script_name) {
    candidate <- normalizePath(file.path(dirname(script_hit), "..", ".."), mustWork = FALSE)
    if (any(file.exists(file.path(candidate, markers)))) {
      return(candidate)
    }
  }

  stop(
    "Could not auto-detect project root. Set rstudio_project_root at the top of this script."
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(args$help)) {
  cat("Usage from command line:\n")
  cat("  Rscript scripts/visualization/04_server_important_variable_paleoclimate_background.R \\\n")
  cat("    --teyf-climate-dir \"/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min\" \\\n")
  cat("    --project-root \"/path/to/proj_sargentodoxa_niche_sim\" \\\n")
  cat("    --top-n 4 --threshold 0.7\n\n")
  cat("Usage from RStudio:\n")
  cat("  Edit the RStudio/direct-run settings block at the top, then click Source.\n")
  quit(status = 0)
}

project_root <- if (!is.null(args[["project-root"]])) {
  args[["project-root"]]
} else if (nzchar(rstudio_project_root)) {
  rstudio_project_root
} else {
  find_project_root(getwd())
}
project_root <- normalizePath(project_root, mustWork = TRUE)
teyf_climate_dir <- if (!is.null(args[["teyf-climate-dir"]])) {
  args[["teyf-climate-dir"]]
} else if (nzchar(rstudio_teyf_climate_dir)) {
  rstudio_teyf_climate_dir
} else {
  "/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min"
}
teyf_climate_dir <- normalizePath(teyf_climate_dir, mustWork = TRUE)

top_n <- if (!is.null(args[["top-n"]])) as.integer(args[["top-n"]]) else as.integer(rstudio_top_n)
threshold <- if (!is.null(args[["threshold"]])) as.numeric(args[["threshold"]]) else as.numeric(rstudio_threshold)
manual_variables <- if (!is.null(args[["variables"]])) {
  trimws(strsplit(args[["variables"]], ",", fixed = TRUE)[[1]])
} else {
  rstudio_variables
}

projection_dir <- if (!is.null(args[["projection-dir"]])) {
  args[["projection-dir"]]
} else if (nzchar(rstudio_projection_dir)) {
  rstudio_projection_dir
} else {
  ""
}
projection_dir <- if (nzchar(projection_dir)) normalizePath(projection_dir, mustWork = TRUE) else ""

teyf_simulations_file <- if (!is.null(args[["teyf-simulations-file"]])) {
  args[["teyf-simulations-file"]]
} else if (nzchar(rstudio_teyf_simulations_file)) {
  rstudio_teyf_simulations_file
} else {
  file.path(project_root, "data", "teyf_simulations.csv")
}
teyf_simulations_file <- normalizePath(teyf_simulations_file, mustWork = TRUE)

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
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
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

normalize_prediction_raster <- function(r) {
  xmax <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (is.finite(xmax) && xmax > 1) {
    r <- r / 1000
  }
  terra::clamp(r, lower = 0, upper = 1, values = TRUE)
}

projection_dir_candidates <- function() {
  candidates <- c(
    projection_dir,
    path_here("output", "09_paleoclimate_teyf", "projections"),
    path_here("output", "09_paleoclimate", "projections"),
    path_here("output", "09_paleoclimate_teyf"),
    path_here("output", "09_paleoclimate"),
    path_here("Models")
  )
  unique(candidates[nzchar(candidates) & dir.exists(candidates)])
}

list_projection_examples <- function(search_dirs) {
  examples <- unlist(lapply(search_dirs, function(d) {
    list.files(d, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  examples <- unique(examples)
  if (length(examples) == 0) {
    return("No .tif files found in candidate projection directories.")
  }
  paste(head(examples, 8), collapse = "\n    ")
}

find_projection_file <- function(period) {
  search_dirs <- projection_dir_candidates()
  if (length(search_dirs) == 0) {
    stop("No projection directories found. Set rstudio_projection_dir at the top of this script.")
  }

  hit <- unlist(lapply(search_dirs, function(d) {
    list.files(
      d,
      pattern = paste0("^", period, "_.*\\.tif$"),
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE)

  if (length(hit) == 0) {
    hit <- unlist(lapply(search_dirs, function(d) {
      files <- list.files(d, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      files[grepl(period, basename(files), ignore.case = TRUE) & grepl("EMwmean", basename(files), ignore.case = TRUE)]
    }), use.names = FALSE)
  }

  if (length(hit) == 0) {
    stop(
      "No suitability projection found for period: ", period, "\n",
      "Searched directories:\n    ", paste(search_dirs, collapse = "\n    "), "\n",
      "Example .tif files found:\n    ", list_projection_examples(search_dirs)
    )
  }

  hit <- unique(hit)
  em_hit <- hit[grepl("EMwmean", basename(hit), ignore.case = TRUE)]
  if (length(em_hit) > 0) {
    return(em_hit[1])
  }
  hit[1]
}

period_dir_candidates <- function(period) {
  c(
    file.path(teyf_climate_dir, period),
    file.path(teyf_climate_dir, paste0(period, "_masked")),
    file.path(teyf_climate_dir, toupper(period)),
    file.path(teyf_climate_dir, paste0(toupper(period), "_masked"))
  )
}

find_period_dir <- function(period) {
  candidates <- period_dir_candidates(period)
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) {
    stop("No climate directory found for period: ", period)
  }
  hit[1]
}

bio_regex <- function(variable) {
  number <- sub("^[Bb]io", "", variable)
  paste0("(^|[^[:alnum:]])bio_?0*", number, "([^[:alnum:]]|$)")
}

standardize_bio_names <- function(x) {
  out <- gsub(".*[Bb]io_?0*(\\d+).*", "Bio\\1", x)
  out
}

find_variable_file <- function(period, variable) {
  pdir <- find_period_dir(period)
  files <- list.files(pdir, pattern = "\\.(nc|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No climate rasters found in: ", pdir)
  }
  names_to_match <- tolower(basename(files))
  hit <- files[grepl(bio_regex(variable), names_to_match, perl = TRUE)]
  if (length(hit) > 0) {
    return(hit[1])
  }
  stop("Could not find ", variable, " file for ", period, " in ", pdir)
}

read_variable_raster <- function(period, variable) {
  variable_file <- find_variable_file(period, variable)
  r <- terra::rast(variable_file)
  names(r) <- standardize_bio_names(names(r))
  if (variable %in% names(r)) {
    r <- r[[variable]]
  } else if (terra::nlyr(r) > 1) {
    message("  Warning: ", variable_file, " has multiple layers; using first layer for ", variable)
    r <- r[[1]]
  }
  names(r) <- variable
  r
}

align_to_suitability <- function(variable_raster, suitability_raster) {
  same_geometry <- tryCatch(
    terra::compareGeom(variable_raster, suitability_raster, stopOnError = FALSE),
    error = function(e) FALSE
  )
  if (!isTRUE(same_geometry)) {
    variable_raster <- terra::resample(variable_raster, suitability_raster, method = "bilinear")
  }
  variable_raster
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

summarise_one <- function(period, ma, variable) {
  message("Processing ", period, " / ", ma, " Ma / ", variable)
  suitability <- normalize_prediction_raster(terra::rast(find_projection_file(period)))
  variable_raster <- read_variable_raster(period, variable)
  variable_raster <- align_to_suitability(variable_raster, suitability)

  high_mask <- suitability >= threshold
  high_values <- terra::values(terra::mask(variable_raster, high_mask, maskvalues = FALSE), mat = FALSE)
  all_values <- terra::values(variable_raster, mat = FALSE)

  tibble(
    simulation_name = period,
    Ma = ma,
    variable = variable,
    threshold = threshold,
    n_high_cells = sum(is.finite(high_values)),
    n_all_cells = sum(is.finite(all_values)),
    mean_high = mean(high_values, na.rm = TRUE),
    median_high = median(high_values, na.rm = TRUE),
    q25_high = safe_quantile(high_values, 0.25),
    q75_high = safe_quantile(high_values, 0.75),
    min_high = suppressWarnings(min(high_values, na.rm = TRUE)),
    max_high = suppressWarnings(max(high_values, na.rm = TRUE)),
    mean_all = mean(all_values, na.rm = TRUE),
    median_all = median(all_values, na.rm = TRUE),
    q25_all = safe_quantile(all_values, 0.25),
    q75_all = safe_quantile(all_values, 0.75)
  )
}

importance_file <- path_here("output", "07_ensemble", "ensemble_variable_importance.csv")
if (!file.exists(importance_file)) {
  importance_file <- path_here("output", "05_models", "variable_importance.csv")
}
importance <- read_csv(importance_file, show_col_types = FALSE)

importance_summary <- importance %>%
  mutate(variable = if ("expl.var" %in% names(.)) .data[["expl.var"]] else .data[["variable"]]) %>%
  mutate(importance = if ("var.imp" %in% names(.)) .data[["var.imp"]] else .data[["mean_importance"]]) %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_importance))

variables_to_use <- if (length(manual_variables) > 0) {
  manual_variables
} else {
  head(importance_summary$variable, top_n)
}

message("Project root: ", project_root)
message("TEYF climate dir: ", teyf_climate_dir)
if (nzchar(projection_dir)) {
  message("Projection dir: ", projection_dir)
} else {
  message("Projection dir: auto-search")
  message("Projection search dirs: ", paste(projection_dir_candidates(), collapse = " | "))
}
message("TEYF simulations file: ", teyf_simulations_file)
message("High-suitability threshold: ", threshold)
message("Variables to summarize: ", paste(variables_to_use, collapse = ", "))

teyf_simulations <- read_csv(teyf_simulations_file, show_col_types = FALSE) %>%
  transmute(
    simulation_name = `Simulation Name`,
    Ma = as.numeric(Ma)
  ) %>%
  filter(Ma <= 69) %>%
  arrange(desc(Ma))

stats <- bind_rows(lapply(seq_len(nrow(teyf_simulations)), function(i) {
  bind_rows(lapply(variables_to_use, function(v) {
    tryCatch(
      summarise_one(teyf_simulations$simulation_name[i], teyf_simulations$Ma[i], v),
      error = function(e) {
        message("  Error: ", conditionMessage(e))
        tibble(
          simulation_name = teyf_simulations$simulation_name[i],
          Ma = teyf_simulations$Ma[i],
          variable = v,
          threshold = threshold,
          error = conditionMessage(e)
        )
      }
    )
  }))
}))

stats <- stats %>%
  left_join(importance_summary, by = c("variable" = "variable")) %>%
  arrange(desc(Ma), desc(mean_importance), variable)

write_csv(
  stats,
  path_here("output", "figures_paper", "tables", "R7_important_variable_paleoclimate_background.csv")
)

if (!"mean_high" %in% names(stats)) {
  stop("No successful variable summaries were produced. Check climate file names and period directories.")
}

plot_data <- stats %>%
  filter(!is.na(mean_high), !is.na(mean_all)) %>%
  group_by(variable) %>%
  mutate(
    mean_high_z = as.numeric(scale(mean_high)),
    mean_all_z = as.numeric(scale(mean_all)),
    q25_high_z = (q25_high - mean(mean_high, na.rm = TRUE)) / sd(mean_high, na.rm = TRUE),
    q75_high_z = (q75_high - mean(mean_high, na.rm = TRUE)) / sd(mean_high, na.rm = TRUE),
    period_label = paste0(Ma, " Ma"),
    period_index = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    variable = factor(variable, levels = variables_to_use),
    period_label = factor(period_label, levels = unique(paste0(teyf_simulations$Ma, " Ma")))
  )

if (nrow(plot_data) == 0) {
  stop("No valid rows available for plotting. The summary table was still written for debugging.")
}

p_line <- ggplot(plot_data, aes(x = period_label, group = variable)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
  geom_line(aes(y = mean_all_z), color = "grey60", linewidth = 0.45, linetype = "22") +
  geom_line(aes(y = mean_high_z), color = "#1F6B3A", linewidth = 0.75) +
  geom_point(aes(y = mean_high_z), color = "#1F6B3A", fill = "white", shape = 21, size = 1.6, stroke = 0.35) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y") +
  labs(
    title = "Paleoclimate background of important variables",
    subtitle = paste0("Green line: mean value within suitability >= ", threshold, "; grey dashed line: full valid climate background"),
    x = "Time slice",
    y = "Standardized variable value"
  ) +
  theme_paper(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.spacing = unit(0.8, "lines")
  )

save_figure(
  p_line,
  path_here("output", "figures_paper", "trends", "R7_important_variable_paleoclimate_background_zscore"),
  width = 7.2,
  height = 5.8
)

p_heat <- ggplot(plot_data, aes(x = period_label, y = variable, fill = mean_high_z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#3C65A6",
    mid = "white",
    high = "#B33A3A",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    title = "High-suitability paleoclimate background",
    subtitle = paste0("Cells with suitability >= ", threshold, "; values standardized within each variable"),
    x = "Time slice",
    y = NULL
  ) +
  theme_paper(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom"
  )

save_figure(
  p_heat,
  path_here("output", "figures_paper", "trends", "R7_important_variable_paleoclimate_background_heatmap"),
  width = 7.2,
  height = 2.9
)

message("Done.")
