# Paper supplementary analysis: fossil-point suitability significance test
# (Optimized for speed)
#
# Question
#   The fossil record of Sargentodoxa is sparse (16 occurrences across 8 TEYF
#   periods). This script tests whether those localities are non-random with
#   respect to predicted climatic suitability, i.e. whether they fall in
#   higher-suitability cells than expected by chance. Because the ensemble model
#   is calibrated only on the present distribution of S. cuneata, the fossils
#   constitute an independent validation of the palaeoclimate projections.
#
# Methods
#   1. Per fossil: extract suitability at the palaeo-coordinate from the matched
#      TEYF projection (identical logic to R5 in 03_paper_fossil_comparison_...R).
#   2. Background = non-NA (land) cells of the matched projection, area-weighted
#      by cellSize (km2), consistent with R2b. Primary background is Northern
#      Hemisphere land (lat > 0); global land is reported as a sensitivity check.
#   3. Monte Carlo randomization (9,999 iterations), fully vectorized: for each
#      period, draw k_p area-weighted random cells (inverse-CDF), use each cell's
#      precomputed area-weighted percentile rank, and pool all 16 draws. The null
#      distribution is the pooled mean percentile; observed = pooled mean of the
#      16 fossils. p = (1 + #{null >= observed}) / (N + 1).
#      A raw-mean-suitability null (Design B, NH) is also reported.
#   4. Binomial test: for thresholds 0.1, 0.3, 0.7, compare the number of fossil
#      cells at/above the threshold to the area-weighted background proportion
#      (pooled by fossil count per period).
#
# Outputs
#   output/figures_paper/tables/R8_fossil_significance_per_fossil.csv
#   output/figures_paper/tables/R8_fossil_significance_tests.csv
#   output/figures_paper/tables/R8_fossil_significance_background_by_period.csv
#   output/figures_paper/fossil_comparison/R8_fossil_significance_figures10.pdf|.png
#     (intended as Suppl. Fig. S10; number to be confirmed against current scheme)

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

# ---- helpers (kept consistent with 02/03) -------------------------------

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

suitability_levels <- c("Least suitable", "Slightly suitable", "Moderately suitable", "Highly suitable")

class_suitability <- function(x) {
  cut(x, breaks = c(0, 0.1, 0.3, 0.7, 1), labels = suitability_levels,
      include.lowest = TRUE, right = FALSE)
}

suitability_colors <- c(
  "Least suitable" = "#E8E8E8", "Slightly suitable" = "#BFD8A9",
  "Moderately suitable" = "#66A65B", "Highly suitable" = "#1F6B3A"
)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      axis.title = element_text(color = "grey15"),
      axis.text = element_text(color = "grey20"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}

save_figure <- function(plot, stem, width, height, dpi = 600) {
  pdf_file <- paste0(stem, ".pdf"); png_file <- paste0(stem, ".png")
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Wrote: ", pdf_file); message("Wrote: ", png_file)
}

# ---- data ---------------------------------------------------------------

teyf_simulations <- read_csv(path_here("documents", "teyf_simulations.csv"), show_col_types = FALSE)
fossil_records <- read_csv(path_here("data", "fossil_record", "FossilRecord-paleocor-v2.1.csv"), show_col_types = FALSE) %>%
  mutate(number = as.integer(number), simulation_name = as.character(simulation_name))

# Optimization: Extract suitability grouped by simulation to avoid repetitive file reads
extract_all_fossils <- function(fossil_records) {
  results <- list()
  unique_sims <- unique(fossil_records$simulation_name)
  for (sim in unique_sims) {
    sub_df <- fossil_records[fossil_records$simulation_name == sim, ]
    r <- normalize_prediction_raster(terra::rast(find_teyf_file(sim)))
    points <- terra::vect(
      data.frame(x = sub_df$paleoLong, y = sub_df$paleoLat),
      geom = c("x", "y"), crs = "EPSG:4326"
    )
    vs <- terra::extract(r, points)
    sub_df$extracted_suitability <- as.numeric(vs[, 2])
    results[[sim]] <- sub_df
  }
  bind_rows(results)
}

fossil_df <- extract_all_fossils(fossil_records) %>%
  left_join(
    teyf_simulations %>% transmute(simulation_name = `Simulation Name`, Ma = as.numeric(Ma)),
    by = "simulation_name"
  ) %>%
  mutate(
    Ma = ifelse(is.na(Ma), age, Ma),
    extracted_suitability = pmin(pmax(extracted_suitability, 0), 1),
    suitability_class = factor(class_suitability(extracted_suitability), levels = suitability_levels)
  ) %>%
  arrange(desc(Ma), number)

if (any(is.na(fossil_df$extracted_suitability))) {
  warning("Some fossil points returned NA suitability (fell on non-land cells). Check coordinates.")
}

period_info <- fossil_df %>%
  group_by(simulation_name, Ma) %>%
  summarise(k_fossils = dplyr::n(), .groups = "drop") %>%
  arrange(desc(Ma))

# ---- background per period (area-weighted land cells) -------------------

# Global cache to avoid computing cellSize multiple times if grids match
cached_cellSize <- NULL

build_background <- function(sim_name) {
  r <- normalize_prediction_raster(terra::rast(find_teyf_file(sim_name)))
  names(r) <- "suitability"
  vals <- as.numeric(terra::values(r))
  ys <- terra::yFromCell(r, seq_len(terra::ncell(r)))

  # Compute or reuse cached cellSize
  if (is.null(cached_cellSize) || !terra::compareGeom(cached_cellSize, r, lyrs = FALSE, stopOnError = FALSE)) {
    cached_cellSize <<- terra::cellSize(r, unit = "km")
  }
  areas <- as.numeric(terra::values(cached_cellSize))

  keep <- !is.na(vals) & !is.na(areas)
  list(suit = vals[keep], area = areas[keep], lat = ys[keep])
}

bg <- setNames(lapply(period_info$simulation_name, build_background), period_info$simulation_name)

bg_subset <- function(b, region) {
  if (region == "NH") {
    keep <- b$lat > 0
  } else {
    keep <- rep(TRUE, length(b$suit))
  }
  list(suit = b$suit[keep], area = b$area[keep])
}

# Pre-sort each background and precompute cell percentiles for O(1) simulation lookup.
prep_bg <- function(suit, w) {
  o <- order(suit)
  s <- suit[o]; w <- w[o]
  tot <- sum(w)
  cw <- c(0, cumsum(w))
  cumw <- cumsum(w) / tot
  n <- length(s)

  # Optimization: Precompute the exact percentile rank of each sorted cell
  # to bypass findInterval() completely during Monte Carlo sampling
  idx_left <- findInterval(s, s, left.open = TRUE)
  pct <- cw[idx_left + 1] / tot * 100

  list(s = s, w = w, cw = cw, tot = tot,
       cumw = cumw, n = n, pct = pct)
}

prep_NH     <- lapply(bg, function(b) { s <- bg_subset(b, "NH");     prep_bg(s$suit, s$area) })
prep_global <- lapply(bg, function(b) { s <- bg_subset(b, "global"); prep_bg(s$suit, s$area) })

# Area-weighted percentile rank of v, O(log n) for observed fossil points
weighted_percentile <- function(v, prep) {
  idx <- findInterval(v, prep$s, left.open = TRUE)
  prep$cw[idx + 1] / prep$tot * 100
}

weighted_prop_above <- function(thr, prep) {
  idx <- findInterval(thr, prep$s, left.open = TRUE)
  1 - prep$cw[idx + 1] / prep$tot
}

# Per-fossil percentiles (NH and global)
fossil_df <- fossil_df %>%
  rowwise() %>%
  mutate(
    pct_NH = weighted_percentile(extracted_suitability, prep_NH[[simulation_name]]),
    pct_global = weighted_percentile(extracted_suitability, prep_global[[simulation_name]])
  ) %>%
  ungroup()

# ---- Monte Carlo: pooled mean percentile (primary) ----------------------

n_iter <- 9999
set.seed(20260706)

# Optimization: Fully vectorized Monte Carlo simulation
# Eliminates nested replicate/lapply/vapply. Runs in milliseconds.
run_simulation_vectorized <- function(prep_list, n_iter) {
  total_fossils <- sum(period_info$k_fossils)
  null_pct_matrix <- matrix(nrow = n_iter, ncol = total_fossils)
  null_suit_matrix <- matrix(nrow = n_iter, ncol = total_fossils)

  col_start <- 1
  for (p in period_info$simulation_name) {
    prep <- prep_list[[p]]
    k <- period_info$k_fossils[period_info$simulation_name == p]

    # Generate random probabilities for all iterations at once
    u <- runif(k * n_iter)
    idx <- pmin(findInterval(u, prep$cumw, left.open = TRUE) + 1L, prep$n)

    # Extract precomputed values and percentiles in O(1) time
    sampled_suits <- prep$s[idx]
    sampled_pcts  <- prep$pct[idx]

    col_end <- col_start + k - 1
    null_suit_matrix[, col_start:col_end] <- matrix(sampled_suits, nrow = n_iter, ncol = k)
    null_pct_matrix[, col_start:col_end]  <- matrix(sampled_pcts, nrow = n_iter, ncol = k)

    col_start <- col_end + 1
  }

  list(
    pct_means = rowMeans(null_pct_matrix),
    suit_means = rowMeans(null_suit_matrix)
  )
}

# Run simulation
sim_NH     <- run_simulation_vectorized(prep_NH, n_iter)
sim_global <- run_simulation_vectorized(prep_global, n_iter)

mc_pct_NH     <- sim_NH$pct_means
mc_suit_NH    <- sim_NH$suit_means
mc_pct_global <- sim_global$pct_means

obs_pct_NH     <- mean(fossil_df$pct_NH)
obs_pct_global <- mean(fossil_df$pct_global)
obs_mean_suit  <- mean(fossil_df$extracted_suitability)

# ---- binomial threshold tests -------------------------------------------

binom_test_thr <- function(thr, region) {
  prep_list <- if (region == "NH") prep_NH else prep_global
  n_total <- sum(period_info$k_fossils)
  obs_count <- sum(fossil_df$extracted_suitability >= thr)
  exp_p <- sum(vapply(period_info$simulation_name, function(p) {
    k <- period_info$k_fossils[period_info$simulation_name == p]
    k * weighted_prop_above(thr, prep_list[[p]])
  }, numeric(1))) / n_total
  bt <- binom.test(obs_count, n_total, exp_p, alternative = "greater")
  tibble::tibble(
    threshold = thr, region = region,
    obs_count = obs_count, n_total = n_total,
    exp_prop = exp_p, p_value = bt$p.value
  )
}

binom_results <- bind_rows(
  binom_test_thr(0.1, "NH"), binom_test_thr(0.3, "NH"), binom_test_thr(0.7, "NH"),
  binom_test_thr(0.1, "global"), binom_test_thr(0.3, "global"), binom_test_thr(0.7, "global")
) %>%
  mutate(test = sprintf("binom >= %.1f", threshold))

# ---- assemble test summary ----------------------------------------------

p_one_sided <- function(null, obs) (1 + sum(null >= obs)) / (n_iter + 1)

tests_summary <- bind_rows(
  tibble::tibble(
    test = "mean percentile rank",
    region = "NH",
    observed = obs_pct_NH,
    null_mean = mean(mc_pct_NH),
    null_sd = sd(mc_pct_NH),
    p_value = p_one_sided(mc_pct_NH, obs_pct_NH)
  ),
  tibble::tibble(
    test = "mean percentile rank",
    region = "global",
    observed = obs_pct_global,
    null_mean = mean(mc_pct_global),
    null_sd = sd(mc_pct_global),
    p_value = p_one_sided(mc_pct_global, obs_pct_global)
  ),
  tibble::tibble(
    test = "mean suitability",
    region = "NH",
    observed = obs_mean_suit,
    null_mean = mean(mc_suit_NH),
    null_sd = sd(mc_suit_NH),
    p_value = p_one_sided(mc_suit_NH, obs_mean_suit)
  ),
  binom_results %>%
    transmute(test, region, observed = obs_count,
              null_mean = NA_real_, null_sd = NA_real_, p_value) %>%
    mutate(exp_prop = binom_results$exp_prop)
)

# ---- background-by-period table -----------------------------------------

bg_by_period <- bind_rows(lapply(period_info$simulation_name, function(p) {
  tibble::tibble(
    simulation_name = p,
    Ma = period_info$Ma[period_info$simulation_name == p],
    k_fossils = period_info$k_fossils[period_info$simulation_name == p],
    NH_land_area_km2 = prep_NH[[p]]$tot,
    NH_prop_ge_0.1 = weighted_prop_above(0.1, prep_NH[[p]]),
    NH_prop_ge_0.3 = weighted_prop_above(0.3, prep_NH[[p]]),
    NH_prop_ge_0.7 = weighted_prop_above(0.7, prep_NH[[p]]),
    global_land_area_km2 = prep_global[[p]]$tot,
    global_prop_ge_0.1 = weighted_prop_above(0.1, prep_global[[p]]),
    global_prop_ge_0.3 = weighted_prop_above(0.3, prep_global[[p]]),
    global_prop_ge_0.7 = weighted_prop_above(0.7, prep_global[[p]])
  )
})) %>% arrange(desc(Ma))

# ---- write tables -------------------------------------------------------

write_csv(
  fossil_df %>%
    select(number, fossil_sp, organ, country_region, simulation_name, Ma,
           paleoLong, paleoLat, extracted_suitability, suitability_class,
           pct_NH, pct_global),
  path_here("output", "figures_paper", "tables", "R8_fossil_significance_per_fossil.csv")
)

write_csv(
  tests_summary,
  path_here("output", "figures_paper", "tables", "R8_fossil_significance_tests.csv")
)

write_csv(
  bg_by_period,
  path_here("output", "figures_paper", "tables", "R8_fossil_significance_background_by_period.csv")
)

# ---- publication summary table (Suppl. Table S?) ------------------------
# Standard randomization-test reporting: observed, null mean +/- SD,
# null 95% range, z score, N permutations, one-tailed P, significance code.

sig_code <- function(p) {
  if (is.na(p)) ""
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

fmt_p <- function(p) ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p))

n_total <- sum(period_info$k_fossils)

mc_specs <- list(
  list(stat = "Mean percentile rank of fossil cells (%)",
       bg = "Northern Hemisphere land",
       obs = obs_pct_NH, null = mc_pct_NH, fmt = "%.1f"),
  list(stat = "Mean percentile rank of fossil cells (%)",
       bg = "Global land",
       obs = obs_pct_global, null = mc_pct_global, fmt = "%.1f"),
  list(stat = "Mean suitability of fossil cells (0-1)",
       bg = "Northern Hemisphere land",
       obs = obs_mean_suit, null = mc_suit_NH, fmt = "%.3f")
)

mc_rows <- bind_rows(lapply(mc_specs, function(x) {
  nm <- mean(x$null); nsd <- sd(x$null)
  pval <- p_one_sided(x$null, x$obs)
  ci_lo <- quantile(x$null, 0.025, names = FALSE)
  ci_hi <- quantile(x$null, 0.975, names = FALSE)
  tibble::tibble(
    test_type = "Monte Carlo randomization",
    statistic = x$stat,
    background = x$bg,
    observed = sprintf(x$fmt, x$obs),
    expected_null = sprintf(x$fmt, nm),
    null_sd = sprintf(x$fmt, nsd),
    null_95CI = sprintf("[%s, %s]", sprintf(x$fmt, ci_lo), sprintf(x$fmt, ci_hi)),
    z_score = sprintf("%.2f", (x$obs - nm) / nsd),
    n_perm_or_trials = n_iter,
    p_value = fmt_p(pval),
    significance = sig_code(pval)
  )
}))

binom_rows <- binom_results %>%
  mutate(significance = case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )) %>%
  transmute(
    test_type = "Binomial test",
    statistic = sprintf("Fossil cells >= %.1f suitability", threshold),
    background = ifelse(region == "NH", "Northern Hemisphere land", "Global land"),
    observed = sprintf("%d/%d (%.1f%%)", obs_count, n_total, obs_count / n_total * 100),
    expected_null = sprintf("%.1f/%d (%.1f%%)", exp_prop * n_total, n_total, exp_prop * 100),
    null_sd = "",
    null_95CI = "",
    z_score = "",
    n_perm_or_trials = n_total,
    p_value = fmt_p(p_value),
    significance = significance
  )

suppl_table <- bind_rows(mc_rows, binom_rows)

write_csv(
  suppl_table,
  path_here("output", "figures_paper", "tables", "R8_fossil_significance_suppl_table.csv")
)

# NH-only version for publication (global rows are a sensitivity check only)
suppl_table_NH <- suppl_table %>%
  filter(background == "Northern Hemisphere land")

write_csv(
  suppl_table_NH,
  path_here("output", "figures_paper", "tables", "R8_fossil_significance_suppl_table_NH.csv")
)

message("---- Suppl. Table (R8 summary, NH only) ----")
print(as.data.frame(suppl_table_NH))

message("---- Suppl. Table (R8 summary, full) ----")
print(as.data.frame(suppl_table))

# ---- Suppl. Fig. S10 ----------------------------------------------------

null_df <- tibble::tibble(mean_pct = mc_pct_NH)
binw <- diff(range(null_df$mean_pct)) / 60
bin_origin <- min(null_df$mean_pct)
x_right <- obs_pct_NH + 16

pA <- ggplot() +
  geom_histogram(data = null_df, aes(x = mean_pct), binwidth = binw,
                 boundary = bin_origin, fill = "grey78", color = "white", linewidth = 0.2) +
  annotate("rect", xmin = obs_pct_NH, xmax = x_right, ymin = -Inf, ymax = Inf,
           fill = "#D95F02", alpha = 0.18) +
  geom_vline(xintercept = obs_pct_NH, color = "#1F6B3A", linewidth = 0.8) +
  annotate("text", x = obs_pct_NH + 0.8, y = Inf, vjust = 2, hjust = 0,
           label = sprintf("Observed = %.1f%%\nP = %.4f", obs_pct_NH,
                           p_one_sided(mc_pct_NH, obs_pct_NH)),
           color = "#1F6B3A", size = 3.2, fontface = "bold") +
  scale_x_continuous(limits = c(min(null_df$mean_pct) - binw, x_right)) +
  labs(
    title = "A. Null distribution of mean percentile rank",
    subtitle = "9,999 area-weighted random draws from Northern Hemisphere land; observed = 16 fossils",
    x = "Mean percentile rank (NH land background, %)", y = "Count"
  ) +
  theme_paper(base_size = 9)

fossil_plot <- fossil_df %>% mutate(fossil_label = number)
mean_pct <- mean(fossil_df$pct_NH)

pB <- ggplot(fossil_plot, aes(x = number, y = pct_NH, fill = suitability_class)) +
  geom_hline(yintercept = 50, color = "grey45", linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = mean_pct, color = "#1F6B3A", linewidth = 0.5) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.3) +
  scale_fill_manual(
    values = suitability_colors, name = "Suitability class",
    labels = c("Least (0-0.1)", "Slight (0.1-0.3)", "Moderate (0.3-0.7)", "High (0.7-1)"),
    drop = FALSE
  ) +
  scale_x_continuous(breaks = 1:16) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "B. Per-fossil percentile rank (NH land background)",
    subtitle = sprintf("Dashed = 50%% (null median); solid = observed mean %.1f%%", mean_pct),
    x = "Fossil record number", y = "Percentile rank (%)"
  ) +
  theme_paper(base_size = 9) +
  theme(legend.position = "right")

p_s10 <- pA / pB + plot_layout(heights = c(1, 1))

save_figure(
  p_s10,
  path_here("output", "figures_paper", "fossil_comparison", "R8_fossil_significance_figures10"),
  width = 7.4,
  height = 6.8
)

# ---- console summary ----------------------------------------------------

message("============================================================")
message("Fossil-point suitability significance (R8)")
message("------------------------------------------------------------")
message(sprintf("N fossils = %d across %d periods", nrow(fossil_df), nrow(period_info)))
message(sprintf("Observed mean suitability = %.3f", obs_mean_suit))
message(sprintf("Mean percentile rank  NH     : observed = %.1f%%  null = %.1f%% (SD %.1f)  P = %.4f",
                obs_pct_NH, mean(mc_pct_NH), sd(mc_pct_NH), p_one_sided(mc_pct_NH, obs_pct_NH)))
message(sprintf("Mean percentile rank  global : observed = %.1f%%  null = %.1f%% (SD %.1f)  P = %.4f",
                obs_pct_global, mean(mc_pct_global), sd(mc_pct_global), p_one_sided(mc_pct_global, obs_pct_global)))
message(sprintf("Mean suitability       NH     : observed = %.3f  null = %.3f (SD %.3f)  P = %.4f",
                obs_mean_suit, mean(mc_suit_NH), sd(mc_suit_NH), p_one_sided(mc_suit_NH, obs_mean_suit)))
for (i in seq_len(nrow(binom_results))) {
  r <- binom_results[i, ]
  message(sprintf("Binom >= %.1f  %-7s : %d/%d  expected prop = %.3f  P = %.4f",
                  r$threshold, r$region, r$obs_count, r$n_total, r$exp_prop, r$p_value))
}
message("============================================================")

# Suppl. Fig. S10 caption (paste into docs/PAPER_FIGURE_CAPTIONS.md):
#
# Supplementary Figure S10. Fossil localities are non-random with respect to
# predicted climatic suitability. (A) Null distribution of the pooled mean
# percentile rank of 16 fossil occurrences, generated by 9,999 area-weighted
# random draws from the Northern Hemisphere land background of each matched
# TEYF palaeoclimate slice; the green line marks the observed mean (82.2%),
# which lies beyond the null distribution, and the orange shading marks the
# one-sided tail region at or above it (P = 0.0001). (B) Per-fossil percentile
# rank within the matched Northern Hemisphere land background, coloured by
# suitability class; the dashed line is 50% (null median) and the solid line is
# the observed mean. Although few, the fossil localities occupy significantly
# higher-suitability cells than random expectation, independently corroborating
# the palaeoclimate projections.
