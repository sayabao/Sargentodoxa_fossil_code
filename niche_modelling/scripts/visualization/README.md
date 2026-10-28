# Visualization Scripts

Paper-figure component scripts. Each writes paper-ready outputs to
`output/figures_paper/`. The author manually assembles the final article
figures in Photoshop from these components.

Recommended pattern:

- Read existing results from `output/`.
- Keep shared figure style settings in one small config script if needed.
- Write paper-ready outputs to `output/figures_paper/`.
- Avoid changing the numbered modeling scripts unless the figure requires a real workflow fix.

Priority inputs:

- `output/04_formatted/all_points.csv` (M1)
- `output/03_variable_select/` (M2)
- `output/06_evaluation/`, `output/07_ensemble/` (M3, M4)
- `output/08_projection/` (Fig S4)
- `output/09_paleoclimate/` (R2b, R3, R3c, R4, R5, R7, R8)
- `data/fossil_record/FossilRecord-paleocor-v2.1.csv` (R3, R4, R5, R8)

Current scripts:

- `00_paper_occurrence_sampling_map.R`: M1 occurrence and pseudo-absence sampling map (Fig S1) from `output/04_formatted/all_points.csv`.
- `01_paper_model_and_trend_figures.R`: M2 variable filtering/VIF (Fig S2), M3 model evaluation (Fig S3), M4 variable importance (Fig 3 A), R3 suitability area trend (Fig 4 B), and the latitudinal suitability-area-by-class table (Fig 4 C input) from 69 Ma to present. The manuscript's latitudinal HEATMAP (suitability `>= 0.3`) is drawn by `R3c_latitudinal_area_heatmap.R`.
- `02_paper_suitability_maps.R`: R2b all-period TEYF suitability panels from 69 Ma to present (Fig 4 A) using the shared 0-0.1 / 0.1-0.3 / 0.3-0.7 / 0.7-1 classes.
- `03_paper_fossil_comparison_figures.R`: R4 fossil-projection panels (Fig 5 A) and R5 fossil-point suitability extraction (Fig 5 B) using the same suitability classes.
- `04_server_important_variable_paleoclimate_background.R`: server-side R7 script summarizing important Bio variables within high-suitability (`>= 0.7`) cells through time. RStudio-friendly: edit the settings block at the top and click Source. Requires the original TEYF Bio nc/tif files. Its output CSV is the input to `05_…`.
- `05_redraw_r7_contrast_heatmap_trends.R`: redraws the R7 high-suitability vs background contrast heatmap locally from the R7 CSV (Fig 3 B).
- `06_fossil_suitability_significance.R`: R8 fossil-point suitability significance test (Monte Carlo + binomial) and figure (Fig S11).
- `R3c_latitudinal_area_heatmap.R`: draws the manuscript's latitudinal heatmap of suitable habitat for suitability `>= 0.3` (moderately + highly suitable, R3c / Fig 4 C), matching the revised suitability threshold (changed from `0.7` to `0.3`). Reads the pre-computed `latitudinal_suitability_area_by_class.csv`; does not recompute rasters.

See the repo `README.md` for the article-figure ↔ script mapping and reproducibility notes.
