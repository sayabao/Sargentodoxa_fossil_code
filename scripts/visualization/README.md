# Visualization Scripts

Use this directory for new paper-figure scripts.

Recommended pattern:

- Read existing results from `output/`.
- Keep shared figure style settings in one small config script if needed.
- Write paper-ready outputs to `output/figures_paper/`.
- Avoid changing the numbered modeling scripts unless the figure requires a real workflow fix.

Priority inputs:

- `output/08_projection/`
- `output/09_paleoclimate_teyf/`
- `output/10_visualization_teyf/`
- `data/fossil_record/FossilRecord-paleocor-v2.0.csv`
- `output/03_variable_select/`
- `output/06_evaluation/`
- `output/07_ensemble/`

Current scripts:

- `00_paper_occurrence_sampling_map.R`: generates M1 occurrence and pseudo-absence sampling map from the final modeling point table `output/04_formatted/all_points.csv`.
- `01_paper_model_and_trend_figures.R`: generates M2 variable filtering/VIF, M3 model evaluation, M4 variable importance, R3 suitability trend, and the current R3b high-suitability (`0.7-1`) latitudinal area heatmap from 69 Ma to present. Earlier R3b curve variants remain as exploratory outputs.
- `02_paper_suitability_maps.R`: generates R1 present-day suitability, R2 key-period TEYF map panels, and R2b all-period panels from 69 Ma to present using the shared 0-0.1 / 0.1-0.3 / 0.3-0.7 / 0.7-1 suitability classes.
- `03_paper_fossil_comparison_figures.R`: generates R4 fossil-projection panels and R5 fossil-point suitability extraction plots/tables using the same suitability classes.
- `04_server_important_variable_paleoclimate_background.R`: server-side R7 script for summarizing important Bio variables within high-suitability (`>= 0.7`) cells through time. It is RStudio-friendly: edit the settings block at the top and click Source. It requires the original TEYF Bio nc/tif files, e.g. `/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min/`.
- `05_paper_fossil_geographic_shift.R`: generates R6 modern-coordinate vs paleo-coordinate fossil geographic shift map and its plotting table.
- `06_redraw_r7_paleoclimate_background_heatmap.R`: redraws the synchronized R7 heatmap locally from the R7 CSV table, ordering variables by M4 importance and using a taller 4.8 inch output for M4 pairing.
