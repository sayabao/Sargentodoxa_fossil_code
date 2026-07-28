# Sargentodoxa fossil SDM - code & reproducibility repository

This repository contains the species distribution model (SDM) code, frozen
intermediate outputs, and figure-generating scripts for the palaeoecological
analysis of *Sargentodoxa* reported in the accompanying manuscript (v5.3,
submitted to *Palaeogeography, Palaeoclimatology, Palaeoecology*).

It is a **curated, de-servered snapshot** of the analysis workspace: only the
scripts, data, and outputs that the manuscript actually uses are included, and
all hard-coded server paths have been replaced with relative paths / environment
variables. The repository ships **component figures**; the author manually
assembles the final article figures in Photoshop.

## What is here

| Path | Contents |
|---|---|
| `scripts/01_…_10_*.R` | Numbered biomod2 SDM pipeline (data prep -> ensemble -> modern/paleo projection -> visualization) |
| `scripts/config.R` | Global config; paths de-servered (local data via relative paths, climate rasters via env vars) |
| `scripts/get_paleocoordinates.R` | Generates palaeo-coordinates for fossil records via the GPlates PALEOMAP service |
| `scripts/visualization/` | Paper-figure component scripts (M1-M4, R2b, R3, R3b, R4, R5, R7, R8) |
| `data/occ_data_cvh/` | Occurrence data (CVH) |
| `data/fossil_record/` | Fossil record tables incl. palaeo-coordinates (`FossilRecord-paleocor-v2.1.csv`) |
| `data/teyf_simulations.csv` | TEYF period metadata (Simulation Name, Ma) used to match fossils/projections to periods |
| `output/` | Frozen model/projection outputs that let the figure scripts reproduce every article figure locally |

## Quickstart - reproduce the article component figures

The figure scripts read from the frozen `output/` and `data/` shipped here, so
they run locally without any server data. From the repo root:

```bash
Rscript scripts/visualization/00_paper_occurrence_sampling_map.R    # Fig S1 (M1)
Rscript scripts/visualization/01_paper_model_and_trend_figures.R    # Fig 3 (M4), Fig 4 (R3,R3b), Fig S2 (M2), Fig S3 (M3)
Rscript scripts/visualization/02_paper_suitability_maps.R           # Fig 4 (R2b)
Rscript scripts/visualization/03_paper_fossil_comparison_figures.R  # Fig 5 (R4,R5)
Rscript scripts/visualization/05_redraw_r7_contrast_heatmap_trends.R # Fig 3 (R7 contrast)
Rscript scripts/visualization/06_fossil_suitability_significance.R  # Fig S11 (R8)
```

`Fig S4` is taken from `output/08_projection/projection_modern_map.pdf` (produced
by `scripts/08_projection.R`): panel A = page 1, panel B = page 6
("Presence vs Background: Ensemble (EMwmean)"), with a green dashed line added
manually by the author. Outputs are written under `output/figures_paper/`.

## Reproducibility

**Locally reproducible from frozen outputs:** all `scripts/visualization/`
figure scripts (and the Presence-vs-Background page of `08_projection.R`).

**Requires server data (not shipped):** the upstream pipeline `01-07` needs the
full WorldClim 2.1 bioclim stack; `09_paleoclimate_projection.R` needs the TEYF
palaeoclimate rasters; `04_server_*.R` (R7) needs the raw TEYF Bio variables.
These scripts are kept for methodology; their frozen outputs ship in `output/`
so downstream figures reproduce without re-running them.

**Random seeds:** spatial thinning / CV use `thinning_seed = 42` and
`seed_val = 123` (see `scripts/config.R`); model training, ensembling, and
projection are bound by these, so results are reproducible.

## R environment

R 4.5.1 (ucrt, Windows). versions are recorded here.
Key packages: biomod2 4.3.4.5, terra 1.9.1, raster 3.6.32, sf 1.1.0, dplyr 1.2.0,
readr 2.1.6, ggplot2 4.0.1, patchwork 1.3.2, tidyterra 1.1.0, viridis 0.6.5,
corrplot 0.95, usdm 2.1.7.
