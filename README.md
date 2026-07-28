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
| `documents/` | TEYF period metadata, variable-selection methods, reproducibility & figure-mapping notes |
| `output/` | Frozen model/projection outputs that let the figure scripts reproduce every article figure locally |

## Quickstart - reproduce the article component figures

The figure scripts read from the frozen `output/` and `data/` shipped here, so
they run locally without any server data. From the repo root:

```bash
Rscript scripts/visualization/00_paper_occurrence_sampling_map.R    # Fig S1 (M1)
Rscript scripts/visualization/01_paper_model_and_trend_figures.R    # Fig 3 (M4), Fig 4 (R3,R3b), Fig S2 (M2), Fig S3 (M3)
Rscript scripts/visualization/02_paper_suitability_maps.R           # Fig 4 (R2b)
Rscript scripts/visualization/03_paper_fossil_comparison_figures.R  # Fig 5 (R4,R5)
Rscript scripts/visualization/07_redraw_r7_contrast_heatmap_trends.R # Fig 3 (R7 contrast)
Rscript scripts/visualization/09_fossil_suitability_significance.R  # Fig S11 (R8)
```

`Fig S4` is taken from `output/08_projection/projection_modern_map.pdf` (produced
by `scripts/08_projection.R`): panel A = page 1, panel B = page 6
("Presence vs Background: Ensemble (EMwmean)"), with a green dashed line added
manually by the author. Outputs are written under `output/figures_paper/`. See
`documents/FIGURE_MAPPING.md` for the exact article-figure <-> script <->
output mapping.

## Reproducibility boundary (summary)

- **Locally reproducible from frozen outputs:** all `scripts/visualization/`
  figure scripts (and the Presence-vs-Background page of `08_projection.R`).
- **Requires server data (not shipped):** the upstream pipeline `01-07` needs
  the full WorldClim 2.1 bioclim stack; `09_paleoclimate_projection.R` needs the
  TEYF palaeoclimate rasters; `04_server_*.R` (R7) needs the raw TEYF Bio
  variables. These scripts are kept for methodology; their frozen outputs ship
  in `output/` so downstream figures reproduce without re-running them.

See `documents/REPRODUCIBILITY.md` for full detail (server paths, env vars,
seeds, and what was excluded).

## R environment

R 4.5.1. Package versions are recorded in `documents/R_ENVIRONMENT.md`. Key
packages: biomod2, terra, sf, ggplot2, patchwork, tidyterra, usdm. (renv not
enabled; see that file for how to add a strict lockfile.)

## Note on fossil taxonomy

The analysis file `data/fossil_record/FossilRecord-paleocor-v2.1.csv` labels
the Tajimi (Japan) record as *Sargentodoxa cuneata*, consistent with the
taxonomic naming at the time of modelling. The manuscript (v5.3) subsequently
revises this record to *Sargentodoxa tajimiensis* sp. nov. This is a taxonomic
label update only - coordinates, age, and TEYF period assignment are unchanged,
so it does not affect the niche model or any figure.
