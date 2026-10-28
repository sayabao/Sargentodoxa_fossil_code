# Molecular phylogeny, dating and ancestral-range reconstruction

This module supports Sections 2.3 and 3.3, Supplementary Methods 2–4, Fig. 10, Fig. S5 and Tables S5–S7.

| Directory | Contents |
| --- | --- |
| data/phylogeny | Sequence sources, five-region alignment, partitions and ML topology |
| data/dating | BEAST model configurations for two chains |
| data/ancestral_ranges | Extant ranges, fossil observations, dated tree and six-model fitted results |
| scripts | Phylogenetic, dating, ancestral-range, plotting and numerical-check scripts |
| output | Tree components, node ages, range probabilities and model weights |

## Software

MAFFT 7.526; IQ-TREE 3.1.4; BEAST 2.7.8 with ORC 1.2.1; BioGeoBEARS 1.1.4. Python plotting dependencies are listed in requirements.txt. R fitting scripts require BioGeoBEARS and its dependencies, including ape, optimx, rexpokit and cladoRcpp.

## Numerical checks and figures

Run from the repository root:

```sh
python ancestral_range_reconstruction/scripts/validation/check_frozen_results.py
python ancestral_range_reconstruction/scripts/phylogeny/draw_primary_ml.py
python ancestral_range_reconstruction/scripts/ancestral_ranges/draw_ancestral_probability_tree.py --with-outgroups
```

## Inference

The Python wrappers below print the inference commands; add --execute to run them.

```sh
python ancestral_range_reconstruction/scripts/phylogeny/run_iqtree.py
python ancestral_range_reconstruction/scripts/dating/run_beast.py --repeat 1
python ancestral_range_reconstruction/scripts/dating/run_beast.py --repeat 2
Rscript ancestral_range_reconstruction/scripts/ancestral_ranges/refit_primary_initial.R ancestral_range_reconstruction
Rscript ancestral_range_reconstruction/scripts/ancestral_ranges/refine_primary_fits.R ancestral_range_reconstruction
```

The original continuation XMLs describe the final 2-million-state extension; separate 18-million-state XMLs are provided for complete runs. Node-age summaries are included; full posterior-chain files are not bundled.

The phylogenetic matrix contains 42 accessions and 6,399 positions. Ancestral ranges use 21 sampled family species and 17 fossil records combined into 14 observations. NAM in data files denotes North America (NA in the manuscript). N01 has a posterior median age of 73.25 Ma; its height in the representative tree is 68.91 Ma.
