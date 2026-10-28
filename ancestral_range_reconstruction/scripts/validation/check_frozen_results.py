"""Check packaged data and recompute model weights; no inference or frozen writes."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]


def rows(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


sequences = {}
name = None
for line in (ROOT / "data/phylogeny/concatenated.fasta").read_text().splitlines():
    if line.startswith(">"):
        name = line[1:].strip()
        sequences[name] = ""
    elif name:
        sequences[name] += line.strip()
assert len(sequences) == 42 and {len(s) for s in sequences.values()} == {6399}
sources = rows("data/phylogeny/sequence_sources.csv")
assert len(sources) == 42 and sum(int(r["retained_locus_count"]) for r in sources) == 199
assert sum(r["analysis_group"] == "ingroup" for r in sources) == 23

models = rows("data/ancestral_ranges/frozen_fits/model_summary.csv")
assert len(models) == 6
minimum = min(float(r["AIC"]) for r in models)
unscaled = [math.exp(-0.5 * (float(r["AIC"]) - minimum)) for r in models]
weights = {r["model"]: value / sum(unscaled) for r, value in zip(models, unscaled)}
for r in models:
    assert abs(float(r["AIC"]) - (2 * int(r["k"]) - 2 * float(r["lnL"]))) < 1e-8
    assert abs(float(r["AIC_weight"]) - weights[r["model"]]) < 1e-8
best = min(models, key=lambda r: float(r["AIC"]))
assert best["model"] == "BAYAREALIKE+J" and round(float(best["AIC_weight"]), 4) == 0.3682

ages = rows("output/ancestral_ranges/frozen/node_age_intervals.csv")
n01 = next(r for r in ages if r["node_label"] == "N01")
assert len(n01["clade_key"].split("|")) == 21
for column, expected in [("sampled_age_Ma", 68.91), ("posterior_median_Ma", 73.25),
                         ("HPD95_lower_Ma", 66.00), ("HPD95_upper_Ma", 92.47)]:
    assert round(float(n01[column]), 2) == expected
probabilities = rows("data/ancestral_ranges/frozen_fits/node_probabilities.csv")
weighted = {}
for r in probabilities:
    if r["clade_key"] == n01["clade_key"]:
        state = r["area_state"]
        weighted[state] = weighted.get(state, 0) + weights[r["model"]] * float(r["probability"])
assert abs(sum(weighted.values()) - 1) < 1e-8
assert round(weighted["EA+NAM+EU"], 3) == 0.509
assert round(weighted["EA+SA+NAM+EU"], 3) == 0.224
fossils = rows("data/ancestral_ranges/table_S6.csv")
observations = [r for r in fossils if re.fullmatch(r"O\d+", r["Observation"])]
assert len(observations) == 17 and len({r["Observation"] for r in observations}) == 14

for repeat in [1, 2]:
    original = ROOT / f"data/dating/original_continuation_xml/Walbeck_minimum_posterior_repeat{repeat}.xml"
    fresh = ROOT / f"data/dating/fresh_run_xml/Walbeck_repeat{repeat}_18m.xml"
    assert ET.parse(original).getroot().find("run").get("chainLength") == "2000000"
    assert ET.parse(fresh).getroot().find("run").get("chainLength") == "18000000"

for r in rows("docs/provenance.csv"):
    path = ROOT / r["destination"]
    assert path.is_file()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == r["sha256"], r["destination"]

legacy = ROOT.parent / "niche_modelling"
legacy_checks = "not_present_in_staging"
if legacy.is_dir():
    for filename in ["R7_important_variable_paleoclimate_background.csv",
                     "R7_important_variable_background_contrast.csv"]:
        with (legacy / "output/figures_paper/tables" / filename).open(encoding="utf-8-sig", newline="") as f:
            r7 = list(csv.DictReader(f))
        assert len(r7) == 60 and {float(r["threshold"]) for r in r7} == {0.3}
    with (legacy / "output/figures_paper/tables/R8_fossil_significance_tests.csv").open(encoding="utf-8-sig", newline="") as f:
        tests = list(csv.DictReader(f))
    nh = [r for r in tests if r["region"] == "NH" and r["test"] in ["mean percentile rank", "mean suitability"]]
    assert len(nh) == 2 and all(float(r["p_value"]) == 0.0001 for r in nh)
    legacy_checks = "R7_threshold_0.3_and_R8_NH_Pvalues_verified"

result = {"status": "passed", "accessions": 42, "regions_retained": 199, "sites": 6399,
          "fossil_records": 17, "historical_observations": 14,
          "best_model": best["model"], "AIC_weight": weights[best["model"]],
          "N01_posterior_median_Ma": float(n01["posterior_median_Ma"]),
          "N01_representative_tree_age_Ma": float(n01["sampled_age_Ma"]),
          "N01_EA_NA_EU_probability": weighted["EA+NAM+EU"],
          "new_inference_performed": False,
          "raw_chain_HPD_recalculated": False}
result["legacy_checks"] = legacy_checks
print(json.dumps(result, indent=2))
