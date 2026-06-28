#!/usr/bin/env bash
# End-to-end benchmark for one dataset. No workflow manager — plain Python + R,
# run in order. Activate the conda env (envs/environment.yml) first; the R tools
# need the packages from envs/install_r.R and OPENAI_API_KEY for GPTCelltype.
#
# Usage: ./run_all.sh GSE131907
set -euo pipefail

DATASET="${1:?usage: ./run_all.sh <dataset>}"
S=scripts

echo "== [$DATASET] download =="
python "$S/00_download_data.py" "$DATASET"

echo "== [$DATASET] preprocess (QC, normalise, Harmony, TME subset) =="
python "$S/01_preprocess.py" "$DATASET"

echo "== [$DATASET] Python tools =="
python "$S/02_run_celltypist.py" "$DATASET"
python "$S/03_run_scgpt.py" "$DATASET" || echo "scgpt skipped (needs checkpoint/GPU)"

echo "== [$DATASET] R tools =="
Rscript "$S/run_singler.R"     "$DATASET"
Rscript "$S/run_azimuth.R"     "$DATASET"
Rscript "$S/run_gptcelltype.R" "$DATASET"
Rscript "$S/run_scatomic.R"    "$DATASET"

echo "== [$DATASET] metrics + figures =="
python "$S/04_compute_metrics.py" "$DATASET"
python "$S/05_make_figures.py"

echo "done: results/metrics/${DATASET}_metrics.csv, results/figures/"
