#!/usr/bin/env bash
# End-to-end benchmark for one dataset. No workflow manager — plain Python + R,
# run in order. Activate the conda env (envs/environment.yml) first; the R tools
# need the packages from envs/install_r.R and OPENAI_API_KEY for GPTCelltype.
#
# Usage: ./run_all.sh GSE131907
set -euo pipefail

DATASET="${1:?usage: ./run_all.sh <dataset>}"
PY=scripts/python
R=scripts/r

echo "== [$DATASET] download =="
python "$PY/00_download_data.py" "$DATASET"

echo "== [$DATASET] preprocess (QC, normalise, Harmony, TME subset) =="
python "$PY/01_preprocess.py" "$DATASET"

echo "== [$DATASET] Python tools =="
python "$PY/02_run_celltypist.py" "$DATASET"
python "$PY/03_run_scgpt.py" "$DATASET" || echo "scgpt skipped (needs checkpoint/GPU)"

echo "== [$DATASET] R tools =="
Rscript "$R/run_singler.R"     "$DATASET"
Rscript "$R/run_azimuth.R"     "$DATASET"
Rscript "$R/run_gptcelltype.R" "$DATASET"
Rscript "$R/run_scatomic.R"    "$DATASET"

echo "== [$DATASET] metrics + figures =="
python "$PY/04_compute_metrics.py" "$DATASET"
python "$PY/05_make_figures.py"

echo "done: results/metrics/${DATASET}_metrics.csv, results/figures/"
