"""Run CellTypist on the non-malignant TME cells, zero-shot.

Writes results/predictions/<dataset>__celltypist.csv  (cell_id, predicted_label).
Ground-truth labels are NEVER passed to the model (project guide §7.2).

Usage:
    python scripts/02_run_celltypist.py GSE131907
"""
from __future__ import annotations

import sys

import celltypist
import pandas as pd
import scanpy as sc
from celltypist import models

from common import load_config, prediction_file, repo_path


def main(dataset: str) -> None:
    cfg = load_config()
    tme_path = repo_path(cfg["paths"]["data_raw"], dataset, f"{dataset}_tme.h5ad")
    adata = sc.read_h5ad(tme_path)

    models.download_models(model=cfg["celltypist"]["model"])
    model = models.Model.load(model=cfg["celltypist"]["model"])

    pred = celltypist.annotate(
        adata, model=model,
        majority_voting=cfg["celltypist"]["majority_voting"],
    )
    col = "majority_voting" if cfg["celltypist"]["majority_voting"] else "predicted_labels"
    labels = pred.predicted_labels[col]

    out = prediction_file(cfg, dataset, "celltypist")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": adata.obs_names, "predicted_label": labels.values}).to_csv(
        out, index=False
    )
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
