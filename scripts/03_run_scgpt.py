"""Run scGPT zero-shot reference-mapping annotation on the TME cells.

Writes results/predictions/<dataset>__scgpt.csv (cell_id, predicted_label).

Approach (zero-shot reference mapping, per the user's choice):
  1. Embed the query TME cells with the pretrained scGPT whole-human model.
  2. Embed a LABELLED reference the same way.
  3. Transfer labels by kNN in the shared embedding space.
The query's own ground-truth labels are never used (project guide §7.2) — only
the independent reference provides labels.

IMPORTANT — runs in a SEPARATE environment:
  scGPT pins torch~2.3 / scvi-tools<1.0, incompatible with the main .venv. Use:

    uv venv .venv-scgpt --python 3.10
    uv pip install --python .venv-scgpt scgpt ipython "torch==2.3.1" "torchtext==0.18.0"
    .venv-scgpt/bin/python scripts/03_run_scgpt.py GSE131907

Checkpoint: download the whole-human model (best_model.pt, vocab.json, args.json)
from https://github.com/bowang-lab/scGPT into config scgpt.checkpoint_dir, e.g.
  gdown --folder <whole-human drive link> -O models/scGPT_human

Config keys (config.yaml -> scgpt): checkpoint_dir, reference (labelled .h5ad),
reference_label_column, gene_column, n_neighbors, max_length.

Usage:
    .venv-scgpt/bin/python scripts/03_run_scgpt.py GSE132465
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier

from common import load_config, prediction_file, repo_path


def _embed(adata: sc.AnnData, model_dir: Path, cfg_scgpt: dict, device: str,
           batch_size: int) -> np.ndarray:
    """Return the scGPT cell embedding matrix (n_cells x d) for `adata`.

    Embeds from raw counts (layers["counts"] if present, else X). scGPT does its
    own binning/normalisation internally.
    """
    from scgpt.tasks import embed_data

    a = adata.copy()
    if "counts" in a.layers:
        a.X = a.layers["counts"]
    a.var[cfg_scgpt["gene_column"]] = a.var_names  # scGPT looks up genes by symbol

    embedded = embed_data(
        a,
        model_dir=str(model_dir),
        gene_col=cfg_scgpt["gene_column"],
        max_length=cfg_scgpt["max_length"],
        batch_size=batch_size,
        obs_to_save=None,
        device=device,
        return_new_adata=True,
    )
    # embed_data stores the embedding either in .X or obsm["X_scGPT"].
    if "X_scGPT" in embedded.obsm:
        return np.asarray(embedded.obsm["X_scGPT"])
    return np.asarray(embedded.X)


def main(dataset: str) -> None:
    cfg = load_config()
    s = cfg["scgpt"]

    comp = cfg.get("compute", {})
    device = comp.get("device", "cpu") if comp.get("gpu", False) else "cpu"
    batch_size = comp.get("scgpt_batch_size", 32)  # RTX 3080 10 GB -> ~32; lower on OOM

    model_dir = repo_path(s["checkpoint_dir"])
    if not (model_dir / "best_model.pt").exists():
        sys.exit(f"scGPT checkpoint not found in {model_dir} — download the "
                 "whole-human model (best_model.pt, vocab.json, args.json) first.")
    ref_path = s.get("reference") or ""
    if not ref_path or not Path(repo_path(ref_path)).exists():
        sys.exit("set scgpt.reference in config.yaml to a labelled reference .h5ad")

    tme_path = repo_path(cfg["paths"]["data_raw"], dataset, f"{dataset}_tme.h5ad")
    query = sc.read_h5ad(tme_path)
    ref = sc.read_h5ad(repo_path(ref_path))
    label_col = s["reference_label_column"]
    if label_col not in ref.obs:
        sys.exit(f"reference has no obs column {label_col!r}")

    print(f"embedding reference ({ref.n_obs} cells)...")
    ref_emb = _embed(ref, model_dir, s, device, batch_size)
    print(f"embedding query ({query.n_obs} cells)...")
    qry_emb = _embed(query, model_dir, s, device, batch_size)

    # kNN label transfer in the shared embedding space.
    knn = KNeighborsClassifier(n_neighbors=s["n_neighbors"],
                               n_jobs=comp.get("n_cores", 1))
    knn.fit(ref_emb, ref.obs[label_col].astype(str).values)
    labels = knn.predict(qry_emb)

    out = prediction_file(cfg, dataset, "scgpt")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": query.obs_names, "predicted_label": labels}).to_csv(
        out, index=False
    )
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
