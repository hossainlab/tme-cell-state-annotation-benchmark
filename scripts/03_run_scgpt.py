"""Run scGPT zero-shot reference-mapping annotation on the TME cells.

Writes results/predictions/<dataset>__scgpt.csv (cell_id, predicted_label).

scGPT installs separately and needs a GPU + a downloaded pretrained checkpoint
(whole-human model). This is a stub of the call shape; fill in once the model
weights are in place. Ground-truth labels are never revealed (§7.2).

Usage:
    python scripts/03_run_scgpt.py GSE131907
"""
from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc

from common import load_config, prediction_file, repo_path


def main(dataset: str) -> None:
    cfg = load_config()
    tme_path = repo_path(cfg["paths"]["data_raw"], dataset, f"{dataset}_tme.h5ad")
    adata = sc.read_h5ad(tme_path)

    comp = cfg.get("compute", {})
    device = comp.get("device", "cpu") if comp.get("gpu", False) else "cpu"
    batch_size = comp.get("scgpt_batch_size", 32)  # RTX 3080 10 GB -> ~32; lower on OOM

    # TODO: implement scGPT zero-shot annotation on `device` with `batch_size`.
    #   from scgpt.tasks import embed_data / reference_mapping
    #   load pretrained whole-human checkpoint, map query -> reference labels.
    raise NotImplementedError(
        f"scGPT zero-shot annotation not yet implemented (device={device}, "
        f"batch_size={batch_size}) — needs pretrained checkpoint. "
        "See https://github.com/bowang-lab/scGPT"
    )

    labels = ...  # noqa: F841  (predicted labels aligned to adata.obs_names)
    out = prediction_file(cfg, dataset, "scgpt")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": adata.obs_names, "predicted_label": labels}).to_csv(
        out, index=False
    )
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
