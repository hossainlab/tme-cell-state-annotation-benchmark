"""Score every tool's predictions against author ground truth, per level.

Reads results/predictions/<dataset>__<tool>.csv for all tools, joins to the
author ground truth in the processed .h5ad, maps both to each granularity
level, and writes a tidy metrics table.

Output: results/metrics/<dataset>_metrics.csv
        columns: tool, level, accuracy, f1_macro, unknown_rate, n_cells

Conventions (project guide §8, §12):
  - "Unknown"/abstain is reported separately, never scored as wrong.
  - Cells whose ground truth does not map to a level are masked out.
  - TME degradation score (vs Zheng68K) is computed in 05/figures, not here.

Usage:
    python scripts/python/04_compute_metrics.py GSE131907
"""
from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score

from common import LEVELS, load_config, map_to_level, prediction_file, repo_path

UNKNOWN_TOKENS = {"unknown", "unassigned", "nan", "none", ""}


def build_level_mapping(dataset: str, level: str) -> dict:
    """Raw author label -> level label.

    TODO: hand-curate one dict per dataset/level from the confirmed author
    labels. Left empty here so the pipeline runs end-to-end; until filled,
    every truth maps to None and is masked (n_cells == 0).
    """
    return {}


def main(dataset: str) -> None:
    cfg = load_config()
    processed = repo_path(cfg["paths"]["data_raw"], dataset, f"{dataset}_processed.h5ad")
    adata = sc.read_h5ad(processed)
    truth_col = cfg["datasets"][dataset]["truth_column"]
    truth_raw = adata.obs[truth_col]

    all_tools = cfg["tools"]["python"] + cfg["tools"]["r"]
    rows = []
    for tool in all_tools:
        pred_path = prediction_file(cfg, dataset, tool)
        if not pred_path.exists():
            print(f"skip {tool}: no predictions at {pred_path}")
            continue
        pred = pd.read_csv(pred_path).set_index("cell_id")["predicted_label"]
        pred = pred.reindex(adata.obs_names)  # align to the same cell order

        for level in LEVELS:
            mapping = build_level_mapping(dataset, level)
            true_lvl = pd.Series(map_to_level(truth_raw, level, mapping),
                                 index=adata.obs_names)
            pred_lvl = pd.Series(map_to_level(pred.fillna("Unknown"), level, mapping),
                                 index=adata.obs_names)

            is_unknown = pred.astype(str).str.lower().isin(UNKNOWN_TOKENS)
            unknown_rate = float(is_unknown.mean())

            # Score only cells with a mappable truth and a non-abstaining pred.
            mask = true_lvl.notna() & pred_lvl.notna() & ~is_unknown.values
            n = int(mask.sum())
            if n == 0:
                rows.append(dict(tool=tool, level=level, accuracy=float("nan"),
                                 f1_macro=float("nan"), unknown_rate=unknown_rate,
                                 n_cells=0))
                continue
            t, p = true_lvl[mask], pred_lvl[mask]
            rows.append(dict(
                tool=tool, level=level,
                accuracy=accuracy_score(t, p),
                f1_macro=f1_score(t, p, average="macro", zero_division=0),
                unknown_rate=unknown_rate,
                n_cells=n,
            ))

    out = repo_path(cfg["paths"]["metrics"]) / f"{dataset}_metrics.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
