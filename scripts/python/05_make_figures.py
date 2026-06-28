"""Generate the benchmark figures from the metrics + prediction tables.

Figures (project guide §10):
  1. Degradation bar chart   — accuracy by level, healthy vs cancer, per tool
  2. Confusion matrices      — what each tool calls an exhausted CD8 T cell (L3)
  3. Per-cell-type F1 heatmap — cell types x tools
  4. Healthy vs tumour paired comparison
  5. Decision-guide table

Output: results/figures/*.png|svg

Stub: wire each figure to results/metrics/ + results/predictions/ once those
tables are populated. Kept thin so the structure is clear.

Usage:
    python scripts/python/05_make_figures.py
"""
from __future__ import annotations

import pandas as pd

from common import load_config, repo_path


def fig1_degradation(cfg: dict) -> None:
    """Bar chart of accuracy at coarse/medium/fine, healthy vs cancer.

    TME degradation score = Accuracy(Zheng68K) - Accuracy(GSE131907, fine).
    """
    metrics_dir = repo_path(cfg["paths"]["metrics"])
    frames = [pd.read_csv(p) for p in metrics_dir.glob("*_metrics.csv")]
    if not frames:
        print("no metrics yet — run 04_compute_metrics.py first")
        return
    # TODO: pivot tool x level, plot healthy vs cancer side by side, save PNG.


def main() -> None:
    cfg = load_config()
    repo_path(cfg["paths"]["figures"]).mkdir(parents=True, exist_ok=True)
    fig1_degradation(cfg)
    # TODO: fig2_confusion, fig3_f1_heatmap, fig4_paired, fig5_decision_table


if __name__ == "__main__":
    main()
