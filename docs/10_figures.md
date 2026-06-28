# 10 — Figures

**Script:** `scripts/05_make_figures.py`
**Input:** `results/metrics/*_metrics.csv` + `results/predictions/*`
**Output:** `results/figures/*.png|svg`
**Language:** Python

## Methods (paper-ready)

Benchmark figures were generated from the per-tool metric tables and prediction
files. The central result is the degradation contrast: annotation accuracy at
coarse, medium and fine granularity, compared between healthy tissue (Zheng68K)
and tumour tissue (GSE131907). Confusion matrices at the fine level visualise the
specific failure mode — what each tool calls a tumour-specific functional state
such as an exhausted CD8 T cell. A per-cell-type F1 heatmap and a healthy-vs-tumour
paired comparison summarise across tools, and a decision-guide table translates the
results into practical recommendations.

## Planned figures (guide §10)

| # | Figure | Source | Status |
|---|--------|--------|--------|
| 1 | Degradation bar chart — accuracy by level, healthy vs cancer, per tool | metrics | ⬜ stub (`fig1_degradation` loads CSVs, TODO plot) |
| 2 | Confusion matrices — what each tool calls an exhausted CD8 (L3) | predictions | ⬜ not implemented |
| 3 | Per-cell-type F1 heatmap (cell types × tools) | metrics | ⬜ not implemented |
| 4 | Healthy vs tumour paired comparison | metrics | ⬜ not implemented |
| 5 | Practical decision-guide table | metrics | ⬜ not implemented |
| S | Colorectal (GSE132465) replication | metrics | ⬜ not implemented |

The TME degradation score (`Accuracy(Zheng68K) − Accuracy(GSE131907, fine)`) is
computed and plotted here (figure 1), not in step 09.

## Current code state

`05_make_figures.py` is a **thin stub**: `fig1_degradation` globs
`results/metrics/*_metrics.csv` and returns early if none exist; all plotting is
`TODO`. `fig2`–`fig5` are placeholders in `main`.

## Dependency chain

Figures depend on step 09 producing non-empty metrics, which depends on curating
the label maps. **Nothing meaningful renders until [09_metrics](09_metrics.md)'s
`build_level_mapping` stub is filled.** Confusion matrices (fig 2) also need a
confusion-matrix routine that 09 does not yet compute.

## Run

```bash
uv run python scripts/05_make_figures.py   # after metrics exist + plotting implemented
```

## Status

⬜ Not run — plotting unimplemented; upstream metrics empty.

## Paper notes

- Figure 2 (L3 confusion, exhausted T cells) is the paper's centrepiece — prioritise.
- Keep a consistent tool colour/order across all panels.
- Save vector (SVG/PDF) for the manuscript; PNG for quick review.
- Suggested libs: matplotlib + seaborn (`seaborn` / `scientific-visualization`).
