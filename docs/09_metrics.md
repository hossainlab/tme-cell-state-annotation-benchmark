# 09 — Metrics

**Script:** `scripts/04_compute_metrics.py <dataset>`
**Input:** all `results/predictions/<ds>__*.csv` + `<ds>_processed.h5ad` (ground truth)
**Output:** `results/metrics/<ds>_metrics.csv`
**Language:** Python

## Methods (paper-ready)

For each tool, predicted labels were aligned to the author ground truth on cell
barcode. Both predictions and ground truth were collapsed to each of three
granularity levels (coarse lineage, medium subtype, fine functional state) via
hand-curated per-dataset label maps. At each level, accuracy and macro-averaged
F1 were computed over cells with a mappable ground-truth label and a non-abstaining
prediction; abstentions ("Unknown"/"Unassigned") were excluded from accuracy and
reported separately as an unknown rate. Cells whose ground truth did not map to a
given level were masked. The per-tool TME degradation score —
`Accuracy(Zheng68K) − Accuracy(GSE131907, fine)` — was derived from the resulting
tables.

## Output schema

`results/metrics/<ds>_metrics.csv` — one row per (tool × level):

| column | meaning |
|--------|---------|
| `tool` | tool name (from `cfg.tools.python + cfg.tools.r`) |
| `level` | `coarse` / `medium` / `fine` |
| `accuracy` | `accuracy_score` on scorable cells |
| `f1_macro` | `f1_score(average="macro")` |
| `unknown_rate` | fraction of predictions that abstained |
| `n_cells` | number of scored cells |

Abstention tokens: `{unknown, unassigned, nan, none, ""}` (case-insensitive).

## Conventions enforced (pitfall #5)

- Unknown/abstain is **never** scored as wrong — masked out and reported as
  `unknown_rate`.
- Truth that doesn't map to a level is masked (not penalised).
- Tools with no prediction file are skipped with a message.

## ⚠️ Critical blocker — the label maps are empty stubs

```python
def build_level_mapping(dataset: str, level: str) -> dict:
    return {}          # <-- every truth maps to None → masked → n_cells == 0
```

`map_to_level` does `mapping.get(label)`; with an empty dict **every** label →
`None`, so every row reports `n_cells=0` and `NaN` accuracy. The pipeline runs
end-to-end but produces no scores until these maps are curated.

**To produce real metrics:** build, per dataset and per level, a dict
`raw_author_label → level_label` from the confirmed `Cell_subtype` vocabulary,
covering both ground-truth labels AND every tool's output vocabulary (CellTypist's
98 immune types, SingleR HPCA/Blueprint fine labels, Azimuth `ann_level_3`,
scATOMIC pan-cancer labels). Target levels are defined in `config.label_levels`
(`coarse`/`medium`/`fine`).

> This is the single highest-value remaining task: without it there are no results.
> Consider externalising the maps to `config.yaml` or a `label_maps/<ds>.yaml` so
> they are version-controlled and citable, rather than hard-coded in the function.

## Other metrics in the spec not yet in code (guide §8)

| Metric | Status |
|--------|--------|
| Accuracy, F1-macro, unknown rate | ✅ implemented |
| TME degradation score | ⬜ computed downstream (figures), not here |
| Confusion matrix (L3) | ⬜ not implemented (needed for centrepiece figure) |
| Functional-state recall (Tex, Tpex) | ⬜ not implemented |
| "What does Tex get called?" | ⬜ not implemented |
| Wilcoxon + Bonferroni across tools | ⬜ not implemented |

## Run

```bash
uv run python scripts/04_compute_metrics.py GSE131907   # after curating label maps
```

## Status

⬜ Not run, and would yield empty results — `build_level_mapping` returns `{}`.
Note metrics score `cfg.tools` only; `singler_blueprint` is written but not listed.

## Paper notes

- Report accuracy + F1 + unknown rate per level, per tool, per dataset.
- The degradation score and L3 confusion matrices are the core results — implement
  them before drafting Results.
- scikit-learn metrics (Pedregosa 2011, *JMLR*).
