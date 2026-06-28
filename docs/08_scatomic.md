# 08 — scATOMIC

**Script:** `scripts/run_scatomic.R <dataset>`
**Input:** `<ds>_processed.h5ad` (**all** cells, incl. malignant) → raw counts
**Output:** `results/predictions/<ds>__scatomic.csv`
**Language:** R

## Methods (paper-ready)

scATOMIC (Nofech-Mozes et al. 2023, *Nat Commun*), a cancer-specialised hierarchical
classifier, was included as a domain-specific comparator. Unlike the general tools,
scATOMIC was run on the **full processed matrix including malignant cells**, because
distinguishing malignant from non-malignant cells is part of its task. It was run on
raw counts (`run_scATOMIC`), and a summary classification produced per cell
(`create_summary_matrix`). For the benchmark comparison, only its predictions on the
cells the authors labelled non-malignant were retained and scored against the same
ground truth as the general tools.

## Why it differs from every other tool

| | General tools | scATOMIC |
|---|---|---|
| Input object | `<ds>_tme.h5ad` (non-malignant) | `<ds>_processed.h5ad` (all cells) |
| Why | cancer cells inflate false positives | must see cancer to do its job |
| Scored cells | all (non-malignant) | non-malignant subset only |

This is loaded via `load_sce(cfg, dataset, "processed")` (not `load_tme_sce`).

## Parameters

| Item | Value | Notes |
|------|-------|-------|
| input assay | `counts` (raw) | `as.matrix(assay(sce,"counts"))` — scATOMIC needs raw |
| parallel | `mc.cores = n_cores` | 16 |
| output column | `scATOMIC_pred` | from `create_summary_matrix` |

## Blockers

1. **R package `scATOMIC` is not installed** (verified missing). Install from GitHub
   via `envs/install_r.R`. It pulls `Rmagic` (off CRAN, from GitHub), which needs a
   Python backend at runtime: `uv pip install magic-impute`.
2. `as.matrix` densifies the full count matrix — memory-heavy on the 200k-cell lung
   dataset. Watch RAM; this is the heaviest R step.

## Run

```bash
Rscript scripts/run_scatomic.R GSE131907   # full matrix, dense — RAM-heavy
```

## Status

⬜ Not run — `scATOMIC` package missing (+ `magic-impute` backend needed).

## Paper notes

- Cite scATOMIC: Nofech-Mozes et al. 2023, *Nat Commun*
  [10.1038/s41467-023-37353-8](https://doi.org/10.1038/s41467-023-37353-8).
- scATOMIC is the "specialist baseline" — the interesting question is whether a
  cancer-aware tool beats the general ones on non-malignant functional states, or
  whether even it struggles at L3.
- Note that scATOMIC's pan-cancer label vocabulary differs from the author labels;
  the level-mapping (step 09) must map its outputs too.
