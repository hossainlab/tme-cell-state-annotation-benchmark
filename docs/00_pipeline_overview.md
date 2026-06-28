# 00 — Pipeline Overview

Methods documentation, one file per pipeline step, written so each maps directly
onto a Methods subsection of the paper. `docs/project_guide.md` remains the
scientific spec; these files document **what the code actually does** and the
current run status.

## Research question

Do general-purpose cell-annotation tools (SingleR, CellTypist, Azimuth,
GPTCelltype, scGPT) correctly identify the *functional states* of non-malignant
tumour-microenvironment (TME) cells — exhausted T cells, TAM subtypes, CAF
subtypes — versus their high accuracy on healthy reference tissue?

Core summary statistic: **TME degradation score** =
`Accuracy(Zheng68K) − Accuracy(GSE131907, Level 3)`, per tool.

## Step map

| Step | Doc | Script | Lang | Output |
|------|-----|--------|------|--------|
| 1 | [01_data_acquisition](01_data_acquisition.md) | `00_download_data.py` | py | `data/raw/<ds>/…` |
| 2 | [02_preprocessing](02_preprocessing.md) | `01_preprocess.py` | py | `<ds>_processed.h5ad`, `<ds>_tme.h5ad` |
| 3 | [03_celltypist](03_celltypist.md) | `02_run_celltypist.py` | py | `<ds>__celltypist.csv` |
| 4 | [04_scgpt](04_scgpt.md) | `03_run_scgpt.py` | py (sep env) | `<ds>__scgpt.csv` |
| 5 | [05_singler](05_singler.md) | `run_singler.R` | R | `<ds>__singler{,_blueprint}.csv` |
| 6 | [06_azimuth](06_azimuth.md) | `run_azimuth.R` | R | `<ds>__azimuth.csv` |
| 7 | [07_gptcelltype](07_gptcelltype.md) | `run_gptcelltype.R` | R | `<ds>__gptcelltype_{lungtumour,lung}.csv` |
| 8 | [08_scatomic](08_scatomic.md) | `run_scatomic.R` | R | `<ds>__scatomic.csv` |
| 9 | [09_metrics](09_metrics.md) | `04_compute_metrics.py` | py | `<ds>_metrics.csv` |
| 10 | [10_figures](10_figures.md) | `05_make_figures.py` | py | `results/figures/*` |

`run_all.sh <dataset>` chains steps 1→10 for one dataset.

## Datasets (three axes of the benchmark)

| Dataset | Role | Tissue | Notes |
|---------|------|--------|-------|
| GSE131907 | primary ground truth | lung adenocarcinoma | Kim et al. 2020, *Nat Commun*, [10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1) |
| GSE132465 | replication | colorectal | Lee et al. 2020, *Nat Genet*, [10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z) |
| Zheng68K | healthy accuracy ceiling | PBMC | Zheng et al. 2017, *Nat Commun*, [10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049) |

Evaluation runs at **3 label-granularity levels**: L1 coarse lineage, L2 medium
subtype, L3 fine/functional state. Hypothesis: accuracy high at L1, collapses at L3.

## Cross-language contract

`01_preprocess.py` writes two `.h5ad` per dataset. The `.h5ad` is the handoff:
Python tools read it directly; R tools read it via `zellkonverter::readH5AD()`.
zellkonverter names the assays `X` (lognorm from `adata.X`) and `counts` (raw from
`layers["counts"]`); `common.R` aliases `X → logcounts` for SingleR. Every tool
writes the same schema: `results/predictions/<dataset>__<tool>.csv` with columns
`cell_id, predicted_label`. `04_compute_metrics.py` consumes those CSVs
language-agnostically.

## Inviolable conventions (experiment validity)

1. **Never reveal ground-truth labels to a tool** — every tool runs zero-shot / with
   default healthy references.
2. **Filter to non-malignant cells** (author labels) before general tools; scATOMIC
   alone runs on all cells.
3. **Use the finest author label as ground truth** (`Cell_subtype`), never the broad
   lineage.
4. **Batch-correct before annotation** (scVI here); UMAP must not cluster by sample.
5. **Preserve raw counts** in `layers["counts"]` before normalisation.
6. **Loose cancer QC** (`mt < 20%`, `genes < 6000`).
7. **"Unknown"/abstain ≠ wrong** — reported as a separate rate.

## Current status (2026-06-28)

| Step | GSE131907 | GSE132465 | Zheng68K | Blocker |
|------|-----------|-----------|----------|---------|
| 1 Download | ✅ | ✅ | ✅ | — |
| 2 Preprocess | ✅ | ✅ | ✅ | — |
| 3 CellTypist | ✅ | ✅ | ✅ | — |
| 4 scGPT | ⬜ | ⬜ | ⬜ | checkpoint + reference unset |
| 5 SingleR | ⬜ | ⬜ | ⬜ | ready to run |
| 6 Azimuth | ⬜ | ⬜ | ⬜ | R pkg `Azimuth` missing |
| 7 GPTCelltype | ⬜ | ⬜ | ⬜ | `OPENAI_API_KEY` unset |
| 8 scATOMIC | ⬜ | ⬜ | ⬜ | R pkg `scATOMIC` missing |
| 9 Metrics | ⬜ | ⬜ | ⬜ | `build_level_mapping` is a stub `{}` |
| 10 Figures | ⬜ | ⬜ | ⬜ | plotting code is a stub |

Legend: ✅ done · ⬜ not run.

> **Two stubs gate the whole result.** Even with all predictions in hand, metrics
> return `n_cells=0` until per-dataset, per-level label maps are curated
> (`04_compute_metrics.py::build_level_mapping`), and figures are unimplemented.
> See [09_metrics](09_metrics.md) and [10_figures](10_figures.md).
