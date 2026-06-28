# 04 — scGPT

**Script:** `scripts/03_run_scgpt.py <dataset>`
**Input:** `<ds>_tme.h5ad` (query) + a labelled reference `.h5ad`
**Output:** `results/predictions/<ds>__scgpt.csv`
**Language:** Python in a **separate** environment (`.venv-scgpt`)

## Methods (paper-ready)

scGPT (Cui et al. 2024, *Nat Methods*), a transformer foundation model pretrained
on tens of millions of cells, was used for zero-shot reference-mapping annotation.
The pretrained whole-human checkpoint embedded both the query TME cells and an
independent labelled reference from raw counts (scGPT performs its own value
binning internally). Labels were transferred from reference to query by k-nearest-
neighbours (k=30) in the shared embedding space. The query's own ground-truth
labels were never used — only the independent reference supplies labels.

## Why a separate environment

scGPT pins `torch~2.3` and `scvi-tools<1.0`, incompatible with the main `.venv`
(which uses modern anndata/scvi). It lives in a throwaway env:

```bash
uv venv .venv-scgpt --python 3.10
uv pip install --python .venv-scgpt scgpt ipython "torch==2.3.1" "torchtext==0.18.0"
.venv-scgpt/bin/python scripts/03_run_scgpt.py GSE131907
```

## Parameters (config `scgpt`)

| Key | Value | Meaning |
|-----|-------|---------|
| `checkpoint_dir` | `models/scGPT_human` | holds `best_model.pt`, `vocab.json`, `args.json` |
| `reference` | `""` **(REQUIRED, unset)** | path to a labelled reference `.h5ad` |
| `reference_label_column` | `cell_type` | obs column with reference labels |
| `gene_column` | `gene_name` | var column with gene symbols |
| `n_neighbors` | 30 | kNN for label transfer |
| `max_length` | 1200 | scGPT input gene budget per cell |
| `compute.scgpt_batch_size` | 32 | fits RTX 3080 10 GB; lower on CUDA OOM |

## Environment state (2026-06-28)

| Component | State |
|-----------|-------|
| `.venv-scgpt` | ✅ exists, `scgpt` importable |
| torch / CUDA | ✅ torch 2.3.1+cu121, CUDA available |
| flash-attn | ⚠️ not installed — warns only, runs without |
| checkpoint `models/scGPT_human/` | ❌ missing — no `models/` dir |
| `scgpt.reference` | ❌ empty in config |

The script hard-exits if `best_model.pt` is absent or `reference` is unset/missing.

## To run — two blockers

1. **Download the whole-human checkpoint** into `models/scGPT_human/`
   (`best_model.pt`, `vocab.json`, `args.json`) from
   <https://github.com/bowang-lab/scGPT>, e.g. `gdown --folder <link> -O models/scGPT_human`.
2. **Set `scgpt.reference`** to a labelled reference `.h5ad` with a `cell_type` obs
   column and gene symbols in `var_names`.

> **Design decision needed:** which reference? For a fair zero-shot test it should
> be an independent healthy/atlas reference (matching the protocol used for SingleR
> and Azimuth), not held-out cells from the same tumour dataset. Document the choice
> — it materially affects the result.

## Status

⬜ Not run — checkpoint + reference both missing.

## Paper notes

- Cite scGPT: Cui et al. 2024, *Nat Methods*
  [10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0).
- Record checkpoint identity (whole-human, training date/version) and the exact
  reference dataset + label column for reproducibility.
