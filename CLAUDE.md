# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Repo Is

A single-cell RNA-seq **benchmarking research project** (preprint, targeting bioRxiv), not a software product. The scientific question: do general-purpose cell-annotation tools (SingleR, CellTypist, Azimuth, GPTCelltype, scGPT) correctly identify the *functional states* of non-malignant tumour-microenvironment (TME) cells — exhausted T cells, TAM subtypes, CAF subtypes — versus their high accuracy on healthy reference tissue.

**Spec:** `docs/project_guide.md` is the authoritative scientific spec — read it before implementing. It defines datasets, tools, experimental design, metrics, the step-by-step analysis plan, and known pitfalls. (Note: the guide mentions Snakemake; the actual pipeline is deliberately plain Python + R scripts, no workflow manager.)

## Repository Layout

- `docs/project_guide.md` — complete project spec, source of truth for all science decisions.
- `config/config.yaml` — single config both languages read: dataset paths, QC thresholds, the 3 label-granularity levels, label maps, per-tool settings. Edit here, not in scripts.
- `scripts/` — all scripts, both languages, flat (no per-language subfolders). Python: numbered `00`→`05` (download, preprocess, CellTypist, scGPT, metrics, figures) + `common.py`. R: `run_singler.R`, `run_azimuth.R`, `run_gptcelltype.R`, `run_scatomic.R` + `common.R` (mirrors `common.py`).
- `run_all.sh <dataset>` — runs the whole chain in order for one dataset.
- `envs/environment.yml` (conda/Python) + `envs/install_r.R` (Bioconductor/GitHub).
- `data/raw/`, `results/{predictions,metrics,figures}/` — gitignored outputs (`.gitkeep` only).

## Pipeline Flow & Cross-Language Handoff

`01_preprocess.py` writes two `.h5ad` per dataset: `<ds>_processed.h5ad` (all cells) and `<ds>_tme.h5ad` (non-malignant only). **The `.h5ad` is the cross-language contract** — Python tools read it directly; R tools read it via `zellkonverter::readH5AD()`. Every tool, in either language, writes the same prediction schema: `results/predictions/<dataset>__<tool>.csv` with columns `cell_id, predicted_label`. `04_compute_metrics.py` consumes those CSVs language-agnostically. General tools run on `_tme.h5ad`; **scATOMIC alone runs on `_processed.h5ad`** (it needs malignant cells to do its job), then its non-malignant predictions are kept.

## Planned Architecture (from project guide §6–9)

The benchmark is a **multi-tool, multi-dataset, multi-granularity comparison**. Understanding the design requires holding several axes at once:

- **Datasets:** GSE131907 (lung adenocarcinoma, primary ground truth), GSE132465 (colorectal, replication), Zheng68K PBMC (healthy accuracy ceiling). Optional: Werba 2023 PDAC.
- **Tools span two languages** (all under `scripts/`, flat). Python: CellTypist, scGPT, plus all preprocessing/metrics/figures (scanpy). R/Bioconductor: SingleR, Azimuth, GPTCelltype, scATOMIC. Handoff is via AnnData `.h5ad` ↔ SingleCellExperiment/Seurat (zellkonverter).
- **Evaluation runs at 3 label-granularity levels** (project guide §7.3): Level 1 coarse lineage, Level 2 medium subtype, Level 3 fine/functional state. The whole hypothesis is that accuracy is high at L1 and collapses at L3 — so metrics must always be computed per-level.
- **Core output metric:** "TME degradation score" = `Accuracy(Zheng68K) − Accuracy(GSE131907, Level 3)` per tool.

## Critical Conventions (do not violate — these are the experiment's validity)

These come from project guide §7.2 and §12 (Common Pitfalls). Getting them wrong invalidates results, not just code:

- **Never reveal ground-truth labels to the annotation tools.** Each tool runs zero-shot / with default healthy references, as a naive user would. That blindness *is* the experiment.
- **Filter to non-malignant cells only** (using author labels) before running general tools. Cancer cells inflate false positives. Run scATOMIC on all cells to identify malignant cells independently.
- **Use the finest available author label as ground truth**, never the broad label — broad labels make the task trivially easy and hide the effect being measured.
- **Batch-correct (Harmony or scVI) before annotation.** Verify the UMAP does not cluster by `sample_id`. Otherwise annotation follows batch, not biology.
- **Preserve raw counts:** run `adata.raw = adata` *before* `sc.pp.normalize_total()`. Several tools need raw counts.
- **Loose QC thresholds for cancer:** `pct_counts_mt < 20` (not 10), upper `n_genes_by_counts < 6000`. Cancer cells are metabolically active.
- **GPTCelltype is non-deterministic:** run ≥3× per cluster, report modal label + inter-run disagreement rate as a metric. Run it both with `tissuename="lung tumour"` and `tissuename="lung"` to measure the cancer-context effect.
- **"Unknown"/abstain is not the same as wrong** — report unknown rate as a separate metric, never fold it into accuracy.

## Environment Setup

Python env is managed by **uv** (`pyproject.toml` + `uv.lock`, Python 3.10):

```bash
uv sync                       # creates .venv with all Python deps
uv run python scripts/01_preprocess.py GSE131907   # run anything via uv run
```

R env (the four R tools) installs from source — slow, compiles Bioconductor + Seurat:

```bash
Rscript envs/install_r.R      # SingleR, celldex, Azimuth, Seurat, zellkonverter,
                              # + GPTCelltype & scATOMIC from GitHub
```

GPTCelltype needs `OPENAI_API_KEY` in the environment.

**System libs (apt, need sudo) before `install_r.R`:** `libgsl-dev libhdf5-dev libfftw3-dev` — `libgsl-dev` is required to compile the Azimuth dependency chain (DirichletMultinomial → TFBSTools → Signac). scATOMIC pulls `Rmagic` (off CRAN) from GitHub, and Rmagic needs a python backend at runtime: `uv pip install magic-impute`.

**Hardware tuning** lives in the `compute:` block of `config.yaml` (`n_cores`, `gpu`, `device`, `scgpt_batch_size`, `read_chunksize`). Every script reads it: Python sets `sc.settings.n_jobs`; R uses it for `BiocParallel`/`future`/`mc.cores`; the text-matrix reader uses `read_chunksize`. Defaults target a 16-core / RTX 3080 (10 GB) box — change `n_cores` to your physical core count (not threads).

**scGPT is deliberately excluded from the uv env** — it pins `scvi-tools<1.0` and would drag everything back to 2023 versions. Install it in a separate throwaway venv only when running `03_run_scgpt.py` (see the note in `pyproject.toml`).

## Citation Discipline

The project guide cites 18 PubMed papers with verified DOIs and keeps an audit log (§14–15) of every literature search. When adding scientific claims or references, follow the same standard: cite a real DOI, and update the audit log rather than asserting facts unsourced.
