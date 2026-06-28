# TME Cell-State Annotation Benchmark

### Can general-purpose annotation tools correctly identify the *functional states* of non-malignant cells in the tumour microenvironment?

![Status](https://img.shields.io/badge/status-in%20development-yellow)
![Python](https://img.shields.io/badge/python-3.10-blue)
![R](https://img.shields.io/badge/R-4.4-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## Overview

Automated cell-type annotation tools report 80–90% accuracy — but almost always on **healthy reference tissue**, and almost always at the **coarse lineage level** (*"T cell", "macrophage", "fibroblast"*). The clinically decisive distinctions in cancer live one level deeper: a *precursor-exhausted* CD8⁺ T cell (Tpex) predicts response to PD-1 blockade, while a *terminally exhausted* one (Tex) does not; an *SPP1⁺* tumour-associated macrophage promotes metastasis, while an *FOLR2⁺* one tracks with better outcomes; an *inflammatory* CAF (iCAF) shapes immune exclusion differently from a *myofibroblastic* CAF (myCAF).

**This project asks a question no existing benchmark answers:** when you hand a standard annotation tool the non-malignant cells of a real tumour, does it recover these *functional states* — or does it quietly collapse them into the healthy-tissue labels it was trained on?

We benchmark **six widely used tools** — SingleR, CellTypist, Azimuth, GPTCelltype (GPT-4), scGPT, and the cancer-specialised scATOMIC — against curated functional-state annotations from pan-cancer single-cell atlases, and quantify exactly where, and by how much, each one degrades.

## Why this matters

- **The gap is real and specific.** The closest prior work, scATOMIC (Nofech-Mozes et al., *Nat Commun* 2023), benchmarks tools at the *malignant vs. non-malignant* boundary and finds general tools score F1 > 0.85 on non-malignant cells. That is precisely the problem: non-malignant cells look fine *in bulk*. No benchmark has tested whether tools recover the functional **sub-states within** those populations — the states that actually drive immunotherapy response.
- **The tooling has outpaced its validation.** LLM-based annotators (GPTCelltype, *Nat Methods* 2024) report strong concordance with manual labels — but their evaluation tissues are predominantly healthy, and the cancer data they do touch is scored at the lineage level. Their behaviour on Tpex-vs-Tex or iCAF-vs-myCAF is simply unknown.
- **The ground truth already exists.** Pan-cancer atlases (Kim et al. *Nat Commun* 2020; Lee et al. *Nat Genet* 2020) ship expert-curated functional-state labels. No new data generation is required — only a rigorous, reproducible evaluation.

## Central hypothesis

> General tools will retain high accuracy at the **lineage level** (Level 1) in cancer tissue, hold up moderately at the **subtype level** (Level 2), and **fail at the functional-state level** (Level 3) — defaulting to canonical healthy-tissue labels. LLM-based tools may partially recover functional states (literature exposure) but at the cost of **run-to-run inconsistency**. Cancer-specialised tools may win at identifying malignant cells without being any better at sub-states *within* the non-malignant compartment.

## Experimental design

Every tool is evaluated across two orthogonal axes:

**Three label-granularity levels** — the heart of the design:

| Level | Granularity | Example labels |
|-------|-------------|----------------|
| **L1 — Coarse** | Lineage | T cell · Myeloid · Fibroblast · Endothelial · B cell |
| **L2 — Medium** | Subtype | CD8 T · CD4 T · NK · Macrophage · Monocyte · DC · CAF |
| **L3 — Fine** | **Functional state** | **Exhausted/Naïve/Effector CD8 · Tpex · Tex · SPP1⁺ TAM · CCL18⁺ TAM · iCAF · myCAF · matCAF** |

**Three datasets** — measuring the *drop*, not just the score:

| Dataset | Role | Tissue |
|---------|------|--------|
| **Zheng68K PBMC** | Healthy accuracy ceiling | Blood |
| **GSE131907** | Primary cancer ground truth | Lung adenocarcinoma (~200k cells) |
| **GSE132465** | Cross-cancer replication | Colorectal (~90k cells) |

The headline metric is the **TME degradation score** = `Accuracy(healthy, Zheng68K) − Accuracy(cancer, GSE131907 @ L3)`, computed per tool. Supporting metrics: per-cell-type F1, abstain/"Unknown" rate (reported separately — honest uncertainty is *not* counted as error), and, for the LLM tools, inter-run disagreement.

## Methodological rigour

The validity of a benchmark lives in its controls. This pipeline enforces them as first-class conventions (see [`CLAUDE.md`](CLAUDE.md) and `docs/project_guide.md` §12):

- **No label leakage.** Tools run zero-shot with default *healthy* references, exactly as a naïve user would invoke them. That blindness *is* the experiment.
- **Non-malignant cells only** for the general tools (malignant cells are filtered using author labels first), so misclassified tumour cells cannot inflate false positives. scATOMIC alone sees the full matrix, since identifying malignant cells is its job.
- **Finest available author label as ground truth** — never the broad label, which makes the task trivially easy and hides the effect being measured.
- **Batch correction (Harmony) before annotation**, verified by checking the UMAP does not cluster by sample.
- **Cancer-appropriate QC** (mtDNA < 20%, not the healthy-tissue 10%) and **raw counts preserved** before normalisation.
- **LLM non-determinism handled explicitly:** GPTCelltype is run 3× per cluster, reporting the modal label *and* the disagreement rate, with and without a cancer-context hint.

## Pipeline

A deliberately lightweight, two-language pipeline — plain scripts, no workflow manager — bridged by a single cross-language contract.

```
config/config.yaml      ── single source of truth: paths, QC, label levels, tool settings
        │
scripts/00_download → 01_preprocess ──► <dataset>_processed.h5ad
                                        <dataset>_tme.h5ad   (non-malignant subset)
        │                                       │
        │   Python tools                        │   R tools (read .h5ad via zellkonverter)
        ├── 02_run_celltypist.py                ├── run_singler.R
        ├── 03_run_scgpt.py                     ├── run_azimuth.R
        │                                       ├── run_gptcelltype.R
        │                                       └── run_scatomic.R
        │                                       │
        └──────────► results/predictions/<dataset>__<tool>.csv  (cell_id, predicted_label)
                                        │
                  04_compute_metrics.py ──► results/metrics/   (per-tool, per-level)
                  05_make_figures.py    ──► results/figures/
```

**The cross-language contract.** Preprocessing emits an AnnData `.h5ad`; R tools read it through `zellkonverter`; *every* tool — Python or R — writes the identical `cell_id, predicted_label` CSV schema. Metrics are therefore computed language-agnostically, and adding a seventh tool means writing one script that emits one CSV.

## Repository structure

```
config/            config.yaml — drives both languages
scripts/           all analysis scripts, flat (Python + R side by side)
envs/              install_r.R — Bioconductor + GitHub R dependencies
docs/              project_guide.md — full scientific protocol, lit review, bibliography
data/raw/          input datasets (gitignored)
results/           predictions · metrics · figures (gitignored)
pyproject.toml     uv-managed Python environment (Python 3.10)
run_all.sh         ordered end-to-end runner
CLAUDE.md          engineering conventions and validity invariants
```

## Getting started

**Python environment** (managed with [uv](https://github.com/astral-sh/uv)):

```bash
uv sync                                    # builds .venv from the lockfile (Python 3.10)
```

**R environment** (the four R-based tools):

```bash
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4   # writable user library
Rscript envs/install_r.R                        # SingleR, celldex, Azimuth, Seurat,
                                                # zellkonverter, GPTCelltype, scATOMIC
export OPENAI_API_KEY=...                        # required by GPTCelltype
```

**Run the benchmark for one dataset:**

```bash
./run_all.sh GSE131907
# → results/metrics/GSE131907_metrics.csv  and  results/figures/
```

> scGPT is installed in a separate environment (it pins an older `scvi-tools`; see `pyproject.toml`).

## Planned deliverables

- **Preprint** (bioRxiv) reporting per-tool, per-level degradation across cancer types.
- **Confusion-matrix atlas** showing *what* each tool calls an exhausted T cell, an SPP1⁺ macrophage, a myCAF.
- **A practical decision guide** — *"which annotation tool should you use for TME work, and when?"* — the most actionable output for the field.
- **Fully reproducible pipeline** (this repository).

## Expected figures

| # | Figure | Question it answers |
|---|--------|--------------------|
| 1 | Degradation bar chart (L1/L2/L3, healthy vs. cancer) | How far does each tool fall, and at which granularity? |
| 2 | Confusion matrices for exhausted T cells | *What* does each tool mistake them for? |
| 3 | Per-cell-type F1 heatmap | Which functional states are hardest? |
| 4 | Healthy-vs-tumour paired comparison | Is degradation specific to the tumour context? |
| 5 | Decision guide | Which tool, for which goal? |

## Documentation

The complete scientific protocol — motivation, systematic literature review with a 12-search audit log, dataset provenance, tool rationale, step-by-step analysis plan, metrics, common pitfalls, and a fully cited bibliography (18 DOI-verified references) — lives in **[`docs/project_guide.md`](docs/project_guide.md)**.

## Author

**Md. Jubayer Hossain** · pre-PhD research project
GitHub: [@hossainlab](https://github.com/hossainlab)

## License

[MIT](LICENSE)
