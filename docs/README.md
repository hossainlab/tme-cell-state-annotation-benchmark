# Docs — Step-by-Step Methods

Methods documentation for the TME cell-state annotation benchmark, one file per
pipeline step. Written to map directly onto Methods subsections of the paper:
each file has a paper-ready methods paragraph, the exact script + parameters,
inputs/outputs, current run status, and citation/TODO notes.

## Read in order

| # | File | Step |
|---|------|------|
| 00 | [00_pipeline_overview.md](00_pipeline_overview.md) | Pipeline map, datasets, conventions, status table |
| 01 | [01_data_acquisition.md](01_data_acquisition.md) | Datasets, sources, layout |
| 02 | [02_preprocessing.md](02_preprocessing.md) | QC, normalise, scVI, TME subset |
| 03 | [03_celltypist.md](03_celltypist.md) | CellTypist (Python) |
| 04 | [04_scgpt.md](04_scgpt.md) | scGPT (Python, sep env) |
| 05 | [05_singler.md](05_singler.md) | SingleR (R) |
| 06 | [06_azimuth.md](06_azimuth.md) | Azimuth (R) |
| 07 | [07_gptcelltype.md](07_gptcelltype.md) | GPTCelltype (R + GPT-4) |
| 08 | [08_scatomic.md](08_scatomic.md) | scATOMIC (R, all cells) |
| 09 | [09_metrics.md](09_metrics.md) | Scoring, granularity levels |
| 10 | [10_figures.md](10_figures.md) | Benchmark figures |

## Companion docs

- [project_guide.md](project_guide.md) — authoritative scientific spec (datasets,
  design, bibliography with verified DOIs, audit log).
- [data_download.md](data_download.md) — practical download how-to.

## Status at a glance (2026-06-28)

✅ Done: data acquisition, preprocessing, CellTypist — all 3 datasets.
⬜ Pending: scGPT (checkpoint+reference), SingleR (ready), Azimuth (pkg+ref),
GPTCelltype (API key), scATOMIC (pkg), metrics (label-map stub), figures (stub).

**Two stubs gate all results:** the per-level label maps in
`04_compute_metrics.py::build_level_mapping` (returns `{}`) and the figure plotting
code. See [09](09_metrics.md) and [10](10_figures.md).
