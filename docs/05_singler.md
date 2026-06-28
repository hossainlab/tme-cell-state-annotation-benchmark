# 05 — SingleR

**Script:** `scripts/run_singler.R <dataset>`
**Input:** `<ds>_tme.h5ad` → SingleCellExperiment (logcounts assay)
**Output:** `results/predictions/<ds>__singler.csv` and `<ds>__singler_blueprint.csv`
**Language:** R / Bioconductor

## Methods (paper-ready)

SingleR (Aran et al. 2019, *Nat Immunol*) was run on the non-malignant TME cells
using two independent healthy reference atlases distributed via celldex: the Human
Primary Cell Atlas (`HumanPrimaryCellAtlasData`) and Blueprint/ENCODE
(`BlueprintEncodeData`). Each reference's fine labels (`label.fine`) were used for
correlation-based per-cell assignment against the log-normalised expression matrix.
Because both references are healthy, any systematic mislabelling of tumour-specific
functional states is the effect under measurement. Ground-truth labels were never
used. Two prediction files are produced, one per reference, and treated as two
tool variants in the benchmark.

## Parameters

| Item | Value | Source |
|------|-------|--------|
| references | HPCA + Blueprint/ENCODE | celldex |
| labels | `label.fine` | finest available |
| test assay | `logcounts` | aliased from `X` by `common.R` |
| threads | `compute.n_cores` (16) | `num.threads` |

## Cross-language detail

`common.R::load_sce` reads the `.h5ad` with zellkonverter and aliases the `X`
(lognorm) assay to `logcounts`, which `SingleR(assay.type.test="logcounts")`
requires. No data copy.

## Run

```bash
Rscript scripts/run_singler.R GSE131907
Rscript scripts/run_singler.R GSE132465
Rscript scripts/run_singler.R Zheng68K
```

## Status

⬜ Not run. **Fully ready** — `SingleR` (2.8.0) and `celldex` (1.16.0) are
installed; no missing dependencies. This is the only tool with no blocker.

## Paper notes

- Cite SingleR: Aran et al. 2019, *Nat Immunol*
  [10.1038/s41590-018-0276-y](https://doi.org/10.1038/s41590-018-0276-y).
- Record celldex reference versions (Bioconductor release).
- Two references = two columns in the results; decide whether to report both or the
  better-performing one as "SingleR" in the main figure.
- `04_compute_metrics.py` iterates `cfg$tools` = `singler` only — the
  `singler_blueprint` variant is written but **not currently scored**. Add it to the
  tool list (or the metrics loop) if you want both in the tables.
