# 02 — Preprocessing

**Script:** `scripts/01_preprocess.py <dataset>`
**Input:** raw matrix + author annotation (step 01)
**Output:** `data/raw/<ds>/<ds>_processed.h5ad` (all cells) and `<ds>_tme.h5ad`
(non-malignant only)

## Methods (paper-ready)

Each dataset was loaded into an AnnData object and author annotations joined on the
cell barcode. Quality control used deliberately permissive thresholds appropriate
for metabolically active cancer cells: cells with >200 and <6,000 detected genes
and <20% mitochondrial counts were retained. Raw UMI counts were preserved in
`layers["counts"]` before any transformation. Counts were library-size normalised
to 10,000 counts per cell and log1p-transformed; the log-normalised matrix was
stored in `adata.raw` prior to feature selection. The 2,000 most highly variable
genes were identified with the Seurat v3 method on raw counts (batch-aware when
>1 sample). For multi-sample datasets, batch effects were corrected with scVI
(30 latent dimensions, 2 layers) trained on the HVG counts over the sample key; a
neighbour graph and UMAP were computed from the scVI latent space. Single-sample
data (Zheng68K) used PCA (50 PCs) with no integration. After embedding, cells were
split into the full processed object and a non-malignant subset by removing the
dataset's malignant/ambiguous lineages (e.g. `Epithelial cells`, `Undetermined`)
using author labels only.

## Two input readers

| Format | Reader | Datasets | Notes |
|--------|--------|----------|-------|
| `geo_text_matrix` | `_read_geo_text_matrix` | GSE131907, GSE132465 | genes × cells `.txt.gz`, streamed line-by-line into a CSR matrix (`np.fromstring`), then transposed to cells × genes. pandas is pathologically slow on the ~208k-cell-wide frame, hence the manual stream. |
| `tenx_mtx` | `_read_tenx_mtx` | Zheng68K | 10x MatrixMarket triplet; gene symbols from column 2 of `genes.tsv`. |

The batch key `sample_id` is standardised in `load_raw` (author `Sample` column, or
the dataset name for single-sample data).

## Parameters (config)

| Block | Key | Value | Meaning |
|-------|-----|-------|---------|
| `qc` | `min_genes` / `max_genes` | 200 / 6000 | gene-count window (loose upper bound for cancer) |
| `qc` | `max_pct_mt` | 20 | mito fraction cutoff (loose, pitfall #6) |
| `preprocess` | `target_sum` | 10000 | normalize_total target |
| `preprocess` | `n_top_genes` | 2000 | HVG count |
| `preprocess` | `n_pcs` | 50 | PCA dims (single-sample fallback) |
| `preprocess` | `batch_key` | `sample_id` | integration/HVG batch key |
| `integration` | `method` | `scvi` | `scvi` (multi-sample) or `pca` |
| `integration` | `n_latent` / `n_layers` | 30 / 2 | scVI architecture |
| `integration` | `batch_size` | 1024 | scVI minibatch |
| `integration` | `hvg_flavor` | `seurat_v3` | HVG on raw counts (needs scikit-misc) |

## Conventions enforced here (and where in code)

- **Raw counts preserved** — `layers["counts"] = adata.X.copy()` before
  `normalize_total` (`normalise_and_embed`). `adata.raw = adata` set after lognorm,
  before HVG subsetting. (pitfall #6)
- **Batch correction before annotation** — scVI latent → neighbours → UMAP.
  UMAP must NOT cluster by `sample_id` (verify interactively). (pitfall #3)
- **Non-malignant subset by author labels** — `~obs[lineage_col].isin(exclude)`
  → `<ds>_tme.h5ad`. (pitfall #1)
- **Finest label kept** — `Cell_subtype` retained as ground truth; broad label
  used only for the malignant filter. (pitfall #2)

## Cross-language output (assay naming)

zellkonverter reads the `.h5ad` in R with assays `X` (lognorm) and `counts` (raw).
`common.R::load_sce` aliases `X → logcounts` so SingleR works; Seurat tools use
`counts` (+ `X` as lognorm `data`).

## Run

```bash
uv run python scripts/01_preprocess.py GSE131907   # ~7 GB processed, ~5 GB tme
uv run python scripts/01_preprocess.py GSE132465
uv run python scripts/01_preprocess.py Zheng68K    # PCA path, no scVI
```

## Status

✅ All six `.h5ad` written:
`GSE131907_{processed,tme}`, `GSE132465_{processed,tme}`, `Zheng68K_{processed,tme}`.

## Paper notes / TODO

- **Doublet removal is not implemented** — `qc_filter` has a `TODO` for per-sample
  Scrublet. The guide's plan calls for it before merging; current pipeline skips it.
  State this as a limitation or add it before final results.
- Report per-dataset cell counts before/after QC and after the TME filter.
- Include the scVI-UMAP-by-sample sanity check as a supplementary figure (evidence
  that annotation follows biology, not batch).
- Methods refs: scanpy (Wolf 2018, [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)),
  scVI (Lopez 2018, [10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2)).
