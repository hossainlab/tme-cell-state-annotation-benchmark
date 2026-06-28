# 01 — Data Acquisition

**Script:** `scripts/00_download_data.py` (GEO datasets) · in-code / 10x for Zheng68K
**Practical how-to:** [`data_download.md`](data_download.md)
**Output:** `data/raw/<dataset>/…` (gitignored)

## Methods (paper-ready)

Three published single-cell RNA-seq datasets were used. Two are tumour datasets
providing the test condition and one is a healthy-tissue dataset providing the
accuracy ceiling:

- **GSE131907** (primary) — metastatic lung adenocarcinoma, 58 samples across
  primary tumour, normal lung, lymph node, brain metastasis and pleural effusion
  (Kim et al. 2020, *Nat Commun*). Distributed as a single genes × cells raw-UMI
  text matrix plus a per-cell author annotation with broad (`Cell_type`) and fine
  (`Cell_subtype`) labels.
- **GSE132465** (replication) — colorectal cancer 10x dataset (Lee et al. 2020,
  *Nat Genet*); same file structure as GSE131907.
- **Zheng68K** (healthy baseline) — ~68,000 human PBMCs sorted into 11 immune
  populations (Zheng et al. 2017, *Nat Commun*), 10x MatrixMarket triplet with a
  separate barcode→label annotation.

GEO supplementary files were retrieved over HTTPS from the NCBI FTP mirror; the
Zheng68K matrix was obtained from the 10x Genomics public sample archive. Author
cell-type annotations were retained verbatim and used only as ground truth at
evaluation time — never exposed to any annotation tool.

## What the script does

`00_download_data.py <dataset>` reads `datasets.<ds>.download_url` from the config
and mirrors the GEO suppl directory with `wget -r -np -nH --cut-dirs=100` into
`data/raw/<dataset>/`. Datasets with no `download_url` (Zheng68K) are fetched
manually or in-code — see `data_download.md`.

```bash
uv run python scripts/00_download_data.py GSE131907
uv run python scripts/00_download_data.py GSE132465
# Zheng68K: manual 10x download or sc.datasets.pbmc68k_reduced() smoke test
```

## File inventory (per config)

| Dataset | Matrix | Annotation | Format |
|---------|--------|-----------|--------|
| GSE131907 | `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` | `GSE131907_Lung_Cancer_cell_annotation.txt.gz` | `geo_text_matrix` |
| GSE132465 | `GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz` | `GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz` | `geo_text_matrix` |
| Zheng68K | `filtered_matrices_mex/hg19/` (mtx triplet) | `68k_pbmc_barcodes_annotation.tsv` | `tenx_mtx` |

## Key annotation columns (config `datasets.<ds>`)

| Role | GSE131907 | GSE132465 | Zheng68K |
|------|-----------|-----------|----------|
| barcode | `Index` | `Index` | `barcodes` |
| ground truth (finest) | `Cell_subtype` | `Cell_subtype` | `celltype` |
| broad lineage (filter) | `Cell_type` | `Cell_type` | — |
| condition | `Sample_Origin` | `Class` | — |
| batch / sample | `Sample` | `Sample` | none (single sample) |
| excluded lineages | `Epithelial cells`, `Undetermined` | `Epithelial cells` | none |

## Verification

```bash
for f in $(find data/raw -name '*.gz'); do gzip -t "$f" && echo "OK $f"; done
zcat data/raw/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz | head -3
```
Confirm header column names match the config; GEO names occasionally differ from
the paper's prose.

## Status

✅ All three datasets present and verified under `data/raw/`.
Disk budget ≈ 40 GB free recommended (GSE131907 alone is several GB compressed,
much larger decompressed). Raw files are gitignored — never commit.

## Paper notes

- Cite Kim 2020 ([10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)),
  Lee 2020 ([10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z)),
  Zheng 2017 ([10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049)).
- Record exact accession URLs and the Zheng68K label provenance for the Data
  Availability statement.
