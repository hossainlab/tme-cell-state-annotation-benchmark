# Data Download Guide

How to obtain the three datasets and lay them out so the pipeline finds them.
All paths below match `config/config.yaml`. Nothing here is run automatically —
follow it manually.

> **Disk budget:** GSE131907 ≈ 5–8 GB compressed (much more decompressed),
> GSE132465 ≈ 1–2 GB, Zheng68K ≈ 0.5 GB. Have ~40 GB free before starting.
> These raw files are gitignored — never commit them.

---

## Target layout

```
data/raw/
├── GSE131907/
│   ├── GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
│   └── GSE131907_Lung_Cancer_cell_annotation.txt.gz
├── GSE132465/
│   ├── GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz
│   └── GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz
└── Zheng68K/
    └── filtered_matrices_mex/            (10x mtx triplet) — or use scanpy shortcut
```

---

## 1. GSE131907 — lung adenocarcinoma (primary ground truth)

GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907>
Paper: Kim et al. 2020, *Nat Commun* — <https://doi.org/10.1038/s41467-020-16164-1>

Download the supplementary files:

```bash
mkdir -p data/raw/GSE131907
cd data/raw/GSE131907

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl"
wget -c "$BASE/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
wget -c "$BASE/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
# optional (already-normalised matrix, if you prefer not to renormalise):
# wget -c "$BASE/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
cd -
```

`-c` resumes interrupted downloads — important for the multi-GB matrix.

To browse the exact filenames first (they occasionally differ):

```bash
wget -qO- "$BASE/" | grep -oE 'GSE131907[^"]+\.(txt|gz|RData)'
```

> ⚠️ **Format note.** The scaffolded `scripts/01_preprocess.py::load_raw` currently
> assumes *per-sample `.h5`* files. GSE131907 actually ships **one big genes×cells
> text matrix**, not per-sample h5. You will need to either (a) adapt `load_raw` to
> read the `.txt.gz` matrix with `pd.read_csv(..., sep="\t")` / `sc.read_text`, or
> (b) split it per sample. Flag this when you start preprocessing.

---

## 2. GSE132465 — colorectal (replication)

GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465>
Paper: Lee et al. 2020, *Nat Genet* — <https://doi.org/10.1038/s41588-020-0636-z>

```bash
mkdir -p data/raw/GSE132465
cd data/raw/GSE132465

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl"
wget -c "$BASE/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
wget -c "$BASE/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
cd -
```

---

## 3. Zheng68K — healthy PBMC (accuracy ceiling)

Two options.

**Option A — quick, in-code (recommended to start).** A subsampled version ships
with scanpy; no manual download:

```python
import scanpy as sc
adata = sc.datasets.pbmc68k_reduced()   # ~700 cells, labelled — fine for a smoke test
```

**Option B — full 68k dataset from 10x Genomics.** For the real baseline:

```bash
mkdir -p data/raw/Zheng68K
cd data/raw/Zheng68K
wget -c "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
tar -xzf fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
cd -
```

The per-cell ground-truth labels for Zheng68K (the 11 immune subtypes) come from
the original SeqWell/Zheng 2017 annotation, distributed with most benchmark
repos (e.g. the CellTypist / scGPT example data). Note where you get them so the
provenance is citable.

---

## Verify downloads

Confirm the gzip files are intact (no truncation from a dropped connection):

```bash
for f in $(find data/raw -name '*.gz'); do
  gzip -t "$f" && echo "OK  $f" || echo "BAD $f  (re-download with wget -c)"
done
```

Peek at a matrix header without decompressing the whole thing:

```bash
zcat data/raw/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz | head -3
```

Check the annotation column names match `config.yaml` (`truth_column`,
`condition_column`, `sample_column`) — GEO column names sometimes differ from
the paper; update the config to whatever the file actually uses.

---

## Next step

Once files are in place and verified:

```bash
uv run python scripts/01_preprocess.py GSE131907
```

(after adapting `load_raw` to the real matrix format, per the note above).
