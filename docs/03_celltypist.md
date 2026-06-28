# 03 — CellTypist

**Script:** `scripts/02_run_celltypist.py <dataset>`
**Input:** `<ds>_tme.h5ad` (non-malignant cells)
**Output:** `results/predictions/<ds>__celltypist.csv`
**Language:** Python (main `.venv`)

## Methods (paper-ready)

CellTypist (Domínguez Conde et al. 2022, *Science*) was run on the non-malignant
TME cells using its pretrained `Immune_All_Low` model — the flagship general
logistic-regression model covering 98 fine-grained immune subtypes, including the
macrophage and T-cell states this benchmark targets. Annotation used the model's
majority-voting mode, which over-clusters the data and assigns a consensus label
per local neighbourhood. The model was applied zero-shot; ground-truth labels were
never provided. Predictions were written per cell as `cell_id, predicted_label`.

## Parameters (config `celltypist`)

| Key | Value | Meaning |
|-----|-------|---------|
| `model` | `Immune_All_Low.pkl` | 98 fine immune subtypes (immune-only) |
| `majority_voting` | `true` | consensus per over-clustered neighbourhood |

The label column read is `majority_voting` when voting is on, else
`predicted_labels`.

## Known limitation (expected finding, not a bug)

`Immune_All_Low` is **immune-only**: fibroblast/CAF and epithelial cells get
mislabelled to the nearest immune type. For a general immune annotator this is an
expected limitation and is itself part of the benchmark result — report it, do not
"fix" it. (`Pan_Human_Atlas` from earlier drafts is **not** a real CellTypist model.)

## Run

```bash
uv run python scripts/02_run_celltypist.py GSE131907
uv run python scripts/02_run_celltypist.py GSE132465
uv run python scripts/02_run_celltypist.py Zheng68K
```

## Status

✅ Done for all three datasets:
`GSE131907__celltypist.csv`, `GSE132465__celltypist.csv`, `Zheng68K__celltypist.csv`.

## Paper notes

- Cite CellTypist: Domínguez Conde et al. 2022, *Science*
  [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197). Record the
  model version/date (`models.download_models` fetches the current release).
- This tool is the healthy-immune specialist: expect high L1/L2 immune accuracy,
  collapse at L3 functional states, and systematic CAF→immune errors feeding the
  confusion-matrix figure.
