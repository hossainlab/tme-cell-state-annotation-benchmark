# 07 — GPTCelltype

**Script:** `scripts/run_gptcelltype.R <dataset>`
**Input:** `<ds>_tme.h5ad` → Seurat (counts + lognorm), clustered
**Output:** `results/predictions/<ds>__gptcelltype_lungtumour.csv` and
`<ds>__gptcelltype_lung.csv`
**Language:** R + OpenAI API (GPT-4)

## Methods (paper-ready)

GPTCelltype (Hou & Ji 2024, *Nat Methods*) annotates clusters from their marker
genes via a large language model (GPT-4). The non-malignant TME cells were
clustered (Seurat `FindNeighbors` on 30 dims, `FindClusters` resolution 1.0) and
positive markers identified per cluster (`FindAllMarkers`, min.pct 0.25, logFC 0.25).
Because the model is non-deterministic, each cluster was annotated **3 independent
times**; the modal label was taken as the cluster call and the inter-run
disagreement rate (`1 − max_vote_fraction`) reported as a stability metric. To
quantify the cancer-context effect, the whole procedure was run twice with
`tissuename = "lung tumour"` and `tissuename = "lung"`, producing two prediction
files. Cluster labels were broadcast to constituent cells. Ground-truth labels were
never used.

## Parameters (config `gptcelltype`)

| Key | Value | Meaning |
|-----|-------|---------|
| `model` | `gpt-4` | OpenAI model |
| `n_runs` | 3 | repeats per cluster (non-determinism, pitfall #4) |
| `tissuenames` | `["lung tumour", "lung"]` | context-effect contrast |

Clustering: `FindNeighbors(dims=1:30)`, `FindClusters(resolution=1.0)`.

## Two reported metrics (beyond accuracy)

- **Modal label** per cluster → the prediction.
- **Mean cluster disagreement rate** across the 3 runs → reproducibility metric
  (printed per tissue; capture it for the paper).

## Blockers

1. **`OPENAI_API_KEY` is unset** — the script stops immediately without it.
2. Cost/rate: 3 runs × 2 tissues × clusters of GPT-4 calls per dataset.

## Dataset note

`tissuenames` is lung-specific. For GSE132465 use `["colorectal tumour", "colon"]`;
for Zheng68K use `["blood"]` (single context — no tumour contrast). Make this
config-driven per dataset; the contrast only makes sense for the tumour datasets.

## Run

```bash
export OPENAI_API_KEY=sk-...
Rscript scripts/run_gptcelltype.R GSE131907
```

## Status

⬜ Not run — no API key; `GPTCelltype` (1.0.1) package is installed.

## Paper notes

- Cite GPTCelltype: Hou & Ji 2024, *Nat Methods*
  [10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4).
- Record GPT-4 model snapshot/date — results are not reproducible across model
  versions. The "lung tumour" vs "lung" delta is a headline result: prompt context
  changing the answer.
- Disagreement rate is itself a finding (LLM annotation instability).
