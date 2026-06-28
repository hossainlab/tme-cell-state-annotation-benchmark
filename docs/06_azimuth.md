# 06 — Azimuth

**Script:** `scripts/run_azimuth.R <dataset>`
**Input:** `<ds>_tme.h5ad` → Seurat object (raw counts)
**Output:** `results/predictions/<ds>__azimuth.csv`
**Language:** R / Seurat

## Methods (paper-ready)

Azimuth (Hao et al. 2021, *Cell*) was used for supervised reference mapping against
a healthy reference atlas. The non-malignant TME cells were converted to a Seurat
object from raw counts, SCTransform-normalised within Azimuth, and mapped onto the
reference; the finest reference annotation level was taken as the predicted label.
The reference is healthy tissue, so degradation on tumour-specific states is the
expected, measurable effect. Ground-truth labels were never used.

## Parameters

| Item | Value | Notes |
|------|-------|-------|
| reference | `lungref` (healthy lung) | **must be installed/pointed to** |
| label level | `predicted.ann_level_3` | finest Azimuth lung level |
| counts | `counts` assay | `as.Seurat(sce, counts="counts", data=NULL)` |
| parallel | `future` multicore | `setup_future` (n_cores) |

## Blockers

1. **R package `Azimuth` is not installed** (verified missing). Reinstall via
   `envs/install_r.R`; needs `libgsl-dev` (apt) to compile the
   DirichletMultinomial → TFBSTools → Signac dependency chain.
2. **Reference `lungref` must be available.** The script has a `TODO` to
   install/point to the Azimuth healthy lung reference; `RunAzimuth(reference="lungref")`
   will fail until it is downloaded.

## Dataset applicability

| Dataset | Reference to use |
|---------|------------------|
| GSE131907 (lung) | `lungref` (as coded) |
| GSE132465 (colorectal) | needs a colon/CRC reference — `lungref` is wrong tissue |
| Zheng68K (PBMC) | needs `pbmcref` |

> The script hard-codes `lungref`. For the replication and healthy-baseline datasets
> the reference must be swapped (ideally config-driven). Document the per-dataset
> reference choice; using a lung reference on PBMC/colon would confound results.

## Run

```bash
Rscript scripts/run_azimuth.R GSE131907   # after installing Azimuth + lungref
```

## Status

⬜ Not run — `Azimuth` package missing; reference not installed; non-lung datasets
need a different reference.

## Paper notes

- Cite Azimuth/Seurat v4: Hao et al. 2021, *Cell*
  [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048).
- Record each Azimuth reference name + version. State that Azimuth's annotation
  granularity is fixed by the reference (`ann_level_*`), unlike the free-form tools.
