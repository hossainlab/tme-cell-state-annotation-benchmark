"""Preprocess a raw dataset into an annotated, batch-corrected .h5ad.

Output: data/raw/<dataset>/<dataset>_processed.h5ad
        + data/raw/<dataset>/<dataset>_tme.h5ad   (non-malignant cells only)

The .h5ad is the cross-language handoff: Python tools read it directly,
R tools read it via zellkonverter::readH5AD().

Pipeline (project guide §9, step 2). Critical conventions enforced here:
  - raw UMI counts kept in layers["counts"]; adata.raw saved after lognorm
  - loose cancer QC thresholds (mt < 20%, genes < 6000)
  - scVI batch integration BEFORE annotation (PCA fallback for single-sample data)
  - filter to non-malignant lineages using AUTHOR labels only

Two input formats are handled (set per dataset in config.yaml -> format):
  - geo_text_matrix : a genes x cells .txt.gz UMI matrix (GSE131907, GSE132465).
                      Too large to hold dense, so read in gene chunks into a
                      sparse matrix, then transpose to cells x genes.
  - tenx_mtx        : a 10x MatrixMarket triplet + a separate label tsv (Zheng68K).

Usage:
    python scripts/01_preprocess.py GSE131907
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from common import load_config, repo_path


def _read_geo_text_matrix(path, chunksize: int = 1000) -> sc.AnnData:
    """Read a genes x cells .txt(.gz) UMI matrix as a cells x genes AnnData.

    These matrices are very WIDE (~208k cell columns). pandas is pathologically
    slow on frames that wide (per-column Python overhead), so we stream the file
    line by line instead: each gene row is parsed with np.fromstring (C speed)
    and its non-zeros are appended directly into a CSR matrix. Stays sparse the
    whole way, so memory is ~the nnz, not genes*cells dense.

    `chunksize` is reused here only as the progress-print interval (in genes).
    """
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    indptr, indices, data, genes = [0], [], [], []
    with opener(path, "rt") as fh:
        cell_names = fh.readline().rstrip("\n").split("\t")[1:]   # drop "Index"
        ncols = len(cell_names)
        for i, line in enumerate(fh):
            tab = line.index("\t")
            genes.append(line[:tab])
            vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
            if vals.shape[0] != ncols:
                raise ValueError(f"row {i} ({line[:tab]}) has {vals.shape[0]} "
                                 f"values, expected {ncols}")
            nz = np.nonzero(vals)[0]
            indices.append(nz.astype(np.int32))
            data.append(vals[nz])
            indptr.append(indptr[-1] + nz.shape[0])
            if (i + 1) % chunksize == 0:
                print(f"  read {i + 1} genes...", flush=True)

    mat = sp.csr_matrix(                                         # genes x cells
        (np.concatenate(data), np.concatenate(indices), np.asarray(indptr)),
        shape=(len(genes), ncols),
    )
    adata = sc.AnnData(
        X=mat.T.tocsr(),                                        # -> cells x genes
        obs=pd.DataFrame(index=pd.Index(cell_names, name="cell_id")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.var_names_make_unique()
    return adata


def _read_tenx_mtx(spec: dict) -> sc.AnnData:
    """Read a 10x MatrixMarket triplet (matrix.mtx + barcodes.tsv + genes.tsv)."""
    mtx_dir = repo_path(spec["mtx_dir"])
    adata = sc.read_mtx(mtx_dir / "matrix.mtx").T            # mtx is genes x cells
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", header=None)[0].to_numpy()
    genes = pd.read_csv(mtx_dir / "genes.tsv", header=None, sep="\t")
    adata.obs_names = barcodes
    adata.var_names = genes[1].to_numpy() if genes.shape[1] > 1 else genes[0].to_numpy()
    adata.var_names_make_unique()
    return adata


def load_raw(cfg: dict, dataset: str) -> sc.AnnData:
    """Load the count matrix and join author annotations on the barcode column."""
    spec = cfg["datasets"][dataset]
    fmt = spec["format"]
    if fmt == "geo_text_matrix":
        chunksize = cfg.get("compute", {}).get("read_chunksize", 1000)
        adata = _read_geo_text_matrix(repo_path(spec["matrix"]), chunksize=chunksize)
    elif fmt == "tenx_mtx":
        adata = _read_tenx_mtx(spec)
    else:
        sys.exit(f"unknown format {fmt!r} for {dataset}")

    ann = pd.read_csv(repo_path(spec["annotation"]), sep=spec["annotation_sep"])
    ann = ann.set_index(spec["barcode_column"])
    keep = [c for c in (spec.get("truth_column"), spec.get("lineage_column"),
                        spec.get("condition_column"), spec.get("sample_column"))
            if c and c in ann.columns]
    adata.obs = adata.obs.join(ann[keep].reindex(adata.obs_names))

    # Standardise the batch key the rest of the pipeline expects.
    sample_col = spec.get("sample_column")
    adata.obs["sample_id"] = (adata.obs[sample_col].astype(str)
                              if sample_col else dataset)

    n_missing = int(adata.obs[spec["truth_column"]].isna().sum())
    if n_missing:
        print(f"warning: {n_missing}/{adata.n_obs} cells lack a ground-truth label "
              "(barcode mismatch between matrix and annotation?)")
    return adata


def qc_filter(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    q = cfg["qc"]
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts > q["min_genes"]]
    adata = adata[adata.obs.n_genes_by_counts < q["max_genes"]]
    adata = adata[adata.obs.pct_counts_mt < q["max_pct_mt"]]
    # TODO: per-sample Scrublet doublet removal could be added here.
    return adata.copy()


def normalise_and_embed(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    p = cfg["preprocess"]
    ig = cfg["integration"]
    batch_key = p["batch_key"]
    n_batches = adata.obs[batch_key].nunique()

    # scVI needs raw UMI counts — stash them before we normalise X in place.
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=p["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata  # keep lognorm for tools/figures BEFORE HVG subsetting

    # HVG on raw counts (seurat_v3); batch-aware only when there is >1 sample.
    sc.pp.highly_variable_genes(
        adata, n_top_genes=p["n_top_genes"], flavor=ig["hvg_flavor"],
        layer="counts", batch_key=batch_key if n_batches > 1 else None,
    )

    # Batch integration BEFORE annotation (pitfall #3).
    if n_batches > 1 and ig["method"] == "scvi":
        import torch
        from scvi.model import SCVI

        accel = ("gpu" if cfg["compute"].get("gpu") and torch.cuda.is_available()
                 else "cpu")
        hvg = adata[:, adata.var.highly_variable].copy()   # scVI trains on HVGs
        SCVI.setup_anndata(hvg, layer="counts", batch_key=batch_key)
        model = SCVI(hvg, n_latent=ig["n_latent"], n_layers=ig["n_layers"])
        model.train(max_epochs=ig["max_epochs"], accelerator=accel, devices=1,
                    batch_size=ig["batch_size"])
        adata.obsm["X_scVI"] = model.get_latent_representation()
        print(f"scVI integration done on {accel} "
              f"({n_batches} batches -> X_scVI {adata.obsm['X_scVI'].shape})")
        rep = "X_scVI"
    else:
        print(f"single batch ({n_batches}) — PCA, no integration")
        sc.pp.pca(adata, n_comps=p["n_pcs"])
        rep = "X_pca"

    sc.pp.neighbors(adata, use_rep=rep)
    sc.tl.umap(adata)
    # Sanity check to run interactively: UMAP must NOT cluster by sample_id.
    return adata


def main(dataset: str) -> None:
    cfg = load_config()
    sc.settings.n_jobs = cfg.get("compute", {}).get("n_cores", 1)  # PCA/neighbors/UMAP threads
    spec = cfg["datasets"][dataset]
    out_dir = repo_path(cfg["paths"]["data_raw"], dataset)

    adata = load_raw(cfg, dataset)
    adata = qc_filter(adata, cfg)
    adata = normalise_and_embed(adata, cfg)

    processed = out_dir / f"{dataset}_processed.h5ad"
    adata.write_h5ad(processed)
    print("wrote", processed, f"({adata.n_obs} cells)")

    # Non-malignant subset for the general tools (pitfall #1): drop the malignant
    # / junk lineages named per dataset; keep everything else.
    lineage_col = spec.get("lineage_column")
    exclude = spec.get("exclude_lineages", []) or []
    if lineage_col and lineage_col in adata.obs:
        mask = ~adata.obs[lineage_col].astype(str).isin(exclude)
    else:
        mask = pd.Series(True, index=adata.obs_names)  # healthy: keep all
    tme = adata[mask.values].copy()
    tme_path = out_dir / f"{dataset}_tme.h5ad"
    tme.write_h5ad(tme_path)
    print("wrote", tme_path, f"({tme.n_obs} non-malignant cells)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
