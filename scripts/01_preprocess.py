"""Preprocess a raw dataset into an annotated, batch-corrected .h5ad.

Output: data/raw/<dataset>/<dataset>_processed.h5ad
        + data/raw/<dataset>/<dataset>_tme.h5ad   (non-malignant cells only)

The .h5ad is the cross-language handoff: Python tools read it directly,
R tools read it via zellkonverter::readH5AD().

Pipeline (project guide §9, step 2). Critical conventions enforced here:
  - adata.raw saved BEFORE normalisation (some tools need raw counts)
  - loose cancer QC thresholds (mt < 20%, genes < 6000)
  - Harmony batch correction BEFORE annotation (skipped for single-sample data)
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

    The full matrix (e.g. ~30k genes x ~208k cells) is far too large to hold
    dense, so stream it in row (gene) chunks, sparsify each block, and vstack.
    """
    # NB: do NOT pass dtype=float to read_csv — it would also try to cast the
    # gene-name index column ("A1BG", ...). Cast the numeric block only.
    reader = pd.read_csv(path, sep="\t", index_col=0, chunksize=chunksize)
    blocks, gene_names, cell_names = [], [], None
    for chunk in reader:
        if cell_names is None:
            cell_names = chunk.columns.to_numpy()
        gene_names.append(chunk.index.to_numpy())
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))   # genes(chunk) x cells
    mat = sp.vstack(blocks).T.tocsr()                       # -> cells x genes
    adata = sc.AnnData(
        X=mat,
        obs=pd.DataFrame(index=pd.Index(cell_names, name="cell_id")),
        var=pd.DataFrame(index=pd.Index(np.concatenate(gene_names), name="gene")),
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
    sc.pp.normalize_total(adata, target_sum=p["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata  # save normalised counts BEFORE HVG subsetting / scaling

    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"])
    sc.pp.pca(adata, n_comps=p["n_pcs"])

    # Harmony batch correction BEFORE annotation (pitfall #3) — only if >1 sample.
    n_batches = adata.obs[p["batch_key"]].nunique()
    if n_batches > 1:
        import harmonypy as hm
        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, p["batch_key"])
        adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        rep = "X_pca_harmony"
    else:
        print(f"single batch ({n_batches}) — skipping Harmony")
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
