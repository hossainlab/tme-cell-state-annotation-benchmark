"""Preprocess a raw dataset into an annotated, batch-corrected .h5ad.

Output: data/raw/<dataset>/<dataset>_processed.h5ad
        + data/raw/<dataset>/<dataset>_tme.h5ad   (non-malignant cells only)

The .h5ad is the cross-language handoff: Python tools read it directly,
R tools read it via zellkonverter::readH5AD().

Pipeline (project guide §9, step 2). Critical conventions enforced here:
  - adata.raw saved BEFORE normalisation (some tools need raw counts)
  - loose cancer QC thresholds (mt < 20%, genes < 6000)
  - Harmony batch correction BEFORE annotation
  - filter to non-malignant lineages using AUTHOR labels only

Usage:
    python scripts/01_preprocess.py GSE131907
"""
from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc

from common import load_config, repo_path


def load_raw(cfg: dict, dataset: str) -> sc.AnnData:
    """Load per-sample matrices and join author annotations.

    TODO: GSE131907 ships per-sample .h5 files; GSE132465 ships one big matrix.
    Branch on dataset spec. This stub shows the GSE131907 shape from the guide.
    """
    spec = cfg["datasets"][dataset]
    raw_dir = repo_path(cfg["paths"]["data_raw"], dataset)

    h5_files = sorted(raw_dir.glob("*.h5"))
    if not h5_files:
        sys.exit(f"no .h5 under {raw_dir} — run 00_download_data.py first")

    adatas = []
    for h5 in h5_files:
        a = sc.read_10x_h5(h5)
        a.var_names_make_unique()
        a.obs["sample_id"] = h5.stem
        adatas.append(a)
    adata = sc.concat(adatas)

    ann = pd.read_csv(spec["annotation"], sep=spec["annotation_sep"], index_col=0)
    # Keep the FINEST author label as ground truth (pitfall #2). Plus condition.
    keep = [c for c in (spec["truth_column"], spec.get("condition_column"),
                        spec.get("sample_column")) if c and c in ann.columns]
    adata.obs = adata.obs.join(ann[keep])
    return adata


def qc_filter(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    q = cfg["qc"]
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts > q["min_genes"]]
    adata = adata[adata.obs.n_genes_by_counts < q["max_genes"]]
    adata = adata[adata.obs.pct_counts_mt < q["max_pct_mt"]]
    # TODO: per-sample Scrublet doublet removal here, BEFORE the concat ideally.
    return adata.copy()


def normalise_and_embed(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    p = cfg["preprocess"]
    sc.pp.normalize_total(adata, target_sum=p["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata  # save normalised+raw BEFORE HVG subsetting / scaling

    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"])
    sc.pp.pca(adata, n_comps=p["n_pcs"])

    # Harmony batch correction BEFORE annotation (pitfall #3).
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, p["batch_key"])
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    # Sanity check to run interactively: UMAP must NOT cluster by sample_id.
    return adata


def main(dataset: str) -> None:
    cfg = load_config()
    out_dir = repo_path(cfg["paths"]["data_raw"], dataset)

    adata = load_raw(cfg, dataset)
    adata = qc_filter(adata, cfg)
    adata = normalise_and_embed(adata, cfg)

    processed = out_dir / f"{dataset}_processed.h5ad"
    adata.write_h5ad(processed)
    print("wrote", processed)

    # Non-malignant subset for the general tools (pitfall #1).
    truth_col = cfg["datasets"][dataset]["truth_column"]
    lineages = cfg["non_malignant_lineages"]
    # TODO: replace with the dataset's actual broad-lineage column once confirmed.
    mask = adata.obs[truth_col].astype(str).str.contains("|".join(lineages), case=False)
    tme = adata[mask].copy()
    tme_path = out_dir / f"{dataset}_tme.h5ad"
    tme.write_h5ad(tme_path)
    print("wrote", tme_path, f"({tme.n_obs} non-malignant cells)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
