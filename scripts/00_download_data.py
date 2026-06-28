"""Download raw datasets into data/raw/<dataset>/.

Usage:
    python scripts/00_download_data.py GSE131907
    python scripts/00_download_data.py GSE132465

Zheng68K is fetched in-code via scanpy (sc.datasets.pbmc68k_reduced) or
downloaded manually from 10x Genomics; no GEO mirror needed.

This is a thin wrapper around wget so downloads are reproducible and logged.
Supplementary files on GEO are large — expect tens of GB for GSE131907.
"""
from __future__ import annotations

import subprocess
import sys

from common import load_config, repo_path


def download(dataset: str) -> None:
    cfg = load_config()
    spec = cfg["datasets"].get(dataset)
    if spec is None:
        sys.exit(f"unknown dataset {dataset!r}; see config/config.yaml")
    url = spec.get("download_url")
    if not url:
        sys.exit(f"{dataset} has no download_url (fetch it in-code, e.g. Zheng68K)")

    out_dir = repo_path(cfg["paths"]["data_raw"], dataset)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["wget", "-r", "-np", "-nH", "--cut-dirs=100", "-P", str(out_dir), url]
    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    download(sys.argv[1])
