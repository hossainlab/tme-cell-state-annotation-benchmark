# R dependencies for the R-language tools (SingleR, Azimuth, GPTCelltype, scATOMIC).
# Run once: Rscript envs/install_r.R
#
# Installs into the per-user library; create it first if it does not exist:
#   mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4

# Rscript is non-interactive and has no default CRAN mirror — set one explicitly,
# otherwise install.packages() errors with "trying to use CRAN without setting a mirror".
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages.
BiocManager::install(c(
  "SingleR",
  "celldex",
  "scuttle",
  "SingleCellExperiment",
  "zellkonverter"    # read AnnData .h5ad written by the Python preprocessing step
), update = FALSE, ask = FALSE)

# CRAN packages.
install.packages(c("Seurat", "remotes", "dplyr"))

# GitHub-only packages. Azimuth is NOT on Bioconductor/CRAN — it ships from
# satijalab/azimuth alongside SeuratData.
remotes::install_github("satijalab/seurat-data")
remotes::install_github("satijalab/azimuth")
remotes::install_github("Winnie09/GPTCelltype")
remotes::install_github("abelson-lab/scATOMIC")
