# R dependencies for the R-language tools (SingleR, Azimuth, GPTCelltype, scATOMIC).
# Run once: Rscript envs/install_r.R
#
# Installs into the per-user library; create it first if it does not exist:
#   mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.4
#
# SYSTEM DEPENDENCIES (apt — need sudo, install BEFORE running this):
#   sudo apt-get install -y libgsl-dev libhdf5-dev libfftw3-dev
# libgsl-dev is required to compile DirichletMultinomial -> TFBSTools -> Signac ->
# Azimuth. Without it the Azimuth chain fails with "gsl/gsl_rng.h: No such file".

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

# scATOMIC depends on Rmagic, which was removed from CRAN — install it (and its
# helper cutoff.scATOMIC) from GitHub before scATOMIC itself.
# Rmagic also needs the python backend at RUNTIME:  uv pip install magic-impute
remotes::install_github("KrishnaswamyLab/MAGIC", subdir = "Rmagic")
remotes::install_github("abelson-lab/scATOMIC")
