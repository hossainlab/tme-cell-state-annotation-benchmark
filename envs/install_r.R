# R dependencies for the R-language tools (SingleR, Azimuth, GPTCelltype, scATOMIC).
# Run once: Rscript envs/install_r.R

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "SingleR",
  "celldex",
  "scuttle",
  "SingleCellExperiment",
  "zellkonverter",   # read AnnData .h5ad written by the Python preprocessing step
  "Azimuth"
))

install.packages(c("Seurat", "remotes", "dplyr"))

# GPTCelltype + scATOMIC are GitHub-only.
remotes::install_github("Winnie09/GPTCelltype")
remotes::install_github("abelson-lab/scATOMIC")
