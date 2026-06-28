# Run Azimuth reference mapping (healthy lung reference) on the TME cells.
# The reference is healthy, so degradation on cancer is the expected, measurable
# effect. Ground-truth labels are never used (project guide §7.2).
#
# Usage: Rscript scripts/r/run_azimuth.R GSE131907

suppressMessages({
  library(Seurat)
  library(Azimuth)
  library(SingleCellExperiment)
})
source(file.path("scripts", "r", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript run_azimuth.R <dataset>")
dataset <- args[1]

cfg <- load_config()
sce <- load_tme_sce(cfg, dataset)
seurat <- as.Seurat(sce, counts = "X", data = NULL)

# TODO: install/point to the Azimuth healthy lung reference ("lungref").
# RunAzimuth annotates predicted.* columns at several reference levels.
seurat <- RunAzimuth(seurat, reference = "lungref")

labels <- seurat$predicted.ann_level_3  # finest Azimuth lung level
write_predictions(cfg, dataset, "azimuth", colnames(seurat), as.character(labels))
