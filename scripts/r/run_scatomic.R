# Run scATOMIC, the cancer-specialised comparison tool.
# UNLIKE the general tools, scATOMIC runs on the FULL tumour matrix (malignant +
# non-malignant) because identifying malignant cells is part of its job. We then
# keep only its non-malignant predictions for the benchmark comparison.
#
# Usage: Rscript scripts/r/run_scatomic.R GSE131907

suppressMessages({
  library(scATOMIC)
  library(zellkonverter)
  library(SingleCellExperiment)
})
source(file.path("scripts", "r", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript run_scatomic.R <dataset>")
dataset <- args[1]

cfg <- load_config()
# Full processed matrix (NOT the TME-only subset) — scATOMIC needs cancer cells.
processed <- file.path(repo_root(), cfg$paths$data_raw, dataset,
                       paste0(dataset, "_processed.h5ad"))
sce <- readH5AD(processed)
counts <- as.matrix(assay(sce, "X"))

results <- run_scATOMIC(counts, mc.cores = 4)
pred <- create_summary_matrix(prediction_list = results, raw_counts = counts)

write_predictions(cfg, dataset, "scatomic",
                  rownames(pred), pred$scATOMIC_pred)
