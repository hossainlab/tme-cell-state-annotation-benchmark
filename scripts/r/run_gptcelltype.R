# Run GPTCelltype (GPT-4) on cluster marker genes.
# Two deliberate design points (project guide §12, pitfall #4):
#   - non-deterministic: run n_runs times, report modal label + disagreement
#   - run with BOTH tissuename = "lung tumour" and "lung" to measure the
#     cancer-context effect
# Needs OPENAI_API_KEY in the environment.
#
# Usage: Rscript scripts/r/run_gptcelltype.R GSE131907

suppressMessages({
  library(Seurat)
  library(GPTCelltype)
  library(dplyr)
})
source(file.path("scripts", "r", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript run_gptcelltype.R <dataset>")
dataset <- args[1]
if (Sys.getenv("OPENAI_API_KEY") == "") stop("set OPENAI_API_KEY")

cfg <- load_config()
sce <- load_tme_sce(cfg, dataset)
seurat <- as.Seurat(sce, counts = "X", data = "X")

seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 1.0)
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)

n_runs <- cfg$gptcelltype$n_runs
for (tissue in cfg$gptcelltype$tissuenames) {
  runs <- replicate(n_runs, gptcelltype(
    input = markers, tissuename = tissue, model = cfg$gptcelltype$model
  ), simplify = FALSE)

  # Modal label per cluster + inter-run disagreement rate.
  mat <- do.call(cbind, runs)
  modal <- apply(mat, 1, function(r) names(sort(table(r), decreasing = TRUE))[1])
  disagree <- apply(mat, 1, function(r) 1 - max(table(r)) / length(r))

  cluster_label <- setNames(modal, rownames(mat))
  cell_labels <- cluster_label[as.character(Idents(seurat))]
  tag <- gsub("[^a-z]", "", tolower(tissue))  # lungtumour / lung
  write_predictions(cfg, dataset, paste0("gptcelltype_", tag),
                    colnames(seurat), cell_labels)

  message(sprintf("%s mean cluster disagreement: %.3f", tissue, mean(disagree)))
}
