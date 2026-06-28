# Run SingleR on the non-malignant TME cells with HEALTHY references.
# Two references (HPCA + Blueprint/Encode) -> two prediction files.
# Ground-truth labels are never used (project guide §7.2).
#
# Usage: Rscript scripts/r/run_singler.R GSE131907

suppressMessages({
  library(SingleR)
  library(celldex)
})
source(file.path("scripts", "r", "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript run_singler.R <dataset>")
dataset <- args[1]

cfg <- load_config()
sce <- load_tme_sce(cfg, dataset)   # expects logcounts assay from preprocessing

ref_hpca <- celldex::HumanPrimaryCellAtlasData()
ref_blueprint <- celldex::BlueprintEncodeData()

pred_hpca <- SingleR(test = sce, ref = ref_hpca,
                     labels = ref_hpca$label.fine, assay.type.test = "logcounts")
pred_bp <- SingleR(test = sce, ref = ref_blueprint,
                   labels = ref_blueprint$label.fine, assay.type.test = "logcounts")

write_predictions(cfg, dataset, "singler", colnames(sce), pred_hpca$labels)
write_predictions(cfg, dataset, "singler_blueprint", colnames(sce), pred_bp$labels)
