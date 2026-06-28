# Shared R helpers: config loading, .h5ad import, prediction output.
# Mirrors scripts/python/common.py so both languages stay consistent.

suppressMessages({
  library(yaml)
})

repo_root <- function() {
  # this file lives at scripts/r/common.R
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "scripts/r/common.R"),
                          "..", ".."))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

load_config <- function() {
  yaml::read_yaml(file.path(repo_root(), "config", "config.yaml"))
}

# Load the non-malignant TME .h5ad written by 01_preprocess.py as a
# SingleCellExperiment (zellkonverter). Seurat-based tools convert from there.
load_tme_sce <- function(cfg, dataset) {
  suppressMessages(library(zellkonverter))
  path <- file.path(repo_root(), cfg$paths$data_raw, dataset,
                    paste0(dataset, "_tme.h5ad"))
  if (!file.exists(path)) stop("missing ", path, " — run 01_preprocess.py first")
  zellkonverter::readH5AD(path)
}

# Write predictions in the canonical schema: cell_id, predicted_label.
write_predictions <- function(cfg, dataset, tool, cell_ids, labels) {
  out_dir <- file.path(repo_root(), cfg$paths$predictions)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, sprintf("%s__%s.csv", dataset, tool))
  write.csv(data.frame(cell_id = cell_ids, predicted_label = labels),
            out, row.names = FALSE)
  message("wrote ", out)
}
