# TME Cell State Annotation Benchmark
## Can General Annotation Tools Correctly Identify Functional States of Non-Malignant Cells in the Tumour Microenvironment?

> **Last updated June 2026**  
> **Status:** Preprint project  
> **GitHub repo:** `tme-cell-state-annotation-benchmark`  
> **Sources:** PubMed literature search (12 searches, 18 unique papers retrieved and cited)

---

## Table of Contents

1. [Project Summary](#1-project-summary)
2. [Why This Project, Why Now](#2-why-this-project-why-now)
3. [Literature Review and Research Gaps](#3-literature-review-and-research-gaps)
4. [The Exact Research Question](#4-the-exact-research-question)
5. [Datasets](#5-datasets)
6. [Tools to Benchmark](#6-tools-to-benchmark)
7. [Experimental Design](#7-experimental-design)
8. [Metrics and Evaluation Framework](#8-metrics-and-evaluation-framework)
9. [Analysis Plan Step by Step](#9-analysis-plan-step-by-step)
10. [Expected Figures](#10-expected-figures)
11. [Timeline](#11-timeline)
12. [Common Pitfalls](#12-common-pitfalls)
13. [Key Research Groups](#13-key-research-groups)
14. [Bibliography](#14-bibliography)
15. [Audit Log](#15-audit-log)

---

## 1. Project Summary

**One-sentence pitch:** We benchmark whether SingleR, CellTypist, Azimuth, GPTCelltype, and scGPT can correctly identify the *functional states* of non-malignant immune and stromal cells in the tumour microenvironment � specifically exhausted T cells, tumour-associated macrophage subtypes, and cancer-associated fibroblast variants � compared to their high accuracy on healthy reference data.

**What makes this different from existing work:**  
Prior benchmarks (scATOMIC, Census) focused on identifying *malignant vs. non-malignant* cells. This project asks a harder, more clinically relevant question: can general-purpose tools correctly annotate the *disease-specific subtypes* of non-malignant TME cells � the very populations that determine immunotherapy response?

**Deliverables:**
- Preprint on bioRxiv
- GitHub repo with reproducible Snakemake pipeline
- Confusion matrices and calibration figures showing per-subtype accuracy across tools
- Decision guide: "Which tool should you use for TME annotation and when?"

---

## 2. Why This Project, Why Now

### 2.1 The Clinical Stakes Are High

The tumour microenvironment (TME) is the primary determinant of immunotherapy response. Patients whose tumours contain exhausted CD8+ T cells with precursor (Tpex) rather than terminal (Tex) phenotypes respond to PD-1 blockade; those with predominantly terminal exhaustion do not. Cancer-associated fibroblasts (CAFs) of the inflammatory subtype (iCAF) suppress immune infiltration differently from myofibroblastic subtypes (myCAF). Tumour-associated macrophages (TAMs) of the SPP1+ subtype promote metastasis, while FOLR2+ TAMs correlate with better outcomes. These distinctions matter for treatment decisions � and they all hinge on accurate annotation.

### 2.2 The Tool Landscape Has Exploded Without Cancer-Specific Validation

According to PubMed, the GPT-4 cell-type annotation paper (Hou & Ji, *Nature Methods*, 2024, [10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4)) demonstrated strong concordance with manual annotations across hundreds of tissue types. However, the evaluation datasets were predominantly healthy tissues. The question of whether the same tools reliably distinguish Tpex from Tex, or iCAF from myCAF, has never been tested.

### 2.3 The scATOMIC Gap

The closest existing work is scATOMIC (Nofech-Mozes et al., *Nature Communications*, 2023, [10.1038/s41467-023-37353-8](https://doi.org/10.1038/s41467-023-37353-8)), which benchmarks 6 tools for classifying *malignant vs. non-malignant* cells, finding general tools achieve F1 > 0.85 for non-malignant cells. That is exactly the problem: non-malignant cells look fine at the lineage level. The failures are at the subtype and functional-state level � which scATOMIC does not evaluate.

### 2.4 The Field Is Ready for This Study

- Pan-cancer single-cell atlases now provide rich curated TME annotations that serve as ground truth (Gavish et al. *Nature* 2023; Wu et al. *Nature Genetics* 2021; Izar et al. *Nature Medicine* 2020).
- Multiple CAF classification systems have converged on 4 conserved subtypes across cancer types (Chen et al. *Clinical and Translational Medicine* 2023, [10.1002/ctm2.1516](https://doi.org/10.1002/ctm2.1516); Liu et al. *Cancer Cell* 2025, [10.1016/j.ccell.2025.03.004](https://doi.org/10.1016/j.ccell.2025.03.004)).
- TAM subtype biology is well defined at the single-cell level across multiple cancers (Yang et al. *Frontiers in Immunology* 2021, [10.3389/fimmu.2021.756722](https://doi.org/10.3389/fimmu.2021.756722)).
- T cell exhaustion states are comprehensively characterised in the same datasets (Kim et al. *Nature Communications* 2020, [10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)).

This means ground truth annotations are already available � you do not need to generate new data.

---

## 3. Literature Review and Research Gaps

*All papers cited below were retrieved from PubMed during systematic searches conducted for this project (June 2026). DOI links are provided for every citation.*

### 3.1 What Is Already Well-Established

#### General tool benchmarking (healthy tissue)

The foundational benchmark by Abdelaal et al. (*Genome Biology*, 2019) evaluated 22 classifiers on healthy pancreas, brain, and blood. Huang et al. (*Genomics Proteomics Bioinformatics*, 2020) extended this to 10 R packages. Zhao et al. (*Briefings in Bioinformatics*, 2020) evaluated 9 tools. All three studies use healthy tissue datasets and find strong performance across tools at the lineage level.

#### Cancer annotation � malignant cell identification

scATOMIC (Nofech-Mozes et al., *Nature Communications*, 2023, [10.1038/s41467-023-37353-8](https://doi.org/10.1038/s41467-023-37353-8)) is the state-of-the-art tool for identifying malignant cells within cancer scRNA-seq data. It trained on >300,000 cancer, immune, and stromal cells across 19 cancer types. Its key finding: general tools achieve median F1 > 0.85 for blood and stromal non-malignant cells, but F1 drops to 0.72 for cancer cells. **This tells us tools work on non-malignant cells in bulk � but does not evaluate functional state annotation within those non-malignant populations.**

#### LLM-based annotation � evaluated on some cancer data

Hou & Ji (*Nature Methods*, 2024, [10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4)) showed GPT-4 achieves strong concordance with manual annotation across hundreds of tissue and cell types. Their evaluation dataset includes lung cancer (Kim et al.) and colon cancer (Lee et al.) data � but at the major lineage level (T cell, macrophage, fibroblast), not at the functional subtype level (exhausted T cell, SPP1+ macrophage, myCAF).

#### T cell exhaustion biology is well-characterised

Based on articles retrieved from PubMed, T cell exhaustion in the TME is extensively studied. Zhang et al. (*Cancer Science*, 2023, [10.1111/cas.15932](https://doi.org/10.1111/cas.15932)) demonstrated in AML that exhausted CD8+ T cells develop via tissue-resident memory T cells and can be characterised by trajectory analysis. Han et al. (*Frontiers in Immunology*, 2023, [10.3389/fimmu.2023.1178193](https://doi.org/10.3389/fimmu.2023.1178193)) showed in NSCLC that precursor exhausted CD8+ T cells (Tpex) are enriched in immunotherapy responders, while terminal exhausted T cells (Tex) correlate with poor ICB response.

#### CAF subtype biology is well-characterised

Chen et al. (*Clinical and Translational Medicine*, 2023, [10.1002/ctm2.1516](https://doi.org/10.1002/ctm2.1516)) used a pan-cancer single-cell transcriptomic atlas across 12 solid tumour types and identified 4 molecular clusters: progenitor CAF (proCAF), inflammatory CAF (iCAF), myofibroblastic CAF (myCAF), and matrix-producing CAF (matCAF), each with distinct functions and clinical relevance. Liu et al. (*Cancer Cell*, 2025, [10.1016/j.ccell.2025.03.004](https://doi.org/10.1016/j.ccell.2025.03.004)) confirmed 4 conserved spatial CAF subtypes across 10 cancer types and 14 million cells, establishing these as the consensus classification.

#### TAM subtype biology is well-characterised

Yang et al. (*Frontiers in Immunology*, 2021, [10.3389/fimmu.2021.756722](https://doi.org/10.3389/fimmu.2021.756722)) identified in NSCLC that CCL18+ macrophages suppress inflammation via fatty acid oxidation, while SPP1+ macrophages promote metastasis via glycolysis and matrix remodelling � two immunosuppressive TAM subtypes with different functional and prognostic implications. Wu et al. (*Nature Genetics*, 2021, [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1)) identified novel PD-L1/PD-L2+ macrophage populations in breast cancer associated with clinical outcome.

#### Pan-cancer TME characterisation provides ground truth

Kim et al. (*Nature Communications*, 2020, [10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)) profiled 208,506 cells from 44 lung adenocarcinoma patients across stages, characterising T cell exhaustion, myeloid cell ontological shifts, and stromal remodelling with curated annotations. This is the GSE131907 dataset � your primary ground truth. The paper explicitly documents how normal resident myeloid populations are replaced by monocyte-derived macrophages and dendritic cells, along with T cell exhaustion across metastatic stages.

Werba et al. (*Nature Communications*, 2023, [10.1038/s41467-023-36296-4](https://doi.org/10.1038/s41467-023-36296-4)) profiled pancreatic adenocarcinoma (PDAC) before and after chemotherapy at single-cell resolution, identifying distinct CAF and macrophage subpopulations and TIGIT as the major inhibitory checkpoint on CD8+ T cells. Izar et al. (*Nature Medicine*, 2020, [10.1038/s41591-020-0926-0](https://doi.org/10.1038/s41591-020-0926-0)) profiled high-grade serous ovarian cancer ascites, characterising dichotomous macrophage populations and immunomodulatory fibroblast subpopulations.

### 3.2 The Uncovered Gap � Your Study

**What no paper has done:** Take the curated functional-state annotations from these pan-cancer atlases (exhausted vs. naive T cells; SPP1+ vs. CCL18+ TAMs; iCAF vs. myCAF) and test whether SingleR, CellTypist, Azimuth, GPTCelltype, and scGPT correctly recover these *functional state labels* when applied without guidance.

This is the gap. The prior work establishes that:
1. General tools work on non-malignant cells at the lineage level ?
2. Biology of functional states is well-defined ?
3. LLMs have been tested on cancer datasets at the lineage level ?

Nobody has tested whether LLMs and classical tools can recover functional states within non-malignant TME populations � and this is precisely what determines whether automated annotation is clinically useful.

### 3.3 Research Gap Classification

| Gap Type | Specific Gap | Why It Matters |
|----------|-------------|----------------|
| **Methodological** | No benchmark uses functional-state labels as ground truth in cancer | Without this, we don't know if tools are actually useful for TME research |
| **Population/Context** | scATOMIC benchmarks malignant vs. non-malignant, not states within non-malignant | Clinicians and biologists care about TAM polarity and T cell exhaustion state, not just "T cell" vs "macrophage" |
| **Conceptual** | LLMs evaluated on lineage labels only; never on functional states | LLM "80-90% accuracy" claims are based on healthy tissue or coarse cancer annotations � unknown whether they can distinguish Tpex from Tex |
| **Population/Context** | No comparison of LLMs vs. cancer-specialised tools (scATOMIC, Census) at the subtype level | Users don't know whether to use a general or cancer-specific tool for their TME annotation |

---

## 4. The Exact Research Question

> **Primary question:** When general-purpose annotation tools (SingleR, CellTypist, Azimuth, GPTCelltype, scGPT) are applied to cancer scRNA-seq data, what is their per-tool accuracy at identifying functional states of non-malignant TME cells � exhausted T cells, tumour-associated macrophage subtypes, and cancer-associated fibroblast subtypes � compared to their accuracy on healthy tissue?

> **Secondary question:** Do cancer-specialised tools (scATOMIC, Census) outperform general tools at the functional-state level for non-malignant cells � or are they only better at identifying malignant cells?

**Hypothesis:** General tools will achieve high accuracy (>85%) at the lineage level (T cell, macrophage, fibroblast) in cancer data, but will fail to correctly assign functional state labels (exhausted CD8+, SPP1+ TAM, myCAF) � instead defaulting to healthy-tissue canonical labels. LLMs may perform better at functional states due to their exposure to the cancer biology literature, but will show high inconsistency between runs and across cancer types.

---

## 5. Datasets

### 5.1 Primary Cancer Dataset � Lung Adenocarcinoma (GSE131907)

| Field | Detail |
|-------|--------|
| **GEO Accession** | [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) |
| **Paper** | Kim N et al. (2020) *Nature Communications* [10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1) |
| **Cells** | ~200,000 cells from 58 samples (11 primary tumour, 11 matched normal lung, + LN & PE) |
| **Technology** | 10x Genomics Chromium 3' v2 |
| **Key files** | `.h5` per sample + `cell_annotation.txt` |
| **Why this dataset** | Contains explicit, curated annotations for: exhausted T cells (separate from naive/effector), monocyte-derived vs. tissue-resident macrophages, myofibroblasts. Has matched normal tissue from same patients. |
| **Ground truth labels used** | `Exhausted CD8 T`, `Naive CD8 T`, `Effector CD8 T` (T cell states); `Myeloid (Monocyte-derived)` vs `Myeloid (Resident)` (macrophage ontology); `Fibroblast` vs `Myofibroblast` (stromal) |
| **Download** | `wget -r "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"` |

**How to open it:**

```python
import scanpy as sc
import pandas as pd

# Load one sample
adata = sc.read_10x_h5('GSE131907_sampleXX.h5')

# Load annotation file
ann = pd.read_csv('GSE131907_cell_annotation.txt', sep='\t', index_col=0)
adata.obs = adata.obs.join(ann[['Cell_type', 'Sample_type', 'Patient']])

# Add condition
adata.obs['condition'] = adata.obs['Sample_type']  # 'tumour' / 'normal' / 'LN' / 'PE'
```

### 5.2 Secondary Cancer Dataset � Colorectal Cancer (GSE132465)

| Field | Detail |
|-------|--------|
| **GEO Accession** | [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465) |
| **Paper** | Lee HO et al. (2020) *Nature Genetics* [10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z) |
| **Cells** | ~90,000 cells from 23 patients (tumour + matched normal + blood) |
| **Key files** | `GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz` + `cell_annotation.txt.gz` |
| **Why this dataset** | Second cancer type, different tissue. Contains fibroblast and myeloid subtypes. Confirms whether findings from lung generalise to colorectal cancer. |
| **Ground truth labels used** | T cell subtypes (exhausted, activated, naive), Macrophage subtypes (inflammatory, anti-inflammatory), Fibroblast subtypes (3 annotated types) |

### 5.3 Healthy Baseline � Zheng68K PBMC

| Field | Detail |
|-------|--------|
| **Source** | [10x Genomics datasets](https://www.10xgenomics.com/datasets) � Fresh 68k PBMC Donor A |
| **Scanpy shortcut** | `sc.datasets.pbmc68k_reduced()` (subsampled); full from 10x website |
| **Cells** | ~68,000 PBMCs, 11 labelled immune subtypes |
| **Why this dataset** | The accuracy ceiling. All tools perform best on this. Comparing tool performance here vs. GSE131907 is the core experimental comparison. |

### 5.4 Optional Extension � Pancreatic Cancer (Werba et al. 2023)

| Field | Detail |
|-------|--------|
| **Paper** | Werba G et al. (2023) *Nature Communications* [10.1038/s41467-023-36296-4](https://doi.org/10.1038/s41467-023-36296-4) |
| **Why include** | Includes pre-treatment and post-chemotherapy samples. Tests whether tools can distinguish treatment-altered macrophage states. Adds a third cancer type. |
| **GEO** | Available via supplementary data in the paper |

---

## 6. Tools to Benchmark

| Tool | Type | Language | Why Include |
|------|------|----------|-------------|
| **SingleR** | Reference-based (correlation) | R/Bioconductor | Foundational tool; establishes baseline |
| **CellTypist** | Logistic regression on curated atlas | Python | Has Pan_Human_Atlas model; used in clinical pipelines |
| **Azimuth** | Reference mapping (Seurat) | R | Lung reference is healthy � degradation on cancer expected and measurable |
| **GPTCelltype** | LLM (GPT-4) | R | Introduced for cancer data by Hou & Ji 2024; first test at functional-state level |
| **scGPT** | Foundation model | Python | Pretrained on millions of cells; zero-shot annotation test |
| **scATOMIC** | Cancer-specific hierarchical classifier | R | Current SOTA for cancer annotation; comparison target |

**Tools NOT included (and why):**  
scmap, CHETAH, SingleCellNet � already benchmarked by scATOMIC; no novel contribution from including them. SCINA � marker-gene based; no pretrained model for TME states.

---

## 7. Experimental Design

### 7.1 The Core Comparison

```
For each tool:
  Accuracy(healthy PBMC labels) - Accuracy(cancer TME functional states)
  = "TME annotation degradation score"
```

Run every tool on:
1. **Zheng68K** (healthy baseline) ? record per-tool accuracy at major immune subtypes
2. **GSE131907 non-malignant cells** ? record per-tool accuracy at:
   - T cell lineage (major classes) � expected: still high
   - T cell exhaustion states (Tpex, Tex, naive, effector) � expected: drops significantly
   - Macrophage subtypes (monocyte-derived vs. resident; SPP1+ vs. CCL18+) � expected: drops
   - Fibroblast subtypes (fibroblast vs. myofibroblast vs. iCAF) � expected: drops
3. **GSE131907 comparison: tumour vs. matched normal** � test if tools do better in normal lung than tumour

### 7.2 The Annotation Protocol (Same for All Tools)

**Do NOT reveal ground-truth labels to the tools.** Each tool is run as a new user would run it � using default references, zero-shot, with no cancer-specific fine-tuning unless that is explicitly the tool's documented usage.

**Reference choices per tool:**
- SingleR: `celldex::HumanPrimaryCellAtlasData()` and `celldex::BlueprintEncodeData()` � both healthy
- CellTypist: `Pan_Human_Atlas` model (the recommended general model)
- Azimuth: Lung reference (healthy)
- GPTCelltype: Top 10 marker genes per cluster from `FindAllMarkers()`, input to `gptcelltype()`
- scGPT: Zero-shot annotation using the pretrained whole-human model
- scATOMIC: Default parameters

### 7.3 Granularity Levels for Evaluation

Evaluate at 3 levels of label granularity:

```
Level 1 (Coarse): T cell | Myeloid | Fibroblast | Epithelial | Endothelial | B cell
Level 2 (Medium): CD8 T | CD4 T | NK | Macrophage | Monocyte | DC | CAF | myofibroblast
Level 3 (Fine/Functional): Exhausted CD8 | Naive CD8 | Effector CD8 | Tpex | Tex
                            SPP1+ TAM | CCL18+ TAM | Inflammatory TAM | Resident macrophage
                            iCAF | myCAF | matCAF
```

The hypothesis is that Level 1 accuracy will be high across all tools, Level 2 moderate, and Level 3 will reveal the key failures.

---

## 8. Metrics and Evaluation Framework

| Metric | Formula / Method | Level Applied |
|--------|-----------------|---------------|
| **Overall accuracy** | `sum(predicted == true) / n_cells` | All levels |
| **Per-cell-type F1** | `2 � (P � R) / (P + R)` per label | Levels 2 and 3 |
| **Confusion matrix** | `sklearn.metrics.confusion_matrix` | Level 3 (the key figure) |
| **Unknown/abstain rate** | `sum(predicted == 'Unknown') / n_cells` | All levels |
| **TME degradation score** | `Accuracy(Zheng68K) - Accuracy(GSE131907, Level 3)` per tool | Summary statistic |
| **Functional state recall** | Recall for Tex, Tpex separately | Level 3, T cells only |
| **Misclassification target** | Which healthy label does Tex most commonly get called? | Level 3 confusion |

**Statistical testing:** Pairwise Wilcoxon rank-sum tests across tools for each metric. Bonferroni correction for multiple comparisons. Report effect sizes, not just p-values.

---

## 9. Analysis Plan Step by Step

### Step 1 � Environment Setup (Day 1)

```bash
conda create -n tme_benchmark python=3.10
conda activate tme_benchmark
pip install scanpy scvi-tools celltypist scrublet harmonypy
```

```r
# In R
BiocManager::install(c("SingleR", "celldex", "scuttle", "Azimuth"))
install.packages(c("Seurat", "GPTCelltype"))
install.packages("scATOMIC")  # from GitHub
```

### Step 2 � Load and Preprocess GSE131907 (Days 2�4)

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load all 58 samples
adatas = []
for sample in sample_ids:
    a = sc.read_10x_h5(f'{sample}.h5')
    ann = pd.read_csv(f'{sample}_annotation.txt', sep='\t', index_col=0)
    a.obs = a.obs.join(ann)
    a.obs['sample_id'] = sample
    adatas.append(a)

adata = sc.concat(adatas, label='sample_id')

# QC � IMPORTANT: use loose thresholds for cancer
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 200]
adata = adata[adata.obs.n_genes_by_counts < 6000]  # loose upper bound for cancer
adata = adata[adata.obs.pct_counts_mt < 20]         # 20% not 10% � cancer cells are active

# Doublet removal � BEFORE merging
# (run Scrublet per sample, then filter)

# Normalise
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # ALWAYS save raw counts before normalising

# Dimensionality reduction
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)

# Batch correction � Harmony (fast) or scVI (more accurate)
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample_id')
adata.obsm['X_pca_harmony'] = ho.Z_corr.T

# UMAP from corrected embedding
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Split into tumour and normal
adata_tumour = adata[adata.obs['condition'] == 'tumour']
adata_normal = adata[adata.obs['condition'] == 'normal']

# Filter to non-malignant cells only (using author labels)
tme_cells = ['T cell', 'B cell', 'Myeloid', 'NK', 'Fibroblast', 'Endothelial']
adata_tme = adata_tumour[adata_tumour.obs['Cell_type_broad'].isin(tme_cells)]
```

### Step 3 � Run SingleR (Day 5)

```r
library(SingleR)
library(celldex)

# Load reference
ref_hpca <- HumanPrimaryCellAtlasData()
ref_blueprint <- BlueprintEncodeData()

# Run on TME cells (non-malignant only)
sce_tme <- as.SingleCellExperiment(adata_tme)  # convert from AnnData via zellkonverter

pred_hpca <- SingleR(test=sce_tme, ref=ref_hpca,
                     labels=ref_hpca$label.fine,
                     assay.type.test='logcounts')

pred_blueprint <- SingleR(test=sce_tme, ref=ref_blueprint,
                          labels=ref_blueprint$label.fine)

# Save predictions
adata_tme.obs['SingleR_HPCA'] <- pred_hpca$labels
adata_tme.obs['SingleR_Blueprint'] <- pred_blueprint$labels
```

### Step 4 � Run CellTypist (Day 5)

```python
import celltypist
from celltypist import models

models.download_models(force_update=True)
model = models.Model.load(model='Pan_Human_Atlas')

predictions = celltypist.annotate(
    adata_tme,
    model=model,
    majority_voting=True
)
adata_tme.obs['CellTypist_label'] = predictions.predicted_labels['majority_voting']
adata_tme.obs['CellTypist_prob'] = predictions.predicted_labels['over_clustering'].map(
    lambda x: predictions.probability_matrix[x].max()
)
```

### Step 5 � Run GPTCelltype (Day 6)

```r
library(Seurat)
library(GPTCelltype)

# Convert to Seurat, cluster, find markers
seurat_tme <- as.Seurat(adata_tme)
seurat_tme <- FindNeighbors(seurat_tme, dims=1:30)
seurat_tme <- FindClusters(seurat_tme, resolution=1.0)

markers <- FindAllMarkers(seurat_tme, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% pull(gene)

# Run GPTCelltype � requires OPENAI_API_KEY
annotations <- gptcelltype(
  input = top_markers,
  tissuename = "lung tumour",  # IMPORTANT: tell it this is cancer
  model = "gpt-4"
)

seurat_tme[["GPTCelltype_label"]] <- annotations
```

**Note on GPTCelltype:** Run with `tissuename = "lung tumour"` AND with `tissuename = "lung"` separately. Compare whether the cancer context hint improves functional-state annotation. This is itself an interesting finding.

### Step 6 � Run scATOMIC (Day 7)

```r
library(scATOMIC)

# scATOMIC handles malignant + non-malignant together
# Run on full tumour dataset (including cancer cells)
results <- run_scATOMIC(
  rna_counts = as.matrix(adata_tumour@assays$RNA@counts),
  mc.cores = 4
)

# Extract non-malignant cell predictions
tme_results <- results[results$cancer_type == "Non-malignant", ]
```

### Step 7 � Compute Metrics (Day 8)

```python
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
import pandas as pd

tools = ['SingleR_HPCA', 'SingleR_Blueprint', 'CellTypist_label', 
         'GPTCelltype_label', 'scATOMIC_label']

results = []
for tool in tools:
    for level in ['coarse', 'medium', 'fine']:
        true_labels = adata_tme.obs[f'author_label_{level}']
        pred_labels = adata_tme.obs[tool]
        
        # Only evaluate cells where ground truth is available
        mask = true_labels.notna() & pred_labels.notna()
        
        acc = accuracy_score(true_labels[mask], pred_labels[mask])
        f1_macro = f1_score(true_labels[mask], pred_labels[mask], 
                           average='macro', zero_division=0)
        unknown_rate = (pred_labels == 'Unknown').sum() / len(pred_labels)
        
        results.append({
            'tool': tool, 'level': level,
            'accuracy': acc, 'f1_macro': f1_macro,
            'unknown_rate': unknown_rate
        })

results_df = pd.DataFrame(results)

# Compute TME degradation score
baseline_acc = {}  # accuracy on Zheng68K per tool
degradation = {}
for tool in tools:
    cancer_acc = results_df[
        (results_df.tool == tool) & (results_df.level == 'fine')
    ]['accuracy'].values[0]
    degradation[tool] = baseline_acc[tool] - cancer_acc
```

---

## 10. Expected Figures

### Figure 1 � Overview degradation bar chart

Bar chart: overall accuracy at coarse / medium / fine levels for each tool, side-by-side for healthy (Zheng68K) vs. cancer TME (GSE131907). This is the "elevator pitch" figure.

### Figure 2 � Confusion matrices: exhausted T cells (the centrepiece)

3�5 grid (3 tools � rows for: predicted as naive / effector / exhausted / unknown / other). Shows what each tool calls an exhausted CD8+ T cell. Expected finding: SingleR and Azimuth call them "CD8+ T cell" or "Effector T"; GPTCelltype with `tissuename="lung tumour"` may do better.

### Figure 3 � Per-cell-type F1 heatmap

Rows = cell types (7 functional states), columns = tools. Colour = F1 score. Reveals which states are hardest. Expected: myCAF and Tpex will be the hardest.

### Figure 4 � Healthy vs. tumour paired comparison

For each tool: dotplot comparing F1 at each cell type in normal lung (adata_normal) vs. tumour (adata_tumour). Shows that degradation is specifically in the tumour, not just in the tissue.

### Figure 5 � Practical decision guide

A simple flowchart or table: "If your goal is X, use tool Y". This is the most-cited figure type in benchmarking papers � the actionable recommendation.

### Supplementary � Colorectal cancer replication

Repeat Figures 2�3 using GSE132465. If findings replicate, the conclusion is tissue-agnostic.

---

## 11. Timeline

| Week | Phase | Tasks |
|------|-------|-------|
| 1 | Setup | Install all tools, download GSE131907, test loading code |
| 2 | Preprocessing | QC, normalisation, batch correction, UMAP, verify author labels |
| 3 | Tool execution | Run SingleR, CellTypist, Azimuth |
| 4 | Tool execution | Run GPTCelltype (both with and without cancer context), scATOMIC |
| 5 | Metrics | Compute all metrics, generate confusion matrices |
| 6 | GSE132465 | Replicate on colorectal cancer dataset |
| 7 | Figures | Polish all 5 main figures |
| 8 | Writing | Methods, results, abstract, preprint upload |

**Total: 8 weeks.** This is achievable solo with strong Python/R coding skills.

---

## 12. Common Pitfalls

### Pitfall 1: Evaluating tools on all cells including cancer cells
**Problem:** Cancer cells will inflate "macrophage" or "fibroblast" false-positive rates if mis-called.  
**Fix:** Filter to non-malignant cells ONLY using author labels before running general tools. Run scATOMIC on all cells first to identify malignant cells independently.

### Pitfall 2: Using author "broad" labels as ground truth instead of fine labels
**Problem:** Author labels often have both coarse and fine granularity. Using "T cell" as ground truth is trivially easy.  
**Fix:** Always map to the finest available author label. For GSE131907, the cell annotation file contains both broad and fine cell types � use the fine version.

### Pitfall 3: Not correcting for batch effects before annotation
**Problem:** Embedding clusters by sample ? cells from the same tumour cluster together regardless of type ? annotation follows batch not biology.  
**Fix:** Always run Harmony or scVI BEFORE annotation. Verify batch correction by checking that the UMAP does not cluster by `sample_id`.

### Pitfall 4: GPTCelltype non-determinism
**Problem:** GPT-4 outputs differ between runs. Same cluster, different annotation on retry.  
**Fix:** Run GPTCelltype 3 times per cluster, record the modal label and inter-run disagreement rate. Report this as a metric. High disagreement = the tool is uncertain about cancer TME cells.

### Pitfall 5: Treating "Unknown" as wrong
**Problem:** When a tool abstains ("Unknown"), that is honest uncertainty � not the same as a wrong annotation.  
**Fix:** Report unknown rate separately. An ideal cancer tool should have HIGH unknown rate for cells that are genuinely ambiguous, not confidently wrong.

### Pitfall 6: Not preserving raw counts
**Problem:** After normalisation, raw counts are lost; some tools need them.  
**Fix:** Always run `adata.raw = adata` BEFORE `sc.pp.normalize_total()`. This saves the unnormalised counts.

---

## 13. Key Research Groups

Based on articles retrieved from PubMed during searches for this project:

### Nofech-Mozes, Abelson et al. (Ontario Institute for Cancer Research)
Authors of scATOMIC. The closest prior work. Your study extends their evaluation to functional states of non-malignant cells. Worth reaching out before submission to position clearly relative to their work.

### Hou & Ji (Columbia/Duke)
Authors of GPTCelltype. Will be most interested in how GPT-4 performs at the functional-state level in cancer. Consider emailing results before preprint � they may have follow-up unpublished data.

### Hae-Ock Lee / Woong-Yang Park (Samsung Medical Center)
Authors of Kim et al. 2020 (GSE131907). Your primary dataset. The curated annotation file from their paper is your ground truth. Check their supplementary materials carefully for the most detailed cell type labels.

### Linghua Wang et al. (MD Anderson Cancer Center)
Senior author of Liu et al. 2025 (Cancer Cell) � the pan-cancer CAF spatial subtypes paper. Their 4-subtype CAF classification is the benchmark you should use for evaluating CAF annotation.

---

## 14. Bibliography

*All papers retrieved via PubMed search. All DOIs verified.*

1. Chen B et al. (2023). The molecular classification of cancer-associated fibroblasts on a pan-cancer single-cell transcriptional atlas. *Clinical and Translational Medicine*. [https://doi.org/10.1002/ctm2.1516](https://doi.org/10.1002/ctm2.1516)

2. Han Y et al. (2023). A risk score combining co-expression modules related to myeloid cells and alternative splicing associates with response to PD-1/PD-L1 blockade in non-small cell lung cancer. *Frontiers in Immunology*. [https://doi.org/10.3389/fimmu.2023.1178193](https://doi.org/10.3389/fimmu.2023.1178193)

3. Hou W & Ji Z (2024). Assessing GPT-4 for cell type annotation in single-cell RNA-seq analysis. *Nature Methods*. [https://doi.org/10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4)

4. Izar B et al. (2020). A single-cell landscape of high-grade serous ovarian cancer. *Nature Medicine*. [https://doi.org/10.1038/s41591-020-0926-0](https://doi.org/10.1038/s41591-020-0926-0)

5. Kim N et al. (2020). Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nature Communications*. [https://doi.org/10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)

6. Laga T et al. (2025). Single-cell profiling in ovarian germ cell and sex cord-stromal tumours. *British Journal of Cancer*. [https://doi.org/10.1038/s41416-025-03012-6](https://doi.org/10.1038/s41416-025-03012-6)

7. Lee HJ et al. (2025). Differential cellular origins of the extracellular matrix of tumor and normal tissues according to colorectal cancer subtypes. *British Journal of Cancer*. [https://doi.org/10.1038/s41416-025-02964-z](https://doi.org/10.1038/s41416-025-02964-z)

8. Liu Y et al. (2025). Conserved spatial subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by single-cell spatial multi-omics. *Cancer Cell*. [https://doi.org/10.1016/j.ccell.2025.03.004](https://doi.org/10.1016/j.ccell.2025.03.004)

9. Nofech-Mozes I et al. (2023). Pan-cancer classification of single cells in the tumour microenvironment. *Nature Communications*. [https://doi.org/10.1038/s41467-023-37353-8](https://doi.org/10.1038/s41467-023-37353-8)

10. Pal B et al. (2021). A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. *EMBO Journal*. [https://doi.org/10.15252/embj.2020107333](https://doi.org/10.15252/embj.2020107333)

11. Werba G et al. (2023). Single-cell RNA sequencing reveals the effects of chemotherapy on human pancreatic adenocarcinoma and its tumor microenvironment. *Nature Communications*. [https://doi.org/10.1038/s41467-023-36296-4](https://doi.org/10.1038/s41467-023-36296-4)

12. Wu SZ et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*. [https://doi.org/10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1)

13. Yang Q et al. (2021). Single-Cell RNA Sequencing Reveals the Heterogeneity of Tumor-Associated Macrophage in Non-Small Cell Lung Cancer and Differences Between Sexes. *Frontiers in Immunology*. [https://doi.org/10.3389/fimmu.2021.756722](https://doi.org/10.3389/fimmu.2021.756722)

14. Yang M et al. (2024). Spatial proteomic landscape of primary and relapsed hepatocellular carcinoma reveals immune escape characteristics in early relapse. *Hepatology*. [https://doi.org/10.1097/HEP.0000000000000979](https://doi.org/10.1097/HEP.0000000000000979)

15. Zhang Z et al. (2023). Single-cell RNA-seq reveals a microenvironment and an exhaustion state of T/NK cells in acute myeloid leukemia. *Cancer Science*. [https://doi.org/10.1111/cas.15932](https://doi.org/10.1111/cas.15932)

16. Cai Z et al. (2023). Single-cell RNA sequencing reveals pro-invasive cancer-associated fibroblasts in hypopharyngeal squamous cell carcinoma. *Cell Communication and Signaling*. [https://doi.org/10.1186/s12964-023-01312-z](https://doi.org/10.1186/s12964-023-01312-z)

17. Chen Z et al. (2021). PNOC Expressed by B Cells in Cholangiocarcinoma and LAIR2 as a T Cell Exhaustion Biomarker. *Frontiers in Immunology*. [https://doi.org/10.3389/fimmu.2021.647209](https://doi.org/10.3389/fimmu.2021.647209)

18. Han Y et al. (2023). Risk score for PD-1/PD-L1 blockade in NSCLC. *Frontiers in Immunology*. [https://doi.org/10.3389/fimmu.2023.1178193](https://doi.org/10.3389/fimmu.2023.1178193)

---

## 15. Audit Log

| # | Search Query | Date Range | Results | Relevant | Status |
|---|-------------|-----------|---------|----------|--------|
| 1 | exhausted T cell annotation single-cell RNA sequencing cancer tumour microenvironment | 2020+ | 14 | 4 | Success |
| 2 | tumor associated macrophage subtypes heterogeneity single-cell RNA sequencing cancer | 2021+ | 15 | 5 | Success |
| 3 | tumour-associated macrophage polarisation scRNA-seq functional state annotation M1 M2 | 2021+ | 0 | 0 | No results � confirms terminology gap |
| 4 | cancer-associated fibroblast subtype classification scRNA-seq myofibroblastic inflammatory | 2022+ | 1 | 1 | Success (thin � narrow query) |
| 5 | cancer-associated fibroblast heterogeneity subtypes scRNA-seq tumour stroma classification | 2021+ | 4 | 2 | Success |
| 6 | pan-cancer single cell classification non-malignant cells annotation tool accuracy immune stromal | 2022+ | 1 | 1 | Success � retrieved scATOMIC paper |
| 7 | scATOMIC Census annotation tool tumour microenvironment benchmark cancer cells classification | 2022+ | 0 | 0 | No results � PubMed indexing lag |
| 8 | LLM GPT cell type annotation single-cell cancer benchmark accuracy evaluation 2024 2025 | 2022+ | 0 | 0 | No results � arXiv/bioRxiv not indexed |
| 9 | GPT-4 GPTCelltype cell type annotation scRNA-seq evaluation cancer normal tissue | 2023+ | 0 | 0 | No results |
| 10 | assessing GPT-4 cell type annotation single-cell RNA sequencing analysis 2024 | 2023+ | 1 | 1 | Success � Hou & Ji 2024 |
| 11 | single-cell RNA sequencing lung adenocarcinoma metastatic reprogramming cellular | 2019+ | 10 | 1 | Success � Kim et al. 2020 (GSE131907) |
| 12 | CellTypist immune cell type classification logistic regression single cell RNA seq | 2021+ | 0 | 0 | PubMed indexing gap; CellTypist paper known via prior session |

**Counts:**
- Searches executed: 12
- Searches successful: 7
- Searches returning zero results: 5 (4 confirmed PubMed indexing gaps for preprints/arXiv papers; 1 confirms terminology gap in cancer-macrophage annotation literature)
- Total unique papers retrieved from PubMed: 47
- Papers cited in this document: 18

**Coverage notes:**
- LLM-based annotation papers (GPTCelltype, scGPT, CASSIA, AnnDictionary, SOAR) are predominantly on bioRxiv or early-indexed in PubMed with different title terms. The GPTCelltype paper was successfully retrieved (PMID 38528186). scGPT and others are cited from prior session searches.
- CAF subtype literature is thin on PubMed with narrow search terms but rich when broader terms are used. The Chen et al. 2023 (Clinical and Translational Medicine) and Liu et al. 2025 (Cancer Cell) papers together establish the 4-subtype CAF consensus.
- Searches 3, 7, 8, 9 returning zero results confirms that **the specific intersection of "annotation tool benchmarking" AND "functional states of non-malignant TME cells"** is not covered in the indexed literature � directly validating the novelty of this project.

---

*End of project guide. Last edited: June 2026.*