# scPharm <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/scPharm)](https://zaoqu-liu.r-universe.dev/scPharm)
[![GitHub release](https://img.shields.io/github/v/release/Zaoqu-Liu/scPharm)](https://github.com/Zaoqu-Liu/scPharm/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**Identification of Pharmacological Subpopulations of Single Cells for Precision Medicine in Cancers**

## Overview

Intratumoral heterogeneity poses a significant challenge to effective cancer treatment, as diverse cellular subpopulations within a tumor may exhibit distinct drug response profiles. **scPharm** is a computational framework that integrates single-cell RNA sequencing (scRNA-seq) data with pharmacogenomics profiles from the Genomics of Drug Sensitivity in Cancer (GDSC2) database to dissect therapeutic heterogeneity at single-cell resolution.

The framework employs Multiple Correspondence Analysis (MCA) for dimensionality reduction and Gene Set Enrichment Analysis (GSEA) to quantify drug sensitivity signatures, enabling the classification of individual cells into pharmacologically distinct subpopulations.

### Key Capabilities

| Module | Function | Description |
|--------|----------|-------------|
| **Cell Classification** | `scPharmIdentify()` | Stratifies cells into drug-sensitive and drug-resistant subpopulations based on pharmacogenomics signatures |
| **Drug Prioritization** | `scPharmDr()` | Ranks therapeutic agents by tumor cell sensitivity profiles using a composite scoring metric |
| **Toxicity Prediction** | `scPharmDse()` | Estimates potential drug side effects by evaluating sensitivity in non-malignant cells |
| **Combination Therapy** | `scPharmCombo()` | Identifies synergistic drug combinations targeting complementary resistant subpopulations |

### Methodological Framework

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           scPharm Workflow                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐  │
│  │  scRNA-seq  │───▶│     MCA     │───▶│    GSEA     │───▶│   Cell      │  │
│  │    Data     │    │  Embedding  │    │  Scoring    │    │   Labels    │  │
│  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘  │
│         │                                                        │          │
│         ▼                                                        ▼          │
│  ┌─────────────┐                                          ┌─────────────┐  │
│  │   GDSC2     │                                          │    Drug     │  │
│  │  Profiles   │                                          │  Ranking    │  │
│  └─────────────┘                                          └─────────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Installation

### From R-universe (Recommended)

```r
install.packages("scPharm", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("Zaoqu-Liu/scPharm")
```

### System Requirements

- **R** ≥ 4.1
- **Seurat** ≥ 4.0.1 (compatible with both V4 and V5)
- Platform: Linux, macOS, Windows

## Quick Start

```r
library(scPharm)
library(Seurat)

# Load example Seurat object
# seurat_obj <- readRDS("your_data.rds")

# Step 1: Identify pharmacological subpopulations
result <- scPharmIdentify(
  seurat_obj,
  type = "tissue",        # "tissue" or "cellline"
  cancer = "LUAD",        # TCGA cancer type abbreviation
  drug = "Erlotinib",     # Drug name(s) or "all"
  nmcs = 50,              # Number of MCA components
  nfeatures = 200,        # Number of features for signature
  cores = 4               # Parallel cores
)

# Step 2: Compute drug prioritization scores
dr_scores <- scPharmDr(result)
head(dr_scores)

# Step 3: Predict drug side effects (tissue samples only)
dse_scores <- scPharmDse(result)

# Step 4: Identify potential drug combinations
combos <- scPharmCombo(result, dr_scores, topN = 5)
```

## Detailed Usage

### Sample Type Selection

| Type | Description | CNV Detection | Use Case |
|------|-------------|---------------|----------|
| `"tissue"` | Primary tumor tissue | Automatic (via CopyKAT algorithm) | Patient-derived samples with stromal/immune infiltration |
| `"cellline"` | Cancer cell lines | Not required | Pure tumor cell populations |

### Cancer Type Specification

scPharm supports all TCGA cancer type abbreviations:

```r
# Single cancer type
result <- scPharmIdentify(seurat_obj, cancer = "BRCA", ...)

# Multiple cancer types
result <- scPharmIdentify(seurat_obj, cancer = c("BRCA", "LUAD"), ...)

# Pan-cancer analysis
result <- scPharmIdentify(seurat_obj, cancer = "pan", ...)
```

**Supported cancer types:** `BRCA`, `LUAD`, `LUSC`, `COAD`, `STAD`, `LIHC`, `KIRC`, `OV`, `PAAD`, `GBM`, `SKCM`, `HNSC`, `BLCA`, `PRAD`, `UCEC`, `ESCA`, `THCA`, and more.

### Pre-annotated Tumor Cells

For samples with known tumor cell annotations:

```r
# Provide tumor cell barcodes directly (skips CNV-based detection)
tumor_barcodes <- c("AAACCTGAGAAACCAT-1", "AAACCTGAGATCACGG-1", ...)
result <- scPharmIdentify(
  seurat_obj,
  type = "tissue",
  tumor.cells = tumor_barcodes,
  ...
)
```

### Threshold Calibration

Generate empirical thresholds from normal tissue samples:

```r
# Using normal tissue reference
thresholds <- scPharmGenNullDist(
  normal_seurat_obj,
  cancer = "LUAD",
  drug = "Erlotinib"
)

# Apply calibrated thresholds
result <- scPharmIdentify(
  seurat_obj,
  threshold.s = thresholds$threshold.s,
  threshold.r = thresholds$threshold.r,
  ...
)
```

## Output Structure

### scPharmIdentify Output

The function returns a Seurat object with additional metadata columns:

| Column | Description |
|--------|-------------|
| `cell.label` | Cell type annotation (`"tumor"` or `"adjacent"`) |
| `scPharm_label_<drug>` | Drug response label (`"sensitive"`, `"resistant"`, or `"other"`) |
| `scPharm_nes_<drug>` | Normalized Enrichment Score (NES) for drug sensitivity |

### scPharmDr Output

| Column | Description |
|--------|-------------|
| `DRUG_ID` | GDSC drug identifier |
| `DRUG_NAME` | Drug name |
| `Dr` | Drug prioritization score (lower = more effective) |

## Computational Considerations

- **Memory**: Large datasets (>50,000 cells) may require 16+ GB RAM
- **Runtime**: ~5-10 minutes per drug for 10,000 cells (single core)
- **Parallelization**: Use `cores` parameter for multi-threading

## Citation

If you use scPharm in your research, please cite:

> Tian P, Zheng J, Wang H, Liu Z. **scPharm: Identification of pharmacological subpopulations of single cells for precision medicine in cancers.** *Briefings in Bioinformatics*, 2024.

## Related Resources

- [GDSC Database](https://www.cancerrxgene.org/) - Genomics of Drug Sensitivity in Cancer
- [Seurat](https://satijalab.org/seurat/) - Single-cell analysis toolkit
- [CelliD](https://github.com/RausellLab/CelliD) - MCA-based cell identity analysis

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Maintainer**: Zaoqu Liu (liuzaoqu@163.com)
- **Issues**: [GitHub Issues](https://github.com/Zaoqu-Liu/scPharm/issues)
- **Homepage**: [https://github.com/Zaoqu-Liu/scPharm](https://github.com/Zaoqu-Liu/scPharm)
