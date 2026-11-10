#' Identify Pharmacological Cell Subpopulations
#'
#' Classify single cells into drug-sensitive, drug-resistant, or other subpopulations
#' based on pharmacogenomics profiles from the GDSC2 database.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param type Data source type. Either \code{"tissue"} for tumor tissue samples or
#'   \code{"cellline"} for cell line samples. When \code{"tissue"}, the function
#'   identifies tumor vs adjacent normal cells using CNV analysis.
#' @param cancer TCGA cancer type(s). A character string or vector specifying cancer
#'   type(s) (e.g., \code{"BRCA"}, \code{c("LUAD", "LUSC")}). Use \code{"pan"} for
#'   pan-cancer analysis.
#' @param drug Drug name to analyze. If \code{NULL} (default), all drugs from
#'   GDSC2 project will be analyzed.
#' @param nmcs Number of MCA components to compute. Default: 50.
#' @param nfeatures Number of genes for cell identity signature. Default: 200.
#' @param cores Number of CPU cores for parallel processing. Default: 1.
#' @param features Character vector of gene names to use. If \code{NULL},
#'   all features are used.
#' @param slot Slot name for Seurat V4 data access. Default: \code{"data"}.
#' @param layer Layer name for Seurat V5 data access. If \code{NULL}, uses \code{slot} value.
#' @param assay Assay to use. Default: \code{"RNA"}.
#' @param threshold.s Threshold for labeling sensitive cells. Cells with NES below
#'   this value are classified as sensitive. Default: -1.751302.
#' @param threshold.r Threshold for labeling resistant cells. Cells with NES above
#'   this value are classified as resistant. Default: 1.518551.
#' @param tumor.cells Character vector of known tumor cell barcodes. If provided
#'   when \code{type="tissue"}, CNV-based detection is skipped.
#' @param normal.cells Character vector of known normal cell barcodes. Used as
#'   reference for CNV analysis when provided.
#' @param bulkdata Bulk RNA-seq data for cell lines. If \code{NULL}, uses
#'   built-in \code{scPharm::bulkdata}.
#' @param gdscdata GDSC pharmacogenomics data. If \code{NULL}, uses
#'   built-in \code{scPharm::gdscdata}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item For tissue samples: identifies tumor cells via CNV analysis (or uses provided annotations)
#'   \item Computes Multiple Correspondence Analysis (MCA) for dimensionality reduction
#'   \item Generates cell identity gene signatures

#'   \item Correlates gene expression with drug response (AUC) across cell lines
#'   \item Performs GSEA to score each cell's drug response profile
#'   \item Classifies cells based on NES thresholds
#' }
#'
#' @return A Seurat object with pharmacological annotations added to metadata:
#' \describe{
#'   \item{cell.label}{Cell type: "tumor" or "adjacent"}
#'   \item{scPharm_label_DRUGID_DRUGNAME}{Drug response: "sensitive", "resistant", or "other"}
#'   \item{scPharm_nes_DRUGID_DRUGNAME}{Normalized enrichment score for the drug}
#' }
#'
#' @import Seurat
#' @importFrom sparseMatrixStats rowVars
#' @importFrom irlba irlba
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom stats cor
#' @importFrom fgsea fgseaMultilevel
#' @import dplyr
#' @importFrom utils data
#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom stats cutree hclust ks.test median pf rnorm sd as.dist setNames var
#' @importFrom utils packageVersion
#'
#' @export
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- scPharmIdentify(seurat_obj, type = "tissue", cancer = "LUAD")
#'
#' # With known tumor cells
#' result <- scPharmIdentify(seurat_obj, type = "tissue", cancer = "LUAD",
#'                           tumor.cells = tumor_barcodes)
#'
#' # Specific drug analysis
#' result <- scPharmIdentify(seurat_obj, type = "cellline", cancer = "BRCA",
#'                           drug = "Erlotinib")
#' }
#'
#' @references
#' Tian P, Zheng J, et al. scPharm: identifying pharmacological subpopulations of
#' single cells for precision medicine in cancers. 2023.
#'
#' @useDynLib scPharm
scPharmIdentify <- function(object,
                            type,
                            cancer,
                            drug = NULL,
                            nmcs = 50,
                            nfeatures = 200,
                            cores = 1,
                            features = NULL,
                            slot = "data",
                            layer = NULL,
                            assay = "RNA",
                            threshold.s = -1.751302,
                            threshold.r = 1.518551,
                            tumor.cells = NULL,
                            normal.cells = NULL,
                            bulkdata = NULL,
                            gdscdata = NULL) {

  # Load reference data
  if (is.null(bulkdata)) bulkdata <- scPharm::bulkdata
  if (is.null(gdscdata)) gdscdata <- scPharm::gdscdata

  # Input validation
  if (!is(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!type %in% c("tissue", "cellline")) {
    stop("type must be 'tissue' or 'cellline'")
  }

  # Seurat version compatibility
  .getAssayData <- function(obj, slot_name, assay_name) {
    seurat_ver <- packageVersion("Seurat")
    if (seurat_ver >= "5.0.0") {
      layer_name <- if (!is.null(layer)) layer else slot_name
      tryCatch(
        GetAssayData(obj, layer = layer_name, assay = assay_name),
        error = function(e) GetAssayData(obj, slot = slot_name, assay = assay_name)
      )
    } else {
      GetAssayData(obj, slot = slot_name, assay = assay_name)
    }
  }

  init_assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  mat <- .getAssayData(object, slot, assay)

  # ========== Step 1: Identify tumor cells ==========
  if (type == "tissue") {
    if (!is.null(tumor.cells)) {
      message("Using provided tumor cell annotations")
      valid_tumor <- intersect(tumor.cells, colnames(object))
      if (length(valid_tumor) == 0) {
        stop("None of the provided tumor.cells found in object")
      }
      if (length(valid_tumor) < length(tumor.cells)) {
        warning(sprintf("%d tumor cells not found", length(tumor.cells) - length(valid_tumor)))
      }
      object@meta.data[, "cell.label"] <- "adjacent"
      object@meta.data[valid_tumor, "cell.label"] <- "tumor"
    } else {
      message("Identifying tumor cells via CNV analysis...")
      data("copykat_DNA.hg20", package = "scPharm", envir = environment())
      data("copykat_cyclegenes", package = "scPharm", envir = environment())
      data("copykat_full.anno.hg20", package = "scPharm", envir = environment())

      exp_rawdata <- as.matrix(.getAssayData(object, "counts", assay))
      copykat_result <- tryCatch({
        .run_copykat_internal(
          rawmat = exp_rawdata, id.type = "S", cell.line = FALSE,
          ngene.chr = 3, min.gene.per.cell = 200, LOW.DR = 0.05, UP.DR = 0.1,
          win.size = 25, norm.cell.names = normal.cells, KS.cut = 0.1,
          distance = "euclidean", n.cores = cores,
          full.anno = copykat_full.anno.hg20,
          DNA.hg20 = copykat_DNA.hg20,
          cyclegenes = copykat_cyclegenes
        )
      }, error = function(e) {
        warning(sprintf("CNV analysis failed: %s. Treating all cells as tumor.", e$message))
        data.frame(cell.names = colnames(object), copykat.pred = "aneuploid",
                   stringsAsFactors = FALSE)
      })

      object@meta.data[, "cell.label"] <- "adjacent"
      tumor_detected <- copykat_result$cell.names[copykat_result$copykat.pred == "aneuploid"]
      # Only update cells that exist in the object
      tumor_detected <- intersect(tumor_detected, rownames(object@meta.data))
      if (length(tumor_detected) > 0) {
        object@meta.data[tumor_detected, "cell.label"] <- "tumor"
      }

      n_tumor <- sum(object@meta.data$cell.label == "tumor")
      n_adj <- sum(object@meta.data$cell.label == "adjacent")
      message(sprintf("  Detected %d tumor cells, %d adjacent cells", n_tumor, n_adj))
    }
  } else {
    object@meta.data[, "cell.label"] <- "tumor"
  }

  # ========== Step 2: MCA dimensionality reduction ==========
  if (!is.null(features)) {
    mat <- mat[features, , drop = FALSE]
  }
  row_var <- rowVars(mat)
  mat <- mat[row_var != 0, , drop = FALSE]
  mat <- mat[nchar(rownames(mat)) > 0, , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

  cells_n <- colnames(mat)
  features_n <- rownames(mat)

  message("Computing MCA...")
  mca_prep <- SparseMCAStep1(mat)
  svd_res <- irlba(A = mca_prep$Z, nv = nmcs + 1, nu = 1)[seq(3)]
  mca_res <- MCAStep2(Z = mca_prep$Z, V = svd_res$v[, -1], Dc = mca_prep$Dc)

  component <- paste0("mca_", seq_len(ncol(mca_res$cellsCoordinates)))
  colnames(mca_res$cellsCoordinates) <- component
  colnames(mca_res$featuresCoordinates) <- component
  rownames(mca_res$cellsCoordinates) <- cells_n
  rownames(mca_res$featuresCoordinates) <- features_n

  dim_reduc <- CreateDimReducObject(
    embeddings = mca_res$cellsCoordinates,
    loadings = mca_res$featuresCoordinates,
    key = "mca_", assay = assay
  )
  object@reductions[["mca"]] <- dim_reduc
  object@reductions[["mca"]]@stdev <- sqrt(svd_res$d[-1])
  object@reductions[["mca"]]@misc[["mca.flag"]] <- TRUE
  DefaultAssay(object) <- init_assay

  # ========== Step 3: Cell identity signatures ==========
  message("Computing cell identity signatures...")
  cell_geneset <- .cellGeneSet(
    object, reduction = "mca", dims = seq_len(nmcs), n.features = nfeatures
  )

  # ========== Step 4: Filter by cancer type ==========
  GDSC <- .filterCancerType(gdscdata, cancer)
  drug_info <- GDSC[, c("DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME")]
  drug_info <- drug_info[!duplicated(drug_info$DRUG_ID), ]
  drug_info <- dplyr::filter(drug_info, !dplyr::if_all(dplyr::everything(), .fns = is.na))

  if (!is.null(drug)) {
    drug_info <- drug_info[drug_info$DRUG_NAME == drug, ]
    if (nrow(drug_info) == 0) {
      stop(sprintf("Drug '%s' not found for specified cancer type", drug))
    }
  }

  # ========== Step 5: GSEA for each drug ==========
  message("Running pharmacogenomics analysis...")
  meta_data <- object@meta.data
  n_drugs <- nrow(drug_info)

  for (i in seq_len(n_drugs)) {
    drug_id <- drug_info$DRUG_ID[i]
    drug_name <- drug_info$DRUG_NAME[i]
    message(sprintf("  [%d/%d] %s", i, n_drugs, drug_name))

    # Get drug response data
    drug_resp <- GDSC[GDSC$DRUG_ID == drug_id, c("CELL_LINE_NAME", "AUC")]
    drug_resp <- drug_resp[drug_resp$CELL_LINE_NAME %in% colnames(bulkdata), ]
    # Remove duplicates - keep first occurrence (usually most recent)
    drug_resp <- drug_resp[!duplicated(drug_resp$CELL_LINE_NAME), ]
    drug_resp <- drug_resp[order(drug_resp$AUC, decreasing = TRUE), ]

    if (nrow(drug_resp) < 3) {
      message(sprintf("    Skipped: insufficient data (n=%d)", nrow(drug_resp)))
      next
    }

    # Compute drug response gene signature
    sub_bulk <- bulkdata[, drug_resp$CELL_LINE_NAME, drop = FALSE]
    sub_bulk <- sub_bulk[rowSums(sub_bulk) > 0, , drop = FALSE]

    gene_cor <- apply(sub_bulk, 1, function(x) cor(x, drug_resp$AUC, method = "pearson"))
    gene_list <- sort(gene_cor, decreasing = TRUE)

    # GSEA
    gsea_res <- tryCatch({
      suppressWarnings(fgseaMultilevel(
        pathways = cell_geneset, stats = gene_list,
        eps = 0, minSize = 10, maxSize = 500,
        nproc = cores, nPermSimple = 10000
      ))
    }, error = function(e) NULL)

    if (is.null(gsea_res)) next
    gsea_res[is.na(gsea_res)] <- 0

    # Assign labels
    labels <- data.frame(
      cell = rownames(meta_data), nes = 0, label = "other",
      stringsAsFactors = FALSE
    )
    rownames(labels) <- labels$cell

    if (nrow(gsea_res) > 0) {
      # Only update cells that exist in both labels and gsea_res
      valid_pathways <- intersect(gsea_res$pathway, rownames(labels))
      if (length(valid_pathways) > 0) {
        gsea_subset <- gsea_res[gsea_res$pathway %in% valid_pathways, ]
        labels[gsea_subset$pathway, "nes"] <- gsea_subset$NES
      }
      labels$label[labels$nes > threshold.r] <- "resistant"
      labels$label[labels$nes < threshold.s] <- "sensitive"
    }

    col_label <- paste("scPharm_label", drug_id, drug_name, sep = "_")
    col_nes <- paste("scPharm_nes", drug_id, drug_name, sep = "_")
    meta_data[[col_label]] <- labels$label
    meta_data[[col_nes]] <- labels$nes
  }

  object@meta.data <- meta_data
  return(object)
}

#' Filter GDSC data by cancer type
#' @param GDSC GDSC data frame
#' @param cancer Cancer type(s)
#' @return Filtered data frame
#' @keywords internal
#' @noRd
.filterCancerType <- function(GDSC, cancer) {
  if (!is.character(cancer)) {
    stop("cancer must be a character or character vector")
  }
  if (length(cancer) == 1 && cancer == "pan") {
    return(GDSC)
  }
  valid_types <- unique(GDSC$TCGA_DESC)
  invalid <- setdiff(cancer, valid_types)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid cancer type(s): %s\nValid: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_types, collapse = ", ")))
  }
  GDSC[GDSC$TCGA_DESC %in% cancer, ]
}
