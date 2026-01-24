#' Generate Null Distribution and Thresholds
#'
#' Generate a null distribution from healthy tissue cells and calculate
#' thresholds for classifying sensitive and resistant cells.
#'
#' @param object A Seurat object containing cells from healthy/normal tissue.
#' @param cancer TCGA cancer type(s) for context. A character string or vector.
#'   Use \code{"pan"} for pan-cancer analysis.
#' @param nmcs Number of MCA components. Default: 50.
#' @param nfeatures Number of genes for cell identity signature. Default: 200.
#' @param cores Number of CPU cores. Default: 1.
#' @param features Character vector of gene names to use. If \code{NULL}, uses all.
#' @param slot Slot for Seurat V4. Default: \code{"data"}.
#' @param layer Layer for Seurat V5. If \code{NULL}, uses \code{slot} value.
#' @param assay Assay to use. Default: \code{"RNA"}.
#' @param bulkdata Bulk RNA-seq data. If \code{NULL}, uses built-in data.
#' @param gdscdata GDSC data. If \code{NULL}, uses built-in data.
#'
#' @details
#' This function computes NES distributions from normal cells and uses a
#' two-component Gaussian mixture model to determine thresholds. The thresholds
#' are calculated as mean +/- 1 standard deviation of each component.
#'
#' @return A list containing:
#' \describe{
#'   \item{NullDist}{Numeric vector of NES values from normal cells}
#'   \item{threshold_s}{Threshold for sensitive cells (NES < threshold_s)}
#'   \item{threshold_r}{Threshold for resistant cells (NES > threshold_r)}
#' }
#'
#' @import Seurat
#' @importFrom sparseMatrixStats rowVars
#' @importFrom irlba irlba
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom CelliD GetCellGeneSet
#' @importFrom stats cor
#' @importFrom fgsea fgseaMultilevel
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom mixtools normalmixEM
#' @importFrom methods is
#' @importFrom utils packageVersion
#'
#' @export
#' @encoding UTF-8
#'
#' @useDynLib scPharm
#'
#' @examples
#' \dontrun{
#' # Using healthy tissue cells
#' thresholds <- scPharmGenNullDist(healthy_seurat, cancer = "BRCA")
#' print(thresholds$threshold_s)
#' print(thresholds$threshold_r)
#' }
#'
#' @seealso \code{\link{scPharmIdentify}}
scPharmGenNullDist <- function(object,
                               cancer,
                               nmcs = 50,
                               nfeatures = 200,
                               cores = 1,
                               features = NULL,
                               slot = "data",
                               layer = NULL,
                               assay = "RNA",
                               bulkdata = NULL,
                               gdscdata = NULL) {

  if (is.null(bulkdata)) bulkdata <- scPharm::bulkdata
  if (is.null(gdscdata)) gdscdata <- scPharm::gdscdata

  if (!is(object, "Seurat")) {
    stop("Input must be a Seurat object")
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

  # ========== MCA ==========
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

  # ========== Cell signatures ==========
  message("Computing cell identity signatures...")
  cell_geneset <- suppressMessages(GetCellGeneSet(
    object, reduction = "mca", dims = seq_len(nmcs), n.features = nfeatures
  ))

  # ========== Filter by cancer type ==========
  GDSC <- gdscdata
  if (!is.character(cancer)) {
    stop("cancer must be a character or character vector")
  }

  if (!(length(cancer) == 1 && cancer == "pan")) {
    valid_types <- unique(GDSC$TCGA_DESC)
    invalid <- setdiff(cancer, valid_types)
    if (length(invalid) > 0) {
      stop(sprintf("Invalid cancer type(s): %s", paste(invalid, collapse = ", ")))
    }
    GDSC <- GDSC[GDSC$TCGA_DESC %in% cancer, ]
  }

  drug_info <- GDSC[, c("DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME")]
  drug_info <- drug_info[!duplicated(drug_info$DRUG_ID), ]
  drug_info <- dplyr::filter(drug_info, !dplyr::if_all(dplyr::everything(), .fns = is.na))

  # ========== GSEA for all drugs ==========
  message("Running pharmacogenomics analysis...")
  meta_data <- object@meta.data
  n_drugs <- nrow(drug_info)

  for (i in seq_len(n_drugs)) {
    drug_id <- drug_info$DRUG_ID[i]
    drug_name <- drug_info$DRUG_NAME[i]
    message(sprintf("  [%d/%d] %s", i, n_drugs, drug_name))

    drug_resp <- GDSC[GDSC$DRUG_ID == drug_id, c("CELL_LINE_NAME", "AUC")]
    drug_resp <- drug_resp[drug_resp$CELL_LINE_NAME %in% colnames(bulkdata), ]
    drug_resp <- drug_resp[!duplicated(drug_resp$CELL_LINE_NAME), ]
    drug_resp <- drug_resp[order(drug_resp$AUC, decreasing = TRUE), ]

    if (nrow(drug_resp) < 3) {
      message(sprintf("    Skipped: insufficient data (n=%d)", nrow(drug_resp)))
      next
    }

    sub_bulk <- bulkdata[, drug_resp$CELL_LINE_NAME, drop = FALSE]
    sub_bulk <- sub_bulk[rowSums(sub_bulk) > 0, , drop = FALSE]

    gene_cor <- apply(sub_bulk, 1, function(x) cor(x, drug_resp$AUC, method = "pearson"))
    gene_list <- sort(gene_cor, decreasing = TRUE)

    gsea_res <- tryCatch({
      suppressWarnings(fgseaMultilevel(
        pathways = cell_geneset, stats = gene_list,
        eps = 0, minSize = 10, maxSize = 500,
        nproc = cores, nPermSimple = 10000
      ))
    }, error = function(e) NULL)

    if (is.null(gsea_res)) next
    gsea_res[is.na(gsea_res)] <- 0

    labels <- data.frame(cell = rownames(meta_data), nes = 0, stringsAsFactors = FALSE)
    rownames(labels) <- labels$cell

    if (nrow(gsea_res) > 0) {
      valid_pathways <- intersect(gsea_res$pathway, rownames(labels))
      if (length(valid_pathways) > 0) {
        gsea_subset <- gsea_res[gsea_res$pathway %in% valid_pathways, ]
        labels[gsea_subset$pathway, "nes"] <- gsea_subset$NES
      }
    }

    col_name <- paste("scPharm_nes", drug_id, drug_name, sep = "_")
    meta_data[[col_name]] <- labels$nes
  }

  # ========== Compute thresholds ==========
  nes_cols <- grep("^scPharm_nes", colnames(meta_data), value = TRUE)
  if (length(nes_cols) == 0) {
    stop("No NES scores computed. Check input data.")
  }

  nes_data <- meta_data[, nes_cols, drop = FALSE]
  nes_long <- gather(nes_data, key = "Drug", value = "NES")
  nes_values <- nes_long$NES[!is.na(nes_long$NES)]

  if (length(nes_values) < 100) {
    warning("Limited data points for threshold estimation")
  }

  # Fit GMM
  gmm <- suppressMessages(normalmixEM(nes_values, k = 2, maxit = 1000))

  # Determine thresholds
  if (gmm$mu[1] < gmm$mu[2]) {
    threshold_s <- gmm$mu[1] - gmm$sigma[1]
    threshold_r <- gmm$mu[2] + gmm$sigma[2]
  } else {
    threshold_s <- gmm$mu[2] - gmm$sigma[2]
    threshold_r <- gmm$mu[1] + gmm$sigma[1]
  }

  message(sprintf("Thresholds: sensitive < %.4f, resistant > %.4f", threshold_s, threshold_r))

  return(list(
    NullDist = nes_values,
    threshold_r = threshold_r,
    threshold_s = threshold_s
  ))
}
