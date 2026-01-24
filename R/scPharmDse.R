#' Predict Drug Side Effects (Dse)
#'
#' Calculate drug side effect scores based on the sensitivity of adjacent
#' (normal) cells to drugs. Higher scores indicate greater potential for
#' off-target toxicity.
#'
#' @param object A Seurat object after running \code{\link{scPharmIdentify}}
#'   with \code{type="tissue"}.
#'
#' @details
#' The Dse score represents the proportion of adjacent (non-tumor) cells that
#' are classified as sensitive to each drug. Drugs with high Dse scores may
#' cause more side effects by affecting normal cells.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{DRUG_ID}{GDSC drug identifier}
#'   \item{DRUG_NAME}{Drug name}
#'   \item{Dse}{Side effect score (0-1, higher = more side effects)}
#' }
#'
#' @export
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' # After running scPharmIdentify with type="tissue"
#' dse_scores <- scPharmDse(result)
#' head(dse_scores)
#' }
#'
#' @seealso \code{\link{scPharmIdentify}}, \code{\link{scPharmDr}}
scPharmDse <- function(object) {
  message("Computing drug side effect scores...")

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  meta_data <- object@meta.data

  if (!"cell.label" %in% colnames(meta_data)) {
    stop("cell.label not found. Run scPharmIdentify first.")
  }

  label_cols <- grep("^scPharm_label_", colnames(meta_data), value = TRUE)
  if (length(label_cols) == 0) {
    stop("No scPharm results found. Run scPharmIdentify first.")
  }

  # Parse drug info
  parse_drug <- function(col) {
    parts <- strsplit(col, "_")[[1]]
    id <- parts[3]
    name <- if (length(parts) == 4) parts[4] else paste(parts[4:length(parts)], collapse = "_")
    c(id, name)
  }

  # Filter adjacent cells
  adj_data <- meta_data[meta_data$cell.label == "adjacent", label_cols, drop = FALSE]

  if (nrow(adj_data) == 0) {
    warning("No adjacent cells found. Returning NA for all drugs.")
    drug_info <- vapply(label_cols, parse_drug, character(2))
    return(data.frame(
      DRUG_ID = drug_info[1, ],
      DRUG_NAME = drug_info[2, ],
      Dse = NA_real_,
      row.names = NULL,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate sensitivity ratio in adjacent cells
  sensi_ratio <- vapply(adj_data, function(x) mean(x == "sensitive", na.rm = TRUE), numeric(1))
  drug_info <- vapply(names(sensi_ratio), parse_drug, character(2))

  score <- data.frame(
    DRUG_ID = drug_info[1, ],
    DRUG_NAME = drug_info[2, ],
    Dse = sensi_ratio,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  score <- score[order(-score$Dse), ]
  message(sprintf("Calculated Dse for %d drugs", nrow(score)))
  return(score)
}
