#' Compute Drug Prioritization Score (Dr)
#'
#' Calculate drug prioritization scores based on the ratio of sensitive and
#' resistant tumor cell populations identified by \code{\link{scPharmIdentify}}.
#'
#' @param object A Seurat object after running \code{\link{scPharmIdentify}}.
#'
#' @details
#' The Dr score is calculated as: \deqn{Dr = S \times (1 - R)}
#' where S is the proportion of sensitive cells and R is the proportion of

#' resistant cells among tumor cells. Higher scores indicate better drug
#' candidates for the patient.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{DRUG_ID}{GDSC drug identifier}
#'   \item{DRUG_NAME}{Drug name}
#'   \item{SENSI_RATIO}{Proportion of sensitive tumor cells}
#'   \item{RESIS_RATIO}{Proportion of resistant tumor cells}
#'   \item{Dr}{Drug prioritization score}
#'   \item{Rank}{Drug ranking (1 = best)}
#' }
#'
#' @export
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' # After running scPharmIdentify
#' dr_scores <- scPharmDr(result)
#' head(dr_scores)
#' }
#'
#' @seealso \code{\link{scPharmIdentify}}, \code{\link{scPharmCombo}}
scPharmDr <- function(object) {
  message("Computing drug prioritization scores...")

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

  # Filter tumor cells
  tumor_data <- meta_data[meta_data$cell.label == "tumor", label_cols, drop = FALSE]
  if (nrow(tumor_data) == 0) {
    stop("No tumor cells found in data")
  }

  # Calculate ratios
  sensi_ratio <- vapply(tumor_data, function(x) mean(x == "sensitive", na.rm = TRUE), numeric(1))
  resis_ratio <- vapply(tumor_data, function(x) mean(x == "resistant", na.rm = TRUE), numeric(1))

  # Parse drug info from column names
  parse_drug <- function(col) {
    parts <- strsplit(col, "_")[[1]]
    id <- parts[3]
    name <- if (length(parts) == 4) parts[4] else paste(parts[4:length(parts)], collapse = "_")
    c(id, name)
  }

  drug_info <- vapply(names(sensi_ratio), parse_drug, character(2))

  score <- data.frame(
    DRUG_ID = drug_info[1, ],
    DRUG_NAME = drug_info[2, ],
    SENSI_RATIO = sensi_ratio,
    RESIS_RATIO = resis_ratio,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  score$Dr <- score$SENSI_RATIO * (1 - score$RESIS_RATIO)
  score <- score[order(-score$Dr), ]
  score$Rank <- seq_len(nrow(score))

  message(sprintf("Ranked %d drugs. Top: %s (Dr=%.4f)", nrow(score), score$DRUG_NAME[1], score$Dr[1]))
  return(score)
}
