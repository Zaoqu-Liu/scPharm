#' Identify Potential Drug Combinations
#'
#' Identify potential drug combinations based on two strategies:
#' (1) Compensation effects: drugs that target cells resistant to the primary drug
#' (2) Booster effects: drugs that enhance sensitivity through different pathways
#'
#' @param object A Seurat object after running \code{\link{scPharmIdentify}}.
#' @param score Output from \code{\link{scPharmDr}}.
#' @param drug Name of the primary drug. If \code{NULL}, uses \code{topN} drugs.
#' @param topN Number of top-ranked drugs to analyze. Default: 1.
#' @param drug_info Drug information table. If \code{NULL}, uses built-in
#'   \code{scPharm::drug_info}.
#'
#' @details
#' \strong{Compensation effects}: Identifies drugs where cells resistant to the
#' primary drug show sensitivity. This suggests the combination could overcome
#' resistance.
#'
#' \strong{Booster effects}: Identifies drugs targeting different pathways that
#' show high sensitivity in cells already sensitive to the primary drug. This
#' suggests synergistic enhancement.
#'
#' @return A named list where each element corresponds to a primary drug and
#'   contains a data frame with:
#' \describe{
#'   \item{DRUG_FIRST}{Primary drug name}
#'   \item{DRUG_ID}{Combination drug ID}
#'   \item{DRUG_NAME}{Combination drug name}
#'   \item{Effect}{Combination effect score}
#'   \item{Strategy}{"compensation effects" or "booster effects"}
#' }
#'
#' @importFrom mixtools normalmixEM
#'
#' @export
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' dr_scores <- scPharmDr(result)
#' combos <- scPharmCombo(result, dr_scores, topN = 3)
#' }
#'
#' @seealso \code{\link{scPharmIdentify}}, \code{\link{scPharmDr}}
scPharmCombo <- function(object, score, drug = NULL, topN = 1, drug_info = NULL) {

  if (is.null(drug_info)) drug_info <- scPharm::drug_info

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  if (!is.data.frame(score) || !"DRUG_NAME" %in% colnames(score)) {
    stop("score must be output from scPharmDr")
  }

  # Select primary drugs
  if (is.null(drug)) {
    if (topN > nrow(score)) {
      warning(sprintf("topN (%d) > available drugs (%d). Using all.", topN, nrow(score)))
      topN <- nrow(score)
    }
    score <- score[seq_len(topN), , drop = FALSE]
  } else {
    score <- score[score$DRUG_NAME == drug, , drop = FALSE]
    if (nrow(score) == 0) {
      stop(sprintf("Drug '%s' not found in score data", drug))
    }
  }

  message(sprintf("Analyzing combinations for %d drug(s)...", nrow(score)))

  meta_data <- object@meta.data
  nes_cols <- grep("^scPharm_nes_", colnames(meta_data), value = TRUE)
  label_cols <- grep("^scPharm_label_", colnames(meta_data), value = TRUE)

  if (length(nes_cols) == 0 || length(label_cols) == 0) {
    stop("scPharm results not found. Run scPharmIdentify first.")
  }

  nes_data <- meta_data[, nes_cols, drop = FALSE]
  label_data <- meta_data[, label_cols, drop = FALSE]

  combo_result <- list()

  for (i in seq_len(nrow(score))) {
    drug_id <- score[i, "DRUG_ID"]
    drug_name <- score[i, "DRUG_NAME"]
    message(sprintf("  [%d/%d] %s", i, nrow(score), drug_name))

    col_label <- paste("scPharm_label", drug_id, drug_name, sep = "_")
    col_nes <- paste("scPharm_nes", drug_id, drug_name, sep = "_")

    if (!col_label %in% colnames(label_data) || !col_nes %in% colnames(nes_data)) {
      warning(sprintf("Data not found for %s, skipping", drug_name))
      next
    }

    # ===== Compensation effects =====
    nes_vals <- nes_data[[col_nes]]
    nes_vals <- nes_vals[!is.na(nes_vals)]

    combo_comp <- NULL
    if (length(nes_vals) > 50) {
      gmm <- tryCatch(
        suppressMessages(normalmixEM(nes_vals, k = 2, maxit = 500)),
        error = function(e) NULL
      )

      if (!is.null(gmm) && max(gmm$mu) > 0) {
        resistant <- label_data[[col_label]] == "resistant"
        if (sum(resistant, na.rm = TRUE) > 0) {
          sub_label <- label_data[resistant, , drop = FALSE]
          overlap <- vapply(sub_label, function(x) mean(x == "sensitive", na.rm = TRUE), numeric(1))
          overlap <- overlap[names(overlap) != col_label]
          overlap <- sort(overlap, decreasing = TRUE)

          n_top <- min(5, length(overlap))
          if (n_top > 0) {
            overlap <- overlap[seq_len(n_top)]
            combo_comp <- .formatComboResult(drug_name, names(overlap), overlap, "compensation effects")
          }
        }
      }
    }

    # ===== Booster effects =====
    target_path <- drug_info[drug_info$DRUG_ID == drug_id, "PATHWAY_NAME"]
    # Handle multiple matches - take first non-NA value
    target_path <- target_path[!is.na(target_path)]
    if (length(target_path) > 0) target_path <- target_path[1]

    combo_boost <- NULL
    if (length(target_path) > 0 && !is.na(target_path) &&
        !grepl("Other|Unclassified", target_path, ignore.case = TRUE)) {

      sensitive <- label_data[[col_label]] == "sensitive"
      if (sum(sensitive, na.rm = TRUE) > 0) {
        sub_label <- label_data[sensitive, , drop = FALSE]
        enhance <- vapply(sub_label, function(x) mean(x == "sensitive", na.rm = TRUE), numeric(1))

        # Remove same pathway drugs
        same_path <- drug_info[drug_info$PATHWAY_NAME == target_path, c("DRUG_ID", "DRUG_NAME")]
        to_remove <- paste("scPharm_label", same_path$DRUG_ID, same_path$DRUG_NAME, sep = "_")
        enhance <- enhance[!names(enhance) %in% to_remove]
        enhance <- sort(enhance, decreasing = TRUE)

        n_top <- min(5, length(enhance))
        if (n_top > 0) {
          enhance <- enhance[seq_len(n_top)]
          combo_boost <- .formatComboResult(drug_name, names(enhance), enhance, "booster effects")
        }
      }
    }

    # Combine results
    combo <- NULL
    if (!is.null(combo_boost) && !is.null(combo_comp)) {
      combo <- rbind(combo_boost, combo_comp)
    } else if (!is.null(combo_boost)) {
      combo <- combo_boost
    } else {
      combo <- combo_comp
    }

    if (!is.null(combo)) {
      combo <- combo[order(-combo$Effect), ]
    }

    combo_result[[paste(drug_id, drug_name, sep = "_")]] <- combo
  }

  return(combo_result)
}

#' Format combination results
#' @keywords internal
#' @noRd
.formatComboResult <- function(drug_first, col_names, effect_vals, strategy) {
  parse_drug <- function(col) {
    parts <- strsplit(col, "_")[[1]]
    id <- parts[3]
    name <- if (length(parts) == 4) parts[4] else paste(parts[4:length(parts)], collapse = "_")
    c(id, name)
  }

  info <- vapply(col_names, parse_drug, character(2))

  data.frame(
    DRUG_FIRST = drug_first,
    DRUG_ID = info[1, ],
    DRUG_NAME = info[2, ],
    Effect = effect_vals,
    Strategy = strategy,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
