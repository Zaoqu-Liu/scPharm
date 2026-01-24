#' GDSC2 Pharmacogenomics Data
#'
#' Drug sensitivity data (IC50 and AUC) for cancer cell lines from the
#' Genomics of Drug Sensitivity in Cancer (GDSC) project.
#'
#' @format A data frame with 196,344 observations and 19 variables including:
#' \describe{
#'   \item{DRUG_ID}{GDSC drug identifier}
#'   \item{DRUG_NAME}{Drug name}
#'   \item{CELL_LINE_NAME}{Cell line identifier}
#'   \item{COSMIC_ID}{COSMIC cell line ID}
#'   \item{SANGER_MODEL_ID}{Sanger model ID}
#'   \item{TCGA_DESC}{TCGA cancer type classification}
#'   \item{DATASET}{Dataset source}
#'   \item{COMPANY_ID}{Company identifier}
#'   \item{NLME_RESULT_ID}{NLME result ID}
#'   \item{NLME_CURVE_ID}{NLME curve ID}
#'   \item{LN_IC50}{Natural log of IC50}
#'   \item{AUC}{Area under the dose-response curve}
#'   \item{RMSE}{Root mean square error}
#'   \item{Z_SCORE}{Z-score}
#'   \item{MAX_CONC}{Maximum concentration tested}
#'   \item{MIN_CONC}{Minimum concentration tested}
#'   \item{PUTATIVE_TARGET}{Known drug target}
#'   \item{PATHWAY_NAME}{Target pathway}
#'   \item{WEBRELEASE}{Web release version}
#' }
#'
#' @source GDSC \url{https://www.cancerrxgene.org/downloads/drug_data}
#'
#' @references
#' Yang W, et al. (2013). Genomics of Drug Sensitivity in Cancer (GDSC):
#' a resource for therapeutic biomarker discovery in cancer cells.
#' Nucleic Acids Research 41, D955-D961.
#'
#' @name gdscdata
#' @keywords datasets
"gdscdata"
