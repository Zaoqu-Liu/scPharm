#' Genomic Bin Coordinates for hg38
#'
#' Contains 220KB genomic bin coordinates for human genome hg38/GRCh38,
#' used for CNV analysis.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{chrom}{Chromosome number}
#'   \item{chrompos}{Chromosome position}
#'   \item{abspos}{Absolute genome position}
#' }
#'
#' @source Gao et al., Nature Biotechnology 2021
#' @name copykat_DNA.hg20
#' @keywords internal
"copykat_DNA.hg20"

#' Cell Cycle Genes for CNV Analysis
#'
#' Cell cycle-related genes excluded from CNV analysis to reduce noise.
#'
#' @format A list containing cell cycle gene symbols
#'
#' @source Gao et al., Nature Biotechnology 2021
#' @name copykat_cyclegenes
#' @keywords internal
"copykat_cyclegenes"

#' Human Gene Annotations for hg38
#'
#' Gene annotations with genomic coordinates for human genes (hg38/GRCh38).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{hgnc_symbol}{HGNC gene symbol}
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome_name}{Chromosome number}
#'   \item{start_position}{Gene start position}
#'   \item{end_position}{Gene end position}
#'   \item{band}{Chromosome band}
#'   \item{abspos}{Absolute genome position}
#' }
#'
#' @source Gao et al., Nature Biotechnology 2021
#' @name copykat_full.anno.hg20
#' @keywords internal
"copykat_full.anno.hg20"

# Suppress R CMD check notes for data loading
utils::globalVariables(c("copykat_DNA.hg20", "copykat_cyclegenes", "copykat_full.anno.hg20"))
