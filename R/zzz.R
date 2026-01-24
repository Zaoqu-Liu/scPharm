# Package initialization and Rcpp registration

#' @useDynLib scPharm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Suppress R CMD check notes for Rcpp functions
utils::globalVariables(c("_scPharm_SparseMCAStep1", "_scPharm_MCAStep2"))
