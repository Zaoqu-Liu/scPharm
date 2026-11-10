# Cell-specific gene signatures from an MCA reduction.
#
# In MCA, a gene lies closer to a cell the more specific it is to that cell.
# The signature of a cell is therefore its `n.features` nearest genes by
# Euclidean distance in the MCA space. This replaces the former dependency on
# CelliD::GetCellGeneSet, because CelliD was removed from Bioconductor 3.23
# and can no longer be installed from a standard repository.
#
# Distances are computed in blocks of cells so that memory stays bounded for
# large datasets instead of materialising the full gene x cell matrix.
.cellGeneSet <- function(object, reduction = "mca", dims, n.features = 200,
                         max.block = 5e7) {
  cells_emb <- Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  genes_emb <- Loadings(object, reduction = reduction)[, dims, drop = FALSE]
  n_genes <- nrow(genes_emb)
  n_keep <- min(as.integer(n.features), n_genes)
  gene_names <- rownames(genes_emb)

  genes_sq <- rowSums(genes_emb^2)
  cells_sq <- rowSums(cells_emb^2)
  block <- max(1L, as.integer(floor(max.block / n_genes)))
  starts <- seq(1L, nrow(cells_emb), by = block)

  geneset <- vector("list", nrow(cells_emb))
  names(geneset) <- rownames(cells_emb)
  for (s in starts) {
    idx <- s:min(s + block - 1L, nrow(cells_emb))
    # squared distance = |g|^2 - 2 g.c + |c|^2; sqrt is monotone, so the
    # ranking equals the Euclidean ranking.
    d2 <- genes_sq - 2 * (genes_emb %*% t(cells_emb[idx, , drop = FALSE]))
    d2 <- sweep(d2, 2, cells_sq[idx], "+")
    for (j in seq_along(idx)) {
      geneset[[idx[j]]] <- gene_names[order(d2[, j], method = "radix")[seq_len(n_keep)]]
    }
  }
  geneset
}
