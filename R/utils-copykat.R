# ============================================================================
# Internal functions for CNV-based tumor cell identification
# Adapted from copykat (Gao et al., Nature Biotechnology 2021)
# ============================================================================

# Windows-compatible parallel distance computation
# parallelDist::parDist can crash R on Windows due to OpenMP threading issues
.safe_parDist <- function(x, threads = 1, method = "euclidean") {
  is_windows <- .Platform$OS.type == "windows"

  if (is_windows || threads <= 1) {
    # Use base R dist() on Windows - slower but stable
    return(dist(x, method = method))
  }

  # On Unix/Mac with multiple threads, use parallelDist
  if (requireNamespace("parallelDist", quietly = TRUE)) {
    tryCatch(
      parallelDist::parDist(x, threads = threads, method = method),
      error = function(e) {
        warning("parallelDist failed, falling back to base dist()")
        dist(x, method = method)
      }
    )
  } else {
    dist(x, method = method)
  }
}

# Windows-compatible parallel lapply
# mclapply uses fork() which is not available on Windows
.safe_mclapply <- function(X, FUN, mc.cores = 1) {
  is_windows <- .Platform$OS.type == "windows"

  if (is_windows || mc.cores <= 1) {
    # Use sequential lapply on Windows
    return(lapply(X, FUN))
  }

  # On Unix/Mac, use mclapply
  parallel::mclapply(X, FUN, mc.cores = mc.cores)
}

#' Annotate genes with genomic coordinates
#' @param mat Expression matrix
#' @param ID.type Gene ID type
#' @param full.anno Annotation data
#' @keywords internal
#' @noRd
.copykat_annotateGenes_hg20 <- function(mat, ID.type = "S", full.anno) {
  if (substring(ID.type, 1, 1) %in% c("E", "e")) {
    shared <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[rownames(mat) %in% shared, , drop = FALSE]
    anno <- full.anno[full.anno$ensembl_gene_id %in% shared, ]
    anno <- anno[!duplicated(anno$hgnc_symbol), ]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))), ]
  } else {
    shared <- intersect(rownames(mat), full.anno$hgnc_symbol)
    mat <- mat[rownames(mat) %in% shared, , drop = FALSE]
    anno <- full.anno[full.anno$hgnc_symbol %in% shared, ]
    anno <- anno[!duplicated(anno$hgnc_symbol), ]
    anno <- anno[order(match(anno$hgnc_symbol, rownames(mat))), ]
  }
  cbind(anno, mat)
}

#' Find diploid baseline using hierarchical clustering
#' @param norm.mat.smooth Smoothed expression matrix
#' @param min.cells Minimum cells per cluster
#' @param n.cores Number of cores
#' @keywords internal
#' @noRd
.copykat_baseline_norm_cl <- function(norm.mat.smooth, min.cells = 5, n.cores = 1) {
  d <- .safe_parDist(t(norm.mat.smooth), threads = n.cores, method = "euclidean")
  km <- 6
  fit <- hclust(d, method = "ward.D2")
  ct <- cutree(fit, k = km)

  while (!all(table(ct) > min.cells)) {
    km <- km - 1
    ct <- cutree(fit, k = km)
    if (km == 2) break
  }

  clusters <- sort(unique(ct))
  SDM <- numeric(length(clusters))

  for (idx in seq_along(clusters)) {
    cols <- which(ct == clusters[idx])
    if (length(cols) == 1) {
      data_c <- norm.mat.smooth[, cols]
    } else {
      data_c <- matrixStats::rowMedians(norm.mat.smooth[, cols, drop = FALSE])
    }
    sx <- max(0.05, 0.5 * sd(data_c))
    gmm <- mixtools::normalmixEM(data_c, lambda = rep(1/3, 3), mu = c(-0.2, 0, 0.2),
                                  sigma = sx, arbvar = FALSE, ECM = FALSE, maxit = 2000)
    SDM[idx] <- gmm$sigma[1]
  }

  wn <- mean(cluster::silhouette(cutree(fit, k = 2), d)[, "sil_width"])
  pdt <- pf(max(SDM)^2 / min(SDM)^2, nrow(norm.mat.smooth), nrow(norm.mat.smooth), lower.tail = FALSE)

  WNS <- if (wn <= 0.15 || !all(table(ct) > min.cells) || pdt > 0.05) "unclassified.prediction" else ""

  best_cluster <- clusters[which.min(SDM)]
  cols <- which(ct == best_cluster)
  basel <- if (length(cols) == 1) norm.mat.smooth[, cols] else matrixStats::rowMedians(norm.mat.smooth[, cols, drop = FALSE])

  list(basel = basel, WNS = WNS, preN = colnames(norm.mat.smooth)[cols], cl = ct)
}

#' Find diploid baseline using GMM
#' @param CNA.mat Expression matrix
#' @param max.normal Maximum normal cells to find
#' @param mu.cut Mean cutoff
#' @param Nfraq.cut Fraction cutoff
#' @param RE.before Previous result for fallback
#' @param n.cores Number of cores
#' @keywords internal
#' @noRd
.copykat_baseline_GMM <- function(CNA.mat, max.normal = 5, mu.cut = 0.05,
                                  Nfraq.cut = 0.99, RE.before = NULL, n.cores = 1) {
  predictions <- character(0)
  for (m in seq_len(ncol(CNA.mat))) {
    sam <- CNA.mat[, m]
    sg <- max(0.05, 0.5 * sd(sam))
    gmm <- mixtools::normalmixEM(sam, lambda = rep(1/3, 3), mu = c(-0.2, 0, 0.2),
                                  sigma = sg, arbvar = FALSE, ECM = FALSE, maxit = 500)
    s <- sum(abs(gmm$mu) <= mu.cut)
    if (s >= 1 && sum(gmm$lambda[abs(gmm$mu) <= mu.cut]) > Nfraq.cut) {
      pred <- "diploid"
    } else {
      pred <- "aneuploid"
    }
    predictions <- c(predictions, pred)
    if (sum(predictions == "diploid") >= max.normal) break
  }

  names(predictions) <- colnames(CNA.mat)[seq_along(predictions)]
  preN <- names(predictions)[predictions == "diploid"]

  d <- .safe_parDist(t(CNA.mat), threads = n.cores, method = "euclidean")
  fit <- hclust(d, method = "ward.D2")
  ct <- cutree(fit, k = 6)

  if (length(preN) > 2) {
    basel <- apply(CNA.mat[, colnames(CNA.mat) %in% preN, drop = FALSE], 1, mean)
    list(basel = basel, WNS = "", preN = preN, cl = ct)
  } else {
    RE.before
  }
}

#' Create synthetic baseline for cell line mode
#' @param norm.mat Normalized matrix
#' @param min.cells Minimum cells per cluster
#' @param n.cores Number of cores
#' @keywords internal
#' @noRd
.copykat_baseline_synthetic <- function(norm.mat, min.cells = 10, n.cores = 1) {
  d <- .safe_parDist(t(norm.mat), threads = n.cores, method = "euclidean")
  fit <- hclust(d, method = "ward.D2")
  ct <- cutree(fit, k = 6)

  while (!all(table(ct) > min.cells)) {
    ct <- cutree(fit, k = max(2, length(unique(ct)) - 1))
    if (length(unique(ct)) == 2) break
  }

  expr_list <- lapply(sort(unique(ct)), function(i) {
    cols <- which(ct == i)
    data_c <- norm.mat[, cols, drop = FALSE]
    sd_vals <- if (ncol(data_c) == 1) rep(0, nrow(data_c)) else matrixStats::rowSds(data_c)
    set.seed(123)
    syn <- sapply(sd_vals, function(x) rnorm(1, mean = 0, sd = x))
    t(data_c - syn)
  })

  list(expr.relat = data.frame(t(do.call(rbind, expr_list))), cl = ct)
}

#' MCMC-based CNV segmentation
#' @param clu Cell clusters
#' @param fttmat Transformed matrix
#' @param bins Window size
#' @param cut.cor KS cutoff
#' @param n.cores Number of cores
#' @keywords internal
#' @noRd
.copykat_CNA_MCMC <- function(clu, fttmat, bins, cut.cor, n.cores) {
  clusters <- sort(unique(clu))
  n <- nrow(fttmat)

  # Compute cluster medians
  CON <- matrix(0, nrow = n, ncol = length(clusters))
  for (idx in seq_along(clusters)) {
    cols <- which(clu == clusters[idx])
    CON[, idx] <- if (length(cols) == 1) fttmat[, cols] else matrixStats::rowMedians(fttmat[, cols, drop = FALSE])
  }

  # Find breakpoints
  norm_sm <- exp(CON)
  n_bins <- max(1, as.integer(n / bins - 1))
  breks <- c(seq(1, n_bins * bins, bins), n)
  BR <- NULL

  for (c in seq_len(ncol(norm_sm))) {
    bre <- NULL
    for (i in seq_len(length(breks) - 2)) {
      idx1 <- breks[i]:breks[i + 1]
      idx2 <- (breks[i + 1] + 1):breks[i + 2]
      post1 <- MCMCpack::MCpoissongamma(norm_sm[idx1, c], max(mean(norm_sm[idx1, c]), 0.001), 1, mc = 1000)
      post2 <- MCMCpack::MCpoissongamma(norm_sm[idx2, c], max(mean(norm_sm[idx2, c]), 0.001), 1, mc = 1000)
      if (ks.test(post1, post2)$statistic[[1]] > cut.cor) {
        bre <- c(bre, breks[i + 1])
      }
    }
    BR <- sort(unique(c(BR, 1, bre, n)))
  }

  # Segment all cells
  norm_sm <- exp(fttmat)
  seg_starts <- BR[-length(BR)]
  seg_ends <- BR[-1]

  seg_fn <- function(z) {
    x <- numeric(n)
    for (i in seq_along(seg_starts)) {
      idx <- seg_starts[i]:seg_ends[i]
      post <- MCMCpack::MCpoissongamma(norm_sm[idx, z], max(mean(norm_sm[idx, z]), 0.001), 1, mc = 1000)
      x[idx] <- mean(post)
    }
    log(x)
  }

  seg_res <- .safe_mclapply(seq_len(ncol(norm_sm)), seg_fn, mc.cores = n.cores)
  logCNA <- matrix(unlist(seg_res), ncol = ncol(norm_sm), byrow = FALSE)

  list(logCNA = logCNA, breaks = BR)
}

#' Convert to genomic bins
#' @param DNA.mat Bin coordinates
#' @param RNA.mat Expression matrix
#' @param full.anno Annotations
#' @param n.cores Number of cores
#' @keywords internal
#' @noRd
.copykat_convert_bins_hg20 <- function(DNA.mat, RNA.mat, full.anno, n.cores = 1) {
  DNA <- DNA.mat[DNA.mat$chrom != 24, , drop = FALSE]
  anno <- full.anno[full.anno$hgnc_symbol %in% RNA.mat$hgnc_symbol &
                      full.anno$chromosome_name != 24,
                    c("hgnc_symbol", "chromosome_name", "start_position", "end_position")]
  anno$center <- (anno$start_position + anno$end_position) / 2

  # Map genes to bins
  bins_by_chr <- split(seq_len(nrow(DNA)), DNA$chrom)
  anno$bin_row <- NA_integer_
  for (chr in names(bins_by_chr)) {
    idx <- anno$chromosome_name == as.numeric(chr)
    if (!any(idx)) next
    chr_bins <- bins_by_chr[[chr]]
    pos <- findInterval(anno$center[idx], DNA$chrompos[chr_bins])
    pos[pos == 0] <- NA_integer_
    anno$bin_row[idx] <- ifelse(is.na(pos), NA_integer_, chr_bins[pos])
  }
  anno <- anno[!is.na(anno$bin_row), ]

  # Merge with row IDs
  anno <- merge(anno, data.frame(hgnc_symbol = RNA.mat$hgnc_symbol,
                                  row_id = seq_len(nrow(RNA.mat)), stringsAsFactors = FALSE))

  RNA_vals <- as.matrix(RNA.mat[, 8:ncol(RNA.mat), drop = FALSE])
  result <- matrix(NA_real_, nrow = nrow(DNA), ncol = ncol(RNA_vals))
  colnames(result) <- colnames(RNA_vals)

  # Aggregate by bin
  anno <- anno[order(anno$bin_row), ]
  if (nrow(anno) > 0) {
    bin_rle <- rle(anno$bin_row)
    start_i <- 1
    for (i in seq_along(bin_rle$values)) {
      end_i <- start_i + bin_rle$lengths[i] - 1
      rows <- anno$row_id[start_i:end_i]
      result[bin_rle$values[i], ] <- if (length(rows) == 1) {
        RNA_vals[rows, ]
      } else {
        matrixStats::colMedians(RNA_vals[rows, , drop = FALSE])
      }
      start_i <- end_i + 1
    }
  }

  # Fill empty bins
  empty <- which(rowSums(!is.na(result)) == 0)
  if (length(empty) > 0) {
    avail <- which(rowSums(!is.na(result)) > 0)
    fill <- sapply(empty, function(i) avail[which.min(abs(avail - i))])
    result[empty, ] <- result[fill, , drop = FALSE]
  }

  list(DNA.adj = DNA, RNA.adj = cbind(DNA[, 1:3], result))
}

#' Run CNV analysis for tumor/normal classification
#' @param rawmat Raw count matrix
#' @param id.type Gene ID type
#' @param cell.line Cell line mode
#' @param ngene.chr Minimum genes per chromosome
#' @param min.gene.per.cell Minimum genes per cell
#' @param LOW.DR Low detection rate
#' @param UP.DR High detection rate
#' @param win.size Window size
#' @param norm.cell.names Known normal cells
#' @param KS.cut KS test cutoff
#' @param distance Distance metric
#' @param n.cores Number of cores
#' @param full.anno Gene annotations
#' @param DNA.hg20 Genomic bins
#' @param cyclegenes Cell cycle genes
#' @keywords internal
#' @noRd
.run_copykat_internal <- function(rawmat, id.type = "S", cell.line = FALSE,
                                   ngene.chr = 5, min.gene.per.cell = 200,
                                   LOW.DR = 0.05, UP.DR = 0.1, win.size = 25,
                                   norm.cell.names = NULL, KS.cut = 0.1,
                                   distance = "euclidean", n.cores = 1,
                                   full.anno, DNA.hg20, cyclegenes) {
  set.seed(1234)
  original_cells <- colnames(rawmat)

  # Filter cells
  message("  [CNV] Filtering cells...")
  genes_per_cell <- colSums(rawmat > 0)
  if (sum(genes_per_cell > min.gene.per.cell) == 0) {
    stop("No cells pass minimum gene filter")
  }
  rawmat <- rawmat[, genes_per_cell >= min.gene.per.cell, drop = FALSE]

  # Filter genes
  det_rate <- rowSums(rawmat > 0) / ncol(rawmat)
  rawmat <- rawmat[det_rate > LOW.DR, , drop = FALSE]
  if (nrow(rawmat) < 7000) UP.DR <- LOW.DR

  # Annotate
  message("  [CNV] Annotating genes...")
  anno_mat <- .copykat_annotateGenes_hg20(rawmat, id.type, full.anno)
  anno_mat <- anno_mat[order(as.numeric(anno_mat$abspos)), ]

  # Remove cell cycle genes
  hla <- anno_mat$hgnc_symbol[grep("^HLA-", anno_mat$hgnc_symbol)]
  to_remove <- which(anno_mat$hgnc_symbol %in% c(cyclegenes[[1]], hla))
  if (length(to_remove) > 0) anno_mat <- anno_mat[-to_remove, ]

  # Filter cells by chromosome coverage
  to_remove <- character(0)
  for (i in 8:ncol(anno_mat)) {
    cell <- cbind(anno_mat$chromosome_name, anno_mat[, i])
    cell <- cell[cell[, 2] != 0, , drop = FALSE]
    n_expressed <- nrow(cell)
    if (n_expressed < 5 ||
        length(rle(cell[, 1])$length) < length(unique(anno_mat$chromosome_name)) ||
        min(rle(cell[, 1])$length) < ngene.chr) {
      to_remove <- c(to_remove, colnames(anno_mat)[i])
    }
  }
  if (length(to_remove) == ncol(anno_mat) - 7) stop("All cells filtered")
  if (length(to_remove) > 0) {
    anno_mat <- anno_mat[, !colnames(anno_mat) %in% to_remove, drop = FALSE]
  }

  # Normalize
  raw_vals <- data.matrix(anno_mat[, 8:ncol(anno_mat)])
  norm_mat <- log(sqrt(raw_vals) + sqrt(raw_vals + 1))
  norm_mat <- apply(norm_mat, 2, function(x) x - mean(x))
  colnames(norm_mat) <- colnames(raw_vals)

  # Smooth with DLM
  message("  [CNV] Smoothing...")
  smooth_fn <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm_mat[, c], model)$s[-1]
    x - mean(x)
  }
  smooth_res <- .safe_mclapply(seq_len(ncol(norm_mat)), smooth_fn, mc.cores = n.cores)
  norm_smooth <- matrix(unlist(smooth_res), ncol = ncol(norm_mat), byrow = FALSE)
  colnames(norm_smooth) <- colnames(norm_mat)

  # Compute baseline
  message("  [CNV] Computing baseline...")
  if (cell.line) {
    base_res <- .copykat_baseline_synthetic(norm_smooth, min.cells = 10, n.cores = n.cores)
    norm_relat <- base_res$expr.relat
    CL <- base_res$cl
    preN <- NULL
  } else if (!is.null(norm.cell.names) && length(norm.cell.names) > 1) {
    valid <- intersect(colnames(norm_smooth), norm.cell.names)
    if (length(valid) == 0) stop("Known normal cells not found")
    basel <- apply(norm_smooth[, valid, drop = FALSE], 1, median)
    d <- .safe_parDist(t(norm_smooth), threads = n.cores, method = "euclidean")
    fit <- hclust(d, method = "ward.D2")
    CL <- cutree(fit, k = 6)
    preN <- norm.cell.names
    norm_relat <- norm_smooth - basel
  } else {
    base_res <- .copykat_baseline_norm_cl(norm_smooth, min.cells = 5, n.cores = n.cores)
    if (base_res$WNS == "unclassified.prediction") {
      base_res <- .copykat_baseline_GMM(norm_smooth, max.normal = 5, mu.cut = 0.05,
                                         Nfraq.cut = 0.99, RE.before = base_res, n.cores = n.cores)
    }
    CL <- base_res$cl
    preN <- base_res$preN
    norm_relat <- norm_smooth - base_res$basel
  }

  # Higher quality genes for segmentation
  det_rate2 <- rowSums(raw_vals > 0) / ncol(raw_vals)
  norm_relat <- norm_relat[det_rate2 >= UP.DR, , drop = FALSE]
  anno_mat2 <- anno_mat[det_rate2 >= UP.DR, ]

  # Filter cells again
  to_remove <- character(0)
  for (i in 8:ncol(anno_mat2)) {
    cell <- cbind(anno_mat2$chromosome_name, anno_mat2[, i])
    cell <- cell[cell[, 2] != 0, , drop = FALSE]
    n_expressed <- nrow(cell)
    if (n_expressed < 5 ||
        length(rle(cell[, 1])$length) < length(unique(anno_mat$chromosome_name)) ||
        min(rle(cell[, 1])$length) < ngene.chr) {
      to_remove <- c(to_remove, colnames(anno_mat2)[i])
    }
  }
  if (length(to_remove) == ncol(norm_relat)) stop("All cells filtered")
  if (length(to_remove) > 0) {
    norm_relat <- norm_relat[, !colnames(norm_relat) %in% to_remove, drop = FALSE]
  }

  CL <- CL[names(CL) %in% colnames(norm_relat)]
  CL <- CL[match(colnames(norm_relat), names(CL))]

  # Segmentation
  message("  [CNV] Segmenting...")
  seg_res <- .copykat_CNA_MCMC(CL, norm_relat, win.size, KS.cut, n.cores)
  if (length(seg_res$breaks) < 25) {
    seg_res <- .copykat_CNA_MCMC(CL, norm_relat, win.size, KS.cut * 0.5, n.cores)
  }
  if (length(seg_res$breaks) < 25) {
    seg_res <- .copykat_CNA_MCMC(CL, norm_relat, win.size, KS.cut * 0.25, n.cores)
  }

  colnames(seg_res$logCNA) <- colnames(norm_relat)
  seg_centered <- apply(seg_res$logCNA, 2, function(x) x - mean(x))
  rna_cnv <- cbind(anno_mat2[, 1:7], seg_centered)

  # Convert to bins
  message("  [CNV] Converting to genomic bins...")
  bin_res <- .copykat_convert_bins_hg20(DNA.hg20, rna_cnv, full.anno, n.cores)
  uber_mat <- data.matrix(bin_res$RNA.adj[, 4:ncol(bin_res$RNA.adj)])

  # Classify
  message("  [CNV] Classifying cells...")
  classify_cells <- function(mat, ref_normal = preN) {
    if (distance == "euclidean") {
      hc <- hclust(.safe_parDist(t(mat), threads = n.cores, method = distance), method = "ward.D")
    } else {
      hc <- hclust(as.dist(1 - cor(mat, method = distance)), method = "ward.D")
    }
    clusters <- cutree(hc, 2)

    # Identify normal cluster based on reference or variance
    if (!is.null(ref_normal) && length(ref_normal) > 0) {
      cl_prop <- sapply(1:2, function(i) {
        cells <- names(clusters)[clusters == i]
        if (length(cells) == 0) return(0)
        length(intersect(cells, ref_normal)) / length(cells)
      })
    } else {
      # If no reference, use variance - lower variance = more likely diploid
      cl_prop <- sapply(1:2, function(i) {
        cells <- names(clusters)[clusters == i]
        if (length(cells) == 0) return(Inf)
        -mean(apply(mat[, cells, drop = FALSE], 2, var))
      })
    }

    pred <- ifelse(clusters == which.max(cl_prop), "diploid", "aneuploid")
    names(pred) <- names(clusters)
    pred
  }

  # First pass
  pred1 <- classify_cells(uber_mat, preN)

  # Baseline adjustment
  diploid_idx <- pred1 == "diploid"
  if (sum(diploid_idx) == 0) {
    # No diploid cells identified - use all cells for baseline
    warning("No diploid cells identified in first pass, using median baseline")
    diploid_idx <- rep(TRUE, length(pred1))
  }
  adj_mat <- uber_mat - apply(uber_mat[, diploid_idx, drop = FALSE], 1, mean)
  adj_mat <- apply(adj_mat, 2, function(x) x - mean(x))
  norm_vals <- adj_mat[, diploid_idx, drop = FALSE]
  cf <- apply(norm_vals, 1, sd)
  base_vals <- apply(norm_vals, 1, mean)

  adj_fn <- function(j) {
    a <- adj_mat[, j]
    a[abs(a - base_vals) <= 0.25 * cf] <- mean(a)
    a
  }
  adj_res <- .safe_mclapply(seq_len(ncol(adj_mat)), adj_fn, mc.cores = n.cores)
  final_mat <- t(t(matrix(unlist(adj_res), ncol = ncol(adj_mat))) -
                   apply(matrix(unlist(adj_res), ncol = ncol(adj_mat)), 2, mean))
  colnames(final_mat) <- colnames(adj_mat)

  # Final classification
  final_pred <- classify_cells(final_mat, preN)

  # Add back filtered cells
  missing <- setdiff(original_cells, names(final_pred))
  if (length(missing) > 0) {
    final_pred <- c(final_pred, setNames(rep("not.defined", length(missing)), missing))
  }

  data.frame(cell.names = names(final_pred), copykat.pred = final_pred,
             stringsAsFactors = FALSE)
}
