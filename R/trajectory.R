#' Infer UMAP-Free Geodesic Pseudotime Trajectory
#'
#' @param expr_mat A numeric matrix or sparse dgCMatrix (genes x cells) of log-normalized expression.
#' @param n_pcs Number of principal components for k-NN manifold construction. Default is 20.
#' @param k Number of nearest neighbors in the geodesic graph. Default is 15.
#' @param n_hvgs Number of highly variable genes to select for fast PCA. Default is 2000.
#' @param root_cell Optional integer index or cell barcode specifying the trajectory root. If NULL, inferred automatically at the 90th percentile of transcriptional diversity (doublet-robust CytoTRACE principle).
#'
#' @return A list containing \code{pseudotime} (scaled 0 to 1), \code{root_cell}, \code{pca_coords}, and \code{transition_matrix}.
#' @export
infer_geodesic_traj <- function(expr_mat, n_pcs = 20, k = 15, n_hvgs = 2000, root_cell = NULL) {
  n_cells <- ncol(expr_mat)
  if (n_cells < 10) stop("At least 10 cells are required for trajectory inference.")
  
  # Fast HVG selection
  if (nrow(expr_mat) > n_hvgs) {
    if (inherits(expr_mat, "dgCMatrix")) {
      gene_means <- Matrix::rowMeans(expr_mat)
      gene_vars <- Matrix::rowMeans(expr_mat^2) - gene_means^2
    } else {
      gene_means <- rowMeans(expr_mat)
      gene_vars <- apply(expr_mat, 1, var)
    }
    top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvgs)]
    sub_mat <- as.matrix(expr_mat[top_idx, , drop = FALSE])
  } else {
    sub_mat <- as.matrix(expr_mat)
  }
  
  vars <- apply(sub_mat, 1, var)
  sub_mat <- sub_mat[vars > 1e-8, , drop = FALSE]
  
  n_pcs_use <- min(n_pcs, ncol(sub_mat) - 1, nrow(sub_mat) - 1)
  if (requireNamespace("irlba", quietly = TRUE) && n_pcs_use < min(dim(sub_mat)) - 2) {
    pca <- irlba::prcomp_irlba(t(sub_mat), n = n_pcs_use, center = TRUE, scale. = TRUE)
    coords <- pca$x
  } else {
    pca <- prcomp(t(sub_mat), center = TRUE, scale. = TRUE, rank. = n_pcs_use)
    coords <- pca$x[, seq_len(n_pcs_use), drop = FALSE]
  }
  
  # Doublet-robust biological root inference (90th percentile of detected genes)
  if (is.null(root_cell)) {
    if (inherits(expr_mat, "dgCMatrix")) {
      n_det <- diff(expr_mat@p)
    } else {
      n_det <- colSums(expr_mat > 0)
    }
    target_det <- stats::quantile(n_det, probs = 0.90, names = FALSE)
    root_idx <- which.min(abs(n_det - target_det))
  } else if (is.character(root_cell)) {
    root_idx <- match(root_cell, colnames(expr_mat))
    if (is.na(root_idx)) stop("Specified root_cell barcode not found in expr_mat.")
  } else {
    root_idx <- as.integer(root_cell)
  }
  
  k_use <- min(k, n_cells - 1)
  P <- fast_knn_transition_cpp(coords, k = k_use)
  pt <- fast_geodesic_dist_cpp(P, root_idx_1based = root_idx, alpha = 0.9, n_steps = 40L)
  
  list(
    pseudotime = as.numeric(pt),
    root_cell = root_idx,
    pca_coords = coords,
    transition_matrix = P
  )
}

#' Test Dynamic Gene Expression Along Geodesic Pseudotime
#'
#' Fits Generalized Additive Models (mgcv::gam with REML) across pseudotime,
#' incorporating donor/patient random effects when supplied.
#'
#' @param expr_mat Expression matrix (genes x cells).
#' @param pseudotime Numeric vector of pseudotime values.
#' @param genes Character vector of genes to test.
#' @param patient_id Optional vector of donor/patient IDs for random-effect modeling.
#' @param df Degrees of freedom for smooth basis. Default is 3.
#'
#' @return A data.frame sorted by FDR with columns \code{gene}, \code{p_value}, \code{fdr}, and \code{deviance_explained}.
#' @export
test_dynamic_genes <- function(expr_mat, pseudotime, genes, patient_id = NULL, df = 3) {
  valid_genes <- intersect(genes, rownames(expr_mat))
  if (length(valid_genes) == 0) stop("None of the specified genes were found in expr_mat.")
  
  n_unique <- length(unique(pseudotime))
  df_use <- min(df, max(1, n_unique - 2), length(pseudotime) - 2)
  if (df_use < 1) stop("Insufficient unique pseudotime values for dynamic modeling.")
  
  has_mgcv <- requireNamespace("mgcv", quietly = TRUE) && (df_use >= 2)
  has_patient <- !is.null(patient_id) && length(unique(patient_id)) > 1
  if (has_patient) patient_fac <- as.factor(patient_id)
  
  res_list <- lapply(valid_genes, function(g) {
    y <- as.numeric(expr_mat[g, ])
    if (sd(y) < 1e-8) {
      return(data.frame(gene = g, p_value = 1.0, deviance_explained = 0.0, stringsAsFactors = FALSE))
    }
    if (has_mgcv) {
      if (has_patient) {
        fit <- mgcv::gam(y ~ s(pseudotime, k = df_use + 1) + s(patient_fac, bs = "re"), method = "REML")
      } else {
        fit <- mgcv::gam(y ~ s(pseudotime, k = df_use + 1), method = "REML")
      }
      sm <- summary(fit)
      pval <- sm$s.table[1, "p-value"]
      dev_exp <- sm$dev.expl
    } else {
      fit <- lm(y ~ poly(pseudotime, degree = df_use))
      sm <- summary(fit)
      fstat <- sm$fstatistic
      pval <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      dev_exp <- sm$r.squared
    }
    data.frame(gene = g, p_value = as.numeric(pval), deviance_explained = as.numeric(dev_exp), stringsAsFactors = FALSE)
  })
  
  out <- do.call(rbind, res_list)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$fdr, out$p_value), c("gene", "p_value", "fdr", "deviance_explained")]
  rownames(out) <- NULL
  out
}

