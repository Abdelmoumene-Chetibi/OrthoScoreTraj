# Internal helper to extract expression matrix from Seurat, SingleCellExperiment, or matrix
.extract_expr_matrix <- function(object, assay = NULL, slot_layer = "data") {
  if (inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required for Seurat objects.")
    }
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    mat <- tryCatch(
      Seurat::GetAssayData(object, assay = assay, layer = slot_layer),
      error = function(e) Seurat::GetAssayData(object, assay = assay, slot = slot_layer)
    )
    return(mat)
  } else if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required for SingleCellExperiment objects.")
    }
    assay_name <- if (is.null(assay)) "logcounts" else assay
    if (!(assay_name %in% SummarizedExperiment::assayNames(object))) {
      assay_name <- SummarizedExperiment::assayNames(object)[1]
    }
    return(SummarizedExperiment::assay(object, assay_name))
  } else if (is.matrix(object) || inherits(object, "Matrix")) {
    return(object)
  } else {
    stop("Unsupported object class. Provide a Seurat, SingleCellExperiment, or matrix object.")
  }
}

# Internal helper to harmonize gene symbols (e.g., underscore vs dash in Seurat)
.match_gene_names <- function(genes, available_genes) {
  if (is.null(genes)) return(NULL)
  matched <- intersect(genes, available_genes)
  if (length(matched) < 3) {
    alt_genes <- gsub("_", "-", genes)
    matched_alt <- intersect(alt_genes, available_genes)
    if (length(matched_alt) > length(matched)) {
      return(matched_alt)
    }
  }
  return(genes)
}

#' Run OrthoScore Directly on Seurat, SingleCellExperiment, or Matrix Objects
#'
#' @param object A \code{Seurat}, \code{SingleCellExperiment}, or expression matrix.
#' @param target_genes Character vector of target signature genes.
#' @param confound_genes Optional character vector of confounding signature genes.
#' @param score_name Prefix for storing results in object metadata. Default is "OrthoScore".
#' @param assay Assay name to extract (defaults to active assay in Seurat or "logcounts" in SCE).
#' @param adjust_dropout Logical; whether to adjust for library-size dropout. Default is TRUE.
#' @param positive_ortho_only Logical; whether to apply constrained positive orthogonalization. Default is TRUE.
#'
#' @return The updated \code{Seurat} or \code{SingleCellExperiment} object with added metadata columns, or a list if a matrix was supplied.
#' @export
runOrthoScore <- function(object,
                          target_genes,
                          confound_genes = NULL,
                          score_name = "OrthoScore",
                          assay = NULL,
                          adjust_dropout = TRUE,
                          positive_ortho_only = TRUE) {
  expr_mat <- .extract_expr_matrix(object, assay = assay)
  avail <- rownames(expr_mat)
  target_genes <- .match_gene_names(target_genes, avail)
  confound_genes <- .match_gene_names(confound_genes, avail)
  
  res <- calc_ortho_score(
    expr_mat = expr_mat,
    target_genes = target_genes,
    confound_genes = confound_genes,
    adjust_dropout = adjust_dropout,
    positive_ortho_only = positive_ortho_only,
    classify_gmm = TRUE
  )
  
  if (inherits(object, "Seurat")) {
    object[[paste0(score_name, "_Score")]] <- res$ortho_score
    object[[paste0(score_name, "_State")]] <- res$state
    object[[paste0(score_name, "_Prob")]] <- res$posterior_high
    return(object)
  } else if (inherits(object, "SingleCellExperiment")) {
    SummarizedExperiment::colData(object)[[paste0(score_name, "_Score")]] <- res$ortho_score
    SummarizedExperiment::colData(object)[[paste0(score_name, "_State")]] <- res$state
    SummarizedExperiment::colData(object)[[paste0(score_name, "_Prob")]] <- res$posterior_high
    return(object)
  }
  return(res)
}

#' Run UMAP-Free Geodesic Trajectory Directly on Seurat, SingleCellExperiment, or Matrix Objects
#'
#' @param object A \code{Seurat}, \code{SingleCellExperiment}, or expression matrix.
#' @param n_pcs Number of PCs for k-NN manifold construction. Default is 20.
#' @param k Number of nearest neighbors. Default is 15.
#' @param root_cell Optional root cell index or barcode name.
#' @param traj_name Column name for storing pseudotime in metadata. Default is "Geodesic_Pseudotime".
#' @param assay Assay name to extract.
#'
#' @return The updated \code{Seurat} or \code{SingleCellExperiment} object with pseudotime stored in metadata, or a list if a matrix was supplied.
#' @export
runGeodesicTraj <- function(object,
                            n_pcs = 20,
                            k = 15,
                            root_cell = NULL,
                            traj_name = "Geodesic_Pseudotime",
                            assay = NULL) {
  expr_mat <- .extract_expr_matrix(object, assay = assay)
  traj <- infer_geodesic_traj(expr_mat = expr_mat, n_pcs = n_pcs, k = k, root_cell = root_cell)
  
  if (inherits(object, "Seurat")) {
    object[[traj_name]] <- traj$pseudotime
    return(object)
  } else if (inherits(object, "SingleCellExperiment")) {
    SummarizedExperiment::colData(object)[[traj_name]] <- traj$pseudotime
    return(object)
  }
  return(traj)
}