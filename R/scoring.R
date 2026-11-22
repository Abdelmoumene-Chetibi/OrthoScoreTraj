#' Compute Dropout-Adjusted and Orthogonalized Gene Signature Scores
#'
#' @param expr_mat A numeric matrix or sparse dgCMatrix (genes x cells) of log-normalized expression.
#' @param target_genes Character vector of target signature gene names.
#' @param confound_genes Optional character vector of confounding signature gene names.
#' @param adjust_dropout Logical; whether to regress out technical library-size dropout. Default is TRUE.
#' @param positive_ortho_only Logical; whether to apply constrained positive orthogonalization. Default is TRUE.
#' @param classify_gmm Logical; whether to fit a 2-component Gaussian Mixture Model to classify states. Default is TRUE.
#'
#' @return A list containing \code{raw_score}, \code{ortho_score}, \code{confound_score}, \code{beta}, \code{state}, \code{posterior_high}, and \code{n_detected}.
#' @export
calc_ortho_score <- function(expr_mat,
                             target_genes,
                             confound_genes = NULL,
                             adjust_dropout = TRUE,
                             positive_ortho_only = TRUE,
                             classify_gmm = TRUE) {
  if (!is.matrix(expr_mat) && !inherits(expr_mat, "Matrix")) {
    stop("expr_mat must be a matrix or sparse Matrix.")
  }
  all_genes <- rownames(expr_mat)
  if (is.null(all_genes)) stop("expr_mat must have rownames corresponding to gene symbols.")
  
  target_idx <- match(intersect(target_genes, all_genes), all_genes)
  if (length(target_idx) < 2) {
    stop("Fewer than 2 target_genes were found in expr_mat.")
  }
  
  target_sub <- as.matrix(expr_mat[target_idx, , drop = FALSE])
  raw_target <- as.numeric(colMeans(target_sub))
  
  # Compute per-cell detected genes (used both for dropout adjustment and diagnostics)
  if (inherits(expr_mat, "dgCMatrix")) {
    n_detected <- as.numeric(diff(expr_mat@p))
  } else {
    n_detected <- as.numeric(colSums(expr_mat > 0))
  }
  det_rate <- n_detected / nrow(expr_mat)
  
  # Technical dropout adjustment
  if (adjust_dropout && sd(det_rate) > 1e-8 && sd(raw_target) > 1e-8) {
    fit_drop <- lm(raw_target ~ det_rate)
    adj_target <- raw_target - coef(fit_drop)[2] * (det_rate - mean(det_rate))
  } else {
    adj_target <- raw_target
  }
  
  # Zero-inflated constrained orthogonalization
  confound_score <- rep(0, ncol(expr_mat))
  beta <- 0
  if (!is.null(confound_genes)) {
    conf_idx <- match(intersect(confound_genes, all_genes), all_genes)
    if (length(conf_idx) >= 2) {
      conf_sub <- as.matrix(expr_mat[conf_idx, , drop = FALSE])
      confound_score <- as.numeric(colMeans(conf_sub))
      
      pos_mask <- (raw_target > 0)
      if (sum(pos_mask) >= 10 && sd(confound_score[pos_mask]) > 1e-8) {
        cov_tc <- cov(adj_target[pos_mask], confound_score[pos_mask])
        var_c  <- var(confound_score[pos_mask])
        beta <- cov_tc / var_c
        if (positive_ortho_only && beta < 0) beta <- 0
        
        ortho_score <- adj_target
        ortho_score[pos_mask] <- adj_target[pos_mask] - beta * (confound_score[pos_mask] - mean(confound_score[pos_mask]))
      } else {
        ortho_score <- adj_target
      }
      
      zero_mask <- !pos_mask
      if (any(zero_mask) && any(pos_mask)) {
        floor_val <- min(ortho_score[pos_mask]) - 0.1 * sd(ortho_score)
        ortho_score[zero_mask] <- pmin(ortho_score[zero_mask], floor_val)
      }
    } else {
      warning("Fewer than 2 confound_genes found; skipping orthogonalization.")
      ortho_score <- adj_target
    }
  } else {
    ortho_score <- adj_target
  }
  
  # Two-component Gaussian Mixture Model classification
  state <- rep("Unclassified", length(ortho_score))
  posterior_high <- rep(NA_real_, length(ortho_score))
  if (classify_gmm && length(unique(ortho_score)) > 5) {
    if (requireNamespace("mclust", quietly = TRUE)) {
      gmm <- tryCatch(
        mclust::Mclust(ortho_score, G = 2, verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(gmm) && !is.null(gmm$z)) {
        high_comp <- which.max(gmm$parameters$mean)
        posterior_high <- gmm$z[, high_comp]
        state <- ifelse(posterior_high >= 0.5, "High", "Low")
      }
    }
    if (all(state == "Unclassified")) {
      med_val <- median(ortho_score)
      state <- ifelse(ortho_score > med_val, "High", "Low")
      posterior_high <- as.numeric(ortho_score > med_val)
    }
  }
  
  list(
    raw_score = raw_target,
    ortho_score = ortho_score,
    confound_score = confound_score,
    beta = beta,
    state = state,
    posterior_high = posterior_high,
    n_detected = n_detected
  )
}