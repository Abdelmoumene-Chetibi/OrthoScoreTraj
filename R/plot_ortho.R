#' Plot Orthogonalized Signature Separation and Geodesic Trajectory
#'
#' @param ortho_res Output list from \code{calc_ortho_score}.
#' @param traj_res Optional output list from \code{infer_geodesic_traj}.
#'
#' @return Invisibly returns NULL; produces base R publication-ready diagnostic plots.
#' @export
plot_ortho_diagnostics <- function(ortho_res, traj_res = NULL) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (!is.null(traj_res)) {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
  }
  
  # Panel 1: Orthogonalized Score vs Confounding Score colored by GMM state
  cols <- ifelse(ortho_res$state == "High", "#D95F02", "#1B9E77")
  x_val <- if (!is.null(ortho_res$confound_score)) ortho_res$confound_score else ortho_res$n_detected
  x_lab <- if (!is.null(ortho_res$confound_score)) "Confounding Signature Score" else "Detected Genes"
  
  graphics::plot(
    x_val, ortho_res$ortho_score,
    col = grDevices::adjustcolor(cols, alpha.f = 0.75),
    pch = 16, cex = 1.1,
    xlab = x_lab,
    ylab = "Orthogonalized Target Score",
    main = "Orthogonal Signature Separation (GMM)"
  )
  graphics::abline(h = 0, lty = 2, col = "gray40")
  graphics::legend("topright", legend = c("High State", "Low State"), col = c("#D95F02", "#1B9E77"), pch = 16, bty = "n")
  
  # Panel 2: Geodesic Pseudotime vs Orthogonal Score
  if (!is.null(traj_res)) {
    graphics::plot(
      traj_res$pseudotime, ortho_res$ortho_score,
      col = grDevices::adjustcolor(cols, alpha.f = 0.75),
      pch = 16, cex = 1.1,
      xlab = "Diffusion Geodesic Pseudotime (UMAP-Free)",
      ylab = "Orthogonalized Target Score",
      main = "Trajectory Progression"
    )
    if (length(unique(traj_res$pseudotime)) >= 4 && stats::sd(ortho_res$ortho_score) > 1e-8) {
      sm <- tryCatch(
        stats::smooth.spline(traj_res$pseudotime, ortho_res$ortho_score, spar = 0.7),
        error = function(e) NULL
      )
      if (!is.null(sm)) graphics::lines(sm, col = "#7570B3", lwd = 3)
    }
  }
  
  invisible(NULL)
}

