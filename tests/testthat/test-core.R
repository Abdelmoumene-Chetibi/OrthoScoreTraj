library(testthat)
library(OrthoScoreTraj)

test_that("Orthogonal scoring decouples confounded signatures and classifies states accurately", {
  set.seed(42)
  n_cells <- 300
  n_genes <- 60
  expr <- matrix(0, nrow = n_genes, ncol = n_cells)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("C", seq_len(n_cells))
  
  # True target cells: 1..100; Low/background cells: 101..300
  true_high <- c(rep(TRUE, 150), rep(FALSE, 150))
  # Confounding signature active across cells 51..200
  conf_high <- c(rep(FALSE, 50), rep(TRUE, 150), rep(FALSE, 100))
  
  # Confounding signature (G9..G16)
  for (g in 9:16) {
    expr[g, conf_high] <- rnorm(150, mean = 3.5, sd = 0.25)
    expr[g, !conf_high] <- pmax(0, rnorm(150, mean = 0.2, sd = 0.1))
  }
  
  # Target signature (G1..G8): high in 1..100, low baseline in 101..300, plus confounding leakage
  for (g in 1:8) {
    expr[g, true_high] <- rnorm(150, mean = 4.5, sd = 0.25)
    expr[g, !true_high] <- rnorm(150, mean = 0.8, sd = 0.20)
    expr[g, conf_high] <- expr[g, conf_high] + 1.5
  }
  
  # Background genes (G17..G60)
  for (g in 17:60) {
    expr[g, ] <- pmax(0, rnorm(n_cells, mean = 0.5, sd = 0.3))
  }
  
  res <- calc_ortho_score(
    expr_mat = expr,
    target_genes = paste0("G", 1:8),
    confound_genes = paste0("G", 9:16)
  )
  
  acc <- mean((res$state == "High") == true_high)
  pos <- res$raw_score > 0
  expect_true(acc > 0.90)
  expect_true(abs(cor(res$ortho_score[pos], res$confound_score[pos])) < 0.05)
  expect_true(abs(cor(res$ortho_score[pos], res$confound_score[pos])) < abs(cor(res$raw_score[pos], res$confound_score[pos])))
  expect_equal(length(res$n_detected), n_cells)
})

test_that("Geodesic trajectory and dynamic gene testing work accurately", {
  set.seed(123)
  n_cells <- 200
  n_genes <- 80
  pt_true <- seq(0, 1, length.out = n_cells)
  
  expr <- matrix(0, nrow = n_genes, ncol = n_cells)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("C", seq_len(n_cells))
  
  for (g in 1:25) expr[g, ] <- 5 * pt_true + rnorm(n_cells, sd = 0.2)
  for (g in 26:50) expr[g, ] <- 5 * (1 - pt_true) + rnorm(n_cells, sd = 0.2)
  for (g in 51:65) expr[g, ] <- 4 * sin(pi * pt_true) + rnorm(n_cells, sd = 0.2)
  for (g in 66:80) expr[g, ] <- rnorm(n_cells, mean = 1.0, sd = 0.3)
  
  traj <- infer_geodesic_traj(expr, n_pcs = 3, k = 20, root_cell = 1)
  expect_true(abs(cor(traj$pseudotime, pt_true)) > 0.85)
  expect_true(is.matrix(traj$transition_matrix))
  
  dyn <- test_dynamic_genes(expr, traj$pseudotime, genes = c("G1", "G26", "G51", "G70", "G75"))
  expect_true(dyn$gene[1] %in% c("G1", "G26", "G51"))
  expect_true(dyn$fdr[1] < 1e-10)
})
