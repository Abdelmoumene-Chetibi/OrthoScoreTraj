library(Seurat)
library(OrthoScoreTraj)
library(ggplot2)

data_folder <- "C:/Users/hp/Desktop/R PROJECT/OrthoScoreTraj/real_data_test"
data_dir <- list.dirs(data_folder, recursive = TRUE, full.names = TRUE)
data_dir <- data_dir[grepl("hg19$", data_dir)][1]
raw_counts <- Read10X(data.dir = data_dir)

pbmc <- CreateSeuratObject(counts = raw_counts, project = "10x_PBMC_Real", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, verbose = FALSE)

# 1. Benchmark A: Cytotoxic Cell Detection (AUROC & Speed)
cytotoxic_genes <- c("GZMB", "GZMA", "GZMH", "PRF1", "GNLY", "NKG7", "FGFBP2", "CST7")
confound_genes  <- c("RPL10", "RPL13", "RPS18", "RPS27", "ACTB", "GAPDH", "MALAT1")

pbmc <- AddModuleScore(pbmc, features = list(cytotoxic_genes), name = "Seurat_Cyto", verbose = FALSE)

t0 <- Sys.time()
pbmc <- runOrthoScore(pbmc, target_genes = cytotoxic_genes, confound_genes = confound_genes, score_name = "Ortho_Cyto")
t_run <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)

norm_mat <- GetAssayData(pbmc, layer = "data")
is_true_cytotoxic <- (norm_mat["CD8A", ] > 0 | norm_mat["NKG7", ] > 1.5) & (norm_mat["CD14", ] == 0 & norm_mat["MS4A1", ] == 0)

calc_auroc <- function(scores, labels) {
  r <- rank(scores)
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  round((sum(r[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg), 4)
}

auroc_seurat <- calc_auroc(pbmc$Seurat_Cyto1, is_true_cytotoxic)
auroc_ortho  <- calc_auroc(pbmc$Ortho_Cyto_Score, is_true_cytotoxic)

# 2. Benchmark B: Confounded Signature Decoupling (Exhaustion vs Cytotoxicity on T/NK cells)
t_nk_cells <- colnames(pbmc)[norm_mat["CD3D", ] > 0 | norm_mat["NKG7", ] > 1.5]
pbmc_tnk <- subset(pbmc, cells = t_nk_cells)

exhaustion_genes <- c("TIGIT", "LAG3", "HAVCR2", "TOX", "TOX2", "ENTPD1", "GZMB", "CCL5", "HLA-DRA")
cyto_confound    <- c("NKG7", "GNLY", "PRF1", "GZMA", "GZMH", "FGFBP2", "CST7")

pbmc_tnk <- AddModuleScore(pbmc_tnk, features = list(exhaustion_genes), name = "Seurat_Exh", verbose = FALSE)
pbmc_tnk <- runOrthoScore(pbmc_tnk, target_genes = exhaustion_genes, confound_genes = cyto_confound, score_name = "Ortho_Exh")

tnk_mat <- GetAssayData(pbmc_tnk, layer = "data")
pure_cyto_signal <- as.numeric(colMeans(tnk_mat[c("NKG7", "GNLY", "PRF1"), ]))

leakage_seurat <- round(cor(pbmc_tnk$Seurat_Exh1, pure_cyto_signal), 4)
leakage_ortho  <- round(cor(pbmc_tnk$Ortho_Exh_Score, pure_cyto_signal), 4)
var_reduction  <- round((1 - (leakage_ortho^2 / leakage_seurat^2)) * 100, 1)

# 3. Benchmark C: UMAP-Free Geodesic Trajectory
naive_score <- as.numeric(tnk_mat["CCR7", ] + tnk_mat["LEF1", ])
root_naive <- which.max(naive_score)
t1 <- Sys.time()
pbmc_tnk <- runGeodesicTraj(pbmc_tnk, n_pcs = 15, k = 15, root_cell = root_naive, traj_name = "Geo_PT")
t_traj <- round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 3)

# 4. Generate Publication-Ready Comparison Figure
df_plot <- rbind(
  data.frame(Method = paste0("Seurat::AddModuleScore (r = ", leakage_seurat, ")"),
             Cyto = pure_cyto_signal, Score = scale(pbmc_tnk$Seurat_Exh1)),
  data.frame(Method = paste0("OrthoScoreTraj (r = ", leakage_ortho, ")"),
             Cyto = pure_cyto_signal, Score = scale(pbmc_tnk$Ortho_Exh_Score))
)

p <- ggplot(df_plot, aes(x = Cyto, y = Score, color = Method)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  facet_wrap(~Method, scales = "free_y") +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  labs(
    title = "Confounded Signature Decoupling on Human 10x PBMC T/NK Cells (n = 1,719)",
    subtitle = paste0("Spurious cytotoxic variance reduced by ", var_reduction, "% via Constrained Orthogonalization"),
    x = "Confounding Cytotoxic Signature (NKG7 / GNLY / PRF1)",
    y = "Standardized Exhaustion Score (Z-score)"
  )

ggsave("inst/benchmarks/real_pbmc_decoupling_benchmark.png", plot = p, width = 9, height = 4.5, dpi = 300)

cat("\n================ COMPLETE REAL HUMAN DATA BENCHMARK REPORT ================\n")
cat("1. OrthoScore Runtime (2,700 cells)           :", t_run, "sec\n")
cat("2. Cytotoxic Cell AUROC (Seurat vs Ortho)     :", auroc_seurat, "vs", auroc_ortho, "\n")
cat("3. Confounded Leakage r (Seurat vs Ortho)     :", leakage_seurat, "vs", leakage_ortho, "(-", var_reduction, "% spurious variance)\n")
cat("4. Geodesic Trajectory Runtime (1,719 cells)  :", t_traj, "sec\n")
cat("5. Figure saved to: inst/benchmarks/real_pbmc_decoupling_benchmark.png\n")