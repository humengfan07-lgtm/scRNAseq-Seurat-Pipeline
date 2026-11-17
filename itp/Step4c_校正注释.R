# Step4c：用本批聚类的 marker 和参考注释人工审核。
library(Seurat)

biofsci <- readRDS("biofsci04_clustered.rds")
clusters <- as.character(biofsci$seurat_clusters)
cluster_ids <- sort(unique(clusters))

pbmc_markers <- list(
  "CD14_Monocyte"  = c("CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1", "S100A12"),
  "CD16_Monocyte"  = c("FCGR3A", "CDKN1C", "LST1", "MS4A7", "IFITM3"),
  "cDC"            = c("FCER1A", "CLEC10A", "CD1C", "HLA-DQA1", "HLA-DQB1"),
  "pDC"            = c("CLEC4C", "LILRA4", "IRF7", "TCF4"),
  "CD4_T"          = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "SELL", "LEF1"),
  "CD8_T"          = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "NKG7", "CCL5"),
  "NK"             = c("NKG7", "GNLY", "KLRF1", "KLRD1", "PRF1", "NCAM1"),
  "B_cell"         = c("MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "BANK1"),
  "Plasma"         = c("SDC1", "MZB1", "IGHG1", "JCHAIN", "XBP1", "IRF4"),
  "Platelet_MK"    = c("PPBP", "PF4", "GP9", "ITGA2B", "ITGB3", "VWF", "PBX1"),
  "Erythroid"      = c("HBB", "HBA1", "HBA2", "HBD", "AHSP", "GYPA"),
  "HSPC"           = c("CD34", "KIT", "CD38", "GATA2", "PROM1"),
  "Proliferating"  = c("MKI67", "TOP2A", "PCNA", "STMN1", "TYMS")
)

# 各聚类平均表达
all_marker_genes <- unique(unlist(pbmc_markers))
all_marker_genes <- intersect(all_marker_genes, rownames(biofsci))

avg_expr_all <- AverageExpression(biofsci, features = all_marker_genes,
                                   assays = "RNA", group.by = "seurat_clusters")
expr_mat_all <- as.matrix(avg_expr_all$RNA)

cat("\nPBMC 经典 marker 在各聚类的平均表达:\n")
for (ct in names(pbmc_markers)) {
  genes <- intersect(pbmc_markers[[ct]], rownames(expr_mat_all))
  if (length(genes) > 0) {
    scores <- colMeans(expr_mat_all[genes, , drop = FALSE])
    top3 <- names(sort(scores, decreasing = TRUE)[1:3])
    cat(sprintf("  %-18s -> 聚类 %s\n", ct, paste(top3, collapse=", ")))
  }
}


write.csv(data.frame(gene = rownames(expr_mat_all), expr_mat_all, check.names = FALSE),
          "results/12_marker_expression_by_cluster.csv", row.names = FALSE)

review <- data.frame(cluster = cluster_ids,
                     n_cells = as.integer(table(clusters)[cluster_ids]),
                     cell_type = "")
if (!file.exists("cluster_annotations.csv")) {
  write.csv(review, "cluster_annotations.csv", row.names = FALSE)
}
message("对照 results/08_clusterMarkers.csv 和参考注释，人工填写 cluster_annotations.csv 的 cell_type。")
