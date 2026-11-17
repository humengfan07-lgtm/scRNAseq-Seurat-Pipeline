#######################################################
# ITP ITP 外周血 - 补充注释
# 使用本地 ref_Human_all 细粒度标签 + ITP marker 分析
#######################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(data.table)

cat("========================================\n")
cat("ITP 外周血 | 补充注释 (细粒度 + ITP特征)\n")
cat("========================================\n\n")

biofsci <- readRDS("biofsci04_clustered.rds")
cat(sprintf("加载数据: %d cells\n", ncol(biofsci)))

# ==================== 1. 重新注释 - 细粒度标签 ====================
cat("\n=== 细粒度注释 (ref_Human_all label.fine) ===\n")

ref_dir <- Sys.getenv("ITP_REF_DIR", unset = "references")
load(file.path(ref_dir, "ref_Human_all.RData"))

# 查看有哪些标签类型
cat("\nref_Human_all 可用标签:\n")
cat("  label.main:\n")
print(table(ref_Human_all$label.main))
cat("\n  label.fine:\n")
print(table(ref_Human_all$label.fine))

test_data <- GetAssayData(biofsci, assay = "RNA", layer = "data")

# 用 label.fine 做细粒度注释
singler_fine <- SingleR(
  test = test_data,
  ref = ref_Human_all,
  labels = ref_Human_all$label.fine,
  clusters = biofsci$seurat_clusters
)

cat("\n--- 细粒度注释结果 ---\n")
print(singler_fine)

# 用 label.main 做粗粒度注释做对比
singler_main <- SingleR(
  test = test_data,
  ref = ref_Human_all,
  labels = ref_Human_all$label.main,
  clusters = biofsci$seurat_clusters
)

# 汇总比较
label_comp <- data.frame(
  cluster = rownames(singler_main),
  main = singler_main$labels,
  fine = singler_fine$labels,
  nCells = as.integer(table(biofsci$seurat_clusters)[rownames(singler_fine)]),
  stringsAsFactors = FALSE
)
label_comp <- label_comp[order(as.numeric(label_comp$cluster)), ]
cat("\n--- 粗/细标签对比 ---\n")
print(label_comp)

# ==================== 2. ITP 相关基因特征 ====================
cat("\n\n=== ITP 相关基因表达分析 ===\n")

# ITP/PBMC 关键细胞类型 marker 集
itp_marker_sets <- list(
  "CD4_Tn"         = c("CD3D", "CD4", "CCR7", "SELL", "TCF7", "LEF1"),
  "CD4_Tm"         = c("CD3D", "CD4", "IL7R", "S100A4", "LTB"),
  "CD8_T"          = c("CD3D", "CD8A", "CD8B", "GZMK", "NKG7", "CCL5"),
  "CD8_Tem"        = c("CD3D", "CD8A", "GZMH", "GZMB", "PRF1"),
  "Treg"           = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "IL10"),
  "NK"             = c("NKG7", "GNLY", "KLRF1", "KLRD1", "NCAM1", "FCGR3A"),
  "B_naive"        = c("MS4A1", "CD79A", "IGHM", "IGHD", "CD19", "TCL1A"),
  "B_memory"       = c("MS4A1", "CD27", "IGHG1", "CD79A", "CD79B", "TNFRSF13B"),
  "Plasma"         = c("SDC1", "MZB1", "IGHG1", "XBP1", "JCHAIN", "IRF4"),
  "CD14_Mono"      = c("CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1"),
  "CD16_Mono"      = c("FCGR3A", "CDKN1C", "LST1", "MS4A7", "HES4"),
  "cDC"            = c("FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1"),
  "pDC"            = c("CLEC4C", "LILRA4", "IRF7", "TCF4", "GZMB"),
  "Platelets_MK"   = c("PPBP", "PF4", "GP9", "ITGA2B", "ITGB3", "PBX1"),
  "HSPC"           = c("CD34", "KIT", "GATA2", "CD38", "PROM1")
)

# 检查可用性
all_genes <- rownames(biofsci)
for (ct in names(itp_marker_sets)) {
  avail <- intersect(itp_marker_sets[[ct]], all_genes)
  itp_marker_sets[[ct]] <- avail
  if (length(avail) > 0) {
    cat(sprintf("  %-15s: %d/%d genes\n", ct, length(avail),
                length(itp_marker_sets[[ct]])))
  }
}

# 平均表达矩阵
itp_genes_all <- unique(unlist(itp_marker_sets))
expr_avg <- AverageExpression(biofsci, features = itp_genes_all,
                               assays = "RNA", group.by = "seurat_clusters")
expr_mat <- as.matrix(expr_avg$RNA)
expr_scaled <- t(scale(t(expr_mat)))

# 各聚类top ITP细胞类型得分
cat("\n--- 各聚类 ITP 细胞类型得分 ---\n")
cluster_scores <- data.frame(cluster = colnames(expr_mat), stringsAsFactors = FALSE)

for (ct in names(itp_marker_sets)) {
  genes <- itp_marker_sets[[ct]]
  if (length(genes) > 0) {
    scores <- colMeans(expr_mat[genes, , drop = FALSE])
    cluster_scores[[ct]] <- round(scores, 3)
  }
}

print(cluster_scores, row.names = FALSE)
fwrite(cluster_scores, "results/09d_ITP_celltype_scores.csv")
cat("\n已保存: results/09d_ITP_celltype_scores.csv\n")

# ==================== 3. 综合推荐注释 ====================
cat("\n\n========================================\n")
cat("=== 综合注释建议 (基于细粒度标签 + ITP特征) ===\n")
cat("========================================\n\n")

# 基于细粒度+ITP marker的最终建议
annotation_guide <- data.frame(
  cluster = as.integer(rownames(singler_fine)),
  fine_label = singler_fine$labels,
  main_label = singler_main$labels,
  nCells = as.integer(table(biofsci$seurat_clusters)[rownames(singler_fine)]),
  stringsAsFactors = FALSE
)

# 输出建议列
cat("cluster | nCells | fine_label | main_label | 建议\n")
cat("--------|--------|------------|------------|-----\n")

for (i in 1:nrow(annotation_guide)) {
  cl <- annotation_guide$cluster[i]

  # 根据ITP得分微调
  fine <- annotation_guide$fine_label[i]
  main <- annotation_guide$main_label[i]

  suggestion <- fine  # 默认使用细粒度标签

  cat(sprintf("  %2d     | %5d | %-35s | %-25s | %s\n",
              cl, annotation_guide$nCells[i], fine, main, suggestion))
}

# ==================== 4. 输出可视化 - 细粒度注释 UMAP ====================
# 将 fine label 应用到 Seurat 对象
fine_labels <- singler_fine$labels[match(as.character(biofsci$seurat_clusters), rownames(singler_fine))]
biofsci$cell_type_fine <- fine_labels

# 清理标签名 (替换下划线)
biofsci$cell_type_fine_clean <- gsub("_", " ", biofsci$cell_type_fine)

p_fine <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                  group.by = "cell_type_fine", label = TRUE, repel = TRUE,
                  label.size = 3) +
  ggtitle("细粒度注释 (label.fine)")
ggsave("results/11_cellAnn_fine_UMAP.pdf", p_fine, width = 14, height = 10)

p_fine_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                        group.by = "cell_type_fine", label = FALSE,
                        split.by = "group")
ggsave("results/11_cellAnn_fine_UMAP_by_group.pdf", p_fine_split, width = 16, height = 7)

# 细胞类型比例
cell_prop_fine <- as.data.frame(prop.table(
  table(biofsci$cell_type_fine, biofsci$group), margin = 2
))
colnames(cell_prop_fine) <- c("CellType", "Group", "Proportion")

p_prop_fine <- ggplot(cell_prop_fine, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "Proportion", title = "Cell type proportions (fine) - Before vs After")
ggsave("results/11_cell_type_proportion_fine.pdf", p_prop_fine, width = 12, height = 8)

# 保存
fwrite(annotation_guide, "results/11_annotation_guide.csv")
cat("\n已保存: results/11_annotation_guide.csv\n")

# 更新 Seurat 的 cell_type
biofsci$cell_type <- biofsci$cell_type_fine

saveRDS(biofsci, "biofsci04_clustered.rds")
cat("\n已更新: biofsci04_clustered.rds (cell_type = fine labels)\n")

cat("\n========== 补充注释完成 ==========\n")
cat("\n输出文件:\n")
cat("  results/09d_ITP_celltype_scores.csv - ITP特征得分\n")
cat("  results/11_annotation_guide.csv - 综合注释建议\n")
cat("  results/11_cellAnn_fine_UMAP.pdf - 细粒度注释UMAP\n")
cat("  results/11_cell_type_proportion_fine.pdf - 细胞比例\n")
