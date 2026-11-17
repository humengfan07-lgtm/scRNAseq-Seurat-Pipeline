#######################################################
# ITP 外周血样本 - 重新注释 (使用 ref_Hematopoietic)
# 造血细胞参考数据集，适合外周血PBMC
#######################################################

library(Seurat)
library(SingleR)
library(dplyr)
library(data.table)
library(ggplot2)


cat("========================================\n")
cat("ITP | 外周血重新注释 (ref_Hematopoietic)\n")
cat("========================================\n\n")

# ==================== 1. 加载数据 ====================
biofsci <- readRDS("biofsci03_integrated.rds")
cat(sprintf("加载数据: %d cells, %d genes\n", ncol(biofsci), nrow(biofsci)))

# 恢复聚类
biofsci <- FindNeighbors(biofsci, reduction = "harmony", dims = 1:20)
biofsci <- FindClusters(biofsci, resolution = 0.5)
n_clusters <- nlevels(biofsci$seurat_clusters)
cat(sprintf("聚类数 (resolution=0.5): %d\n", n_clusters))

# ==================== 2. 加载本地参考数据集 ====================
ref_dir <- Sys.getenv("ITP_REF_DIR", unset = "references")

# 加载 ref_Hematopoietic (造血细胞参考)
cat("\n加载 ref_Hematopoietic ...\n")
load(file.path(ref_dir, "ref_Hematopoietic.RData"))
cat("ref_Hematopoietic 标签:\n")
print(table(ref_Hematopoietic[["label.main"]]))

# ==================== 3. SingleR 注释 (ref_Hematopoietic) ====================
cat("\n--- SingleR 注释 ---\n")

test_data <- GetAssayData(biofsci, assay = "RNA", layer = "data")

singler_hema <- SingleR(
  test = test_data,
  ref = ref_Hematopoietic,
  labels = ref_Hematopoietic[["label.main"]],
  clusters = biofsci$seurat_clusters
)

cat("\n注释结果:\n")
print(singler_hema)

# ==================== 4. 整理注释 ====================
clusterAnn <- data.frame(
  cluster = rownames(singler_hema),
  label   = singler_hema[["labels"]],
  stringsAsFactors = FALSE
)
cat("\n集群 -> 细胞类型映射:\n")
print(clusterAnn)

cellAnn <- singler_hema$labels[match(as.character(biofsci$seurat_clusters), rownames(singler_hema))]
biofsci[["cell_type"]] <- cellAnn

cat("\n细胞类型统计 (整体):\n")
print(table(cellAnn))

# ==================== 5. 保存注释表 ====================
fwrite(clusterAnn, "results/09_clusterAnn_hema.csv")
fwrite(data.frame(
  cell_barcode = colnames(biofsci),
  cell_type    = cellAnn,
  stringsAsFactors = FALSE
), "results/09_cellAnn_hema.csv")

cat("\n已保存: results/09_clusterAnn_hema.csv\n")
cat("已保存: results/09_cellAnn_hema.csv\n")

# ==================== 6. Before vs After 统计 ====================
cat("\n--- 治疗前后各细胞类型统计 ---\n")
stats <- table(biofsci$cell_type, biofsci$group)
print(stats)

prop <- prop.table(stats, margin = 2) * 100
cat("\n治疗前后各细胞类型比例 (%):\n")
print(round(prop, 1))

# 比例变化
diff_df <- data.frame(
  CellType   = rownames(stats),
  Before     = as.vector(stats[, "Before"]),
  After      = as.vector(stats[, "After"]),
  Before_pct = round(as.vector(prop[, "Before"]), 1),
  After_pct  = round(as.vector(prop[, "After"]), 1),
  stringsAsFactors = FALSE
)
diff_df[["Change"]] <- round(diff_df[["After_pct"]] - diff_df[["Before_pct"]], 1)
diff_df <- diff_df[order(diff_df[["Change"]], decreasing = TRUE), ]

cat("\n治疗前后变化 (百分点):\n")
print(diff_df[, c("CellType", "Before_pct", "After_pct", "Change")])

fwrite(diff_df, "results/10_cell_type_change_hema.csv")

# ==================== 7. 注释可视化 ====================
cat("\n绘制可视化...\n")

# UMAP 按细胞类型
p1 <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
              group.by = "cell_type", label = TRUE, repel = TRUE,
              label.size = 3)
ggsave("results/10_cellAnn_UMAP_hema.pdf", p1, width = 12, height = 8)

# UMAP 按分组拆分
p2 <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
              group.by = "cell_type", label = FALSE,
              split.by = "group")
ggsave("results/10_cellAnn_UMAP_split_hema.pdf", p2, width = 16, height = 7)

# 细胞比例堆叠柱状图
cell_prop <- as.data.frame(prop.table(
  table(biofsci$cell_type, biofsci$group), margin = 2
))
colnames(cell_prop) <- c("CellType", "Group", "Proportion")

p3 <- ggplot(cell_prop, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion", title = "Cell Type Proportion: Before vs After") +
  theme(legend.title = element_blank())
ggsave("results/10_cell_type_proportion_hema.pdf", p3, width = 12, height = 7)

# 每个患者的比例
patient_prop <- as.data.frame(prop.table(
  table(biofsci$cell_type, biofsci$patient, biofsci$group),
  margin = c(2, 3)
))
colnames(patient_prop) <- c("CellType", "Patient", "Group", "Proportion")

p4 <- ggplot(patient_prop,
             aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Patient, nrow = 1) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion", title = "Cell Type Per Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/10_cell_type_proportion_by_patient_hema.pdf",
       p4 + theme(strip.text = element_text(family = "PingFang SC")),
       width = 18, height = 6, device = grDevices::quartz, type = "pdf")

# 变化热图式柱状图
p5 <- ggplot(diff_df, aes(x = reorder(CellType, Change), y = Change, fill = Change)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(x = "", y = "Change (After - Before, %)", title = "Cell Type Proportion Change") +
  theme(legend.position = "none")
ggsave("results/10_cell_type_change_bar_hema.pdf", p5, width = 10, height = 6)

# 每患者每个细胞类型变化
patient_change <- prop.table(
  table(biofsci$cell_type, biofsci$patient, biofsci$group),
  margin = c(2, 3)
) * 100
patient_change_df <- data.frame(
  CellType = rep(rownames(patient_change), 3),
  Patient  = rep(colnames(patient_change), each = nrow(patient_change)),
  After    = as.vector(patient_change[, , "After"]),
  Before   = as.vector(patient_change[, , "Before"]),
  Change   = as.vector(patient_change[, , "After"]) - as.vector(patient_change[, , "Before"]),
  stringsAsFactors = FALSE
)

p6 <- ggplot(patient_change_df, aes(x = CellType, y = Change, fill = Change)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Patient, nrow = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(x = "", y = "Change (After - Before, %)", title = "Per-Patient Cell Type Change") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))
ggsave("results/10_cell_type_change_per_patient_hema.pdf",
       p6 + theme(strip.text = element_text(family = "PingFang SC")),
       width = 20, height = 7, device = grDevices::quartz, type = "pdf")

# ==================== 8. 保存最终对象 ====================
saveRDS(biofsci, "biofsci04_clustered_hema.rds")
cat("\n已保存: biofsci04_clustered_hema.rds\n")

cat("\n========== 重新注释完毕 ==========\n")
cat("\n输出文件:\n")
cat("  results/09_clusterAnn_hema.csv - 聚类→细胞类型\n")
cat("  results/09_cellAnn_hema.csv   - 每个细胞注释\n")
cat("  results/10_cellAnn_UMAP_hema.pdf - UMAP注释图\n")
cat("  results/10_cell_type_proportion_hema.pdf - 比例图\n")
cat("  results/10_cell_type_change_hema.csv - 变化数据\n")
cat("  biofsci04_clustered_hema.rds - Seurat对象\n")
