#######################################################
# ITP 单细胞测序分析 - Step3: 数据整合与降维
# 标准化 → 高变基因 → PCA → Harmony批次校正 → UMAP/tSNE
#######################################################

library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(harmony)

set.seed(2024)

cat("========================================
")
cat("ITP | Step3: 数据整合与降维
")
cat("========================================

")

# ==================== 1. 加载过滤后数据 ====================
biofsci <- readRDS("biofsci02_filtered.rds")
cat(sprintf("加载数据: %d cells, %d genes
", ncol(biofsci), nrow(biofsci)))

# ==================== 2. 数据标准化 ====================
cat("

标准化 (LogNormalize)...
")
biofsci <- NormalizeData(
  biofsci,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ==================== 3. 高变基因筛选 ====================
cat("

筛选高变基因 (VST, n=2000)...
")
biofsci <- FindVariableFeatures(
  biofsci,
  selection.method = "vst",
  nfeatures = 2000
)

top10 <- head(VariableFeatures(biofsci), 10)
cat("Top10 高变基因:
")
print(top10)

p_var <- VariableFeaturePlot(biofsci)
# 平均表达为 0 的基因无法画在 log10 横轴上，只从图中排除。
p_var$data <- p_var$data[is.finite(p_var$data$mean) & p_var$data$mean > 0, , drop = FALSE]
p_var <- LabelPoints(plot = p_var, points = intersect(top10, rownames(p_var$data)), repel = TRUE)
ggsave("results/01_VariableFeatures.pdf", p_var, width = 10, height = 6)

# ==================== 4. 数据缩放与PCA ====================
cat("

数据缩放 & PCA...
")
# 只对高变基因缩放，避免内存不足 (Seurat v5)
biofsci <- ScaleData(biofsci)

biofsci <- RunPCA(
  biofsci,
  features = VariableFeatures(object = biofsci),
  npcs = 50,
  verbose = FALSE
)

# PCA可视化
p_pca_group <- DimPlot(biofsci, reduction = "pca", group.by = "group", pt.size = 0.8)
ggsave("results/02_PCA_by_group.pdf", p_pca_group, width = 8, height = 6)

p_pca_sample <- DimPlot(biofsci, reduction = "pca", group.by = "sample_id", pt.size = 0.5)
ggsave("results/02_PCA_by_sample.pdf", p_pca_sample, width = 10, height = 7)

p_pca_patient <- DimPlot(biofsci, reduction = "pca", group.by = "patient", pt.size = 0.8)
patient_label_theme <- theme(legend.text = element_text(family = "PingFang SC"),
                             strip.text = element_text(family = "PingFang SC"))
ggsave("results/02_PCA_by_patient.pdf", p_pca_patient + patient_label_theme,
       width = 8, height = 6, device = grDevices::quartz, type = "pdf")

# ElbowPlot
p_elbow <- ElbowPlot(biofsci, ndims = 50)
ggsave("results/02_ElbowPlot.pdf", p_elbow, width = 8, height = 5)

# ==================== 5. Harmony 批次校正 ====================
cat("

Harmony 批次校正 (by sample)...
")
biofsci <- RunHarmony(
  biofsci,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  reduction.save = "harmony",
  theta = 2,
  verbose = TRUE
)

# Harmony 后的 PCA
p_harmony_group <- DimPlot(biofsci, reduction = "harmony",
                           group.by = "group", pt.size = 0.8)
ggsave("results/03_Harmony_by_group.pdf", p_harmony_group, width = 8, height = 6)

p_harmony_sample <- DimPlot(biofsci, reduction = "harmony",
                            group.by = "sample_id", pt.size = 0.5)
ggsave("results/03_Harmony_by_sample.pdf", p_harmony_sample, width = 10, height = 7)

# ==================== 6. 非线性降维可视化 ====================
cat("

UMAP & tSNE 降维...
")

# 确定使用的维度数 (默认20, 可根据ElbowPlot调整)
dims_use <- 1:20
cat(sprintf("  使用维度: %d - %d
", min(dims_use), max(dims_use)))

biofsci <- RunTSNE(biofsci, reduction = "harmony", dims = dims_use)
biofsci <- RunUMAP(biofsci, reduction = "harmony", dims = dims_use)

# --- UMAP 可视化 ---
p_umap_group <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                        group.by = "group")
ggsave("results/04_UMAP_by_group.pdf", p_umap_group, width = 8, height = 6)

p_umap_group_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                              split.by = "group")
ggsave("results/04_UMAP_split_by_group.pdf", p_umap_group_split, width = 14, height = 7)

p_umap_sample <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                         group.by = "sample_id")
ggsave("results/04_UMAP_by_sample.pdf", p_umap_sample, width = 10, height = 7)

p_umap_sample_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                               split.by = "sample_id", ncol = 3)
ggsave("results/04_UMAP_split_by_sample.pdf", p_umap_sample_split,
       width = 18, height = 12)

p_umap_patient <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                          group.by = "patient")
ggsave("results/04_UMAP_by_patient.pdf", p_umap_patient + patient_label_theme,
       width = 8, height = 6, device = grDevices::quartz, type = "pdf")

p_umap_patient_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                                split.by = "patient")
ggsave("results/04_UMAP_split_by_patient.pdf", p_umap_patient_split + patient_label_theme,
       width = 16, height = 6, device = grDevices::quartz, type = "pdf")

# --- tSNE 可视化 ---
p_tsne_group <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.5,
                        group.by = "group")
ggsave("results/05_tSNE_by_group.pdf", p_tsne_group, width = 8, height = 6)

p_tsne_group_split <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.5,
                              split.by = "group")
ggsave("results/05_tSNE_split_by_group.pdf", p_tsne_group_split, width = 14, height = 7)

p_tsne_sample <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.3,
                         group.by = "sample_id")
ggsave("results/05_tSNE_by_sample.pdf", p_tsne_sample, width = 10, height = 7)

p_tsne_sample_split <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.3,
                               split.by = "sample_id", ncol = 3)
ggsave("results/05_tSNE_split_by_sample.pdf", p_tsne_sample_split,
       width = 18, height = 12)

# ==================== 7. 保存整合后数据 ====================
saveRDS(biofsci, "biofsci03_integrated.rds")
cat("

已保存: biofsci03_integrated.rds
")
cat(sprintf("  最终数据: %d cells, %d genes
", ncol(biofsci), nrow(biofsci)))

cat("
========== Step3: 数据整合与降维 完毕 ==========
")
cat("输出文件:
")
cat("  results/01_VariableFeatures.pdf - 高变基因图
")
cat("  results/02_PCA_*.pdf - PCA可视化
")
cat("  results/02_ElbowPlot.pdf - 肘图
")
cat("  results/03_Harmony_*.pdf - Harmony校正后
")
cat("  results/04_UMAP_*.pdf - UMAP可视化
")
cat("  results/05_tSNE_*.pdf - tSNE可视化
")
cat("  biofsci03_integrated.rds - 整合后Seurat对象
")
