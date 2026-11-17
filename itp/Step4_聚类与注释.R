#######################################################
# ITP 单细胞测序分析 - Step4: 聚类与注释
# 聚类 → 标记基因 → SingleR自动注释 → 输出供人工审核
#######################################################

library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(clustree)
library(data.table)
library(presto)
library(SingleR)
# 参考库直接从本地 .RData 文件读取。

set.seed(2024)

FCfilter   <- 0.25
adjPFilter <- 0.05

cat("========================================
")
cat("ITP | Step4: 聚类与注释
")
cat("========================================

")

# ==================== 1. 加载整合后数据 ====================
biofsci <- readRDS("biofsci03_integrated.rds")
cat(sprintf("加载数据: %d cells, %d genes
", ncol(biofsci), nrow(biofsci)))

# ==================== 2. 细胞聚类 ====================
cat("

=== 细胞聚类 ===
")

# 使用维度
dims_use <- 1:20
cat(sprintf("  使用维度: %d - %d
", min(dims_use), max(dims_use)))

biofsci <- FindNeighbors(biofsci, reduction = "harmony", dims = dims_use)

# 多分辨率尝试
cat("  尝试多个分辨率: 0.3, 0.5, 0.8, 1.0, 1.2
")
for (res in c(0.3, 0.5, 0.8, 1.0, 1.2)) {
  biofsci <- FindClusters(biofsci, resolution = res)
}

# Clustree 可视化
p_tree <- clustree(biofsci, prefix = "RNA_snn_res.")
ggsave("results/06_Clustree.pdf", p_tree, width = 10, height = 10)
cat("  已保存: results/06_Clustree.pdf
")

# 选择最终分辨率 (默认0.5, 可根据clustree调整)
final_res <- 0.5
biofsci$seurat_clusters <- biofsci[[paste0("RNA_snn_res.", final_res)]]
Idents(biofsci) <- "seurat_clusters"

n_clusters <- nlevels(biofsci$seurat_clusters)
cat(sprintf("  最终聚类数 (resolution=%.1f): %d
", final_res, n_clusters))

# 每聚类的细胞数
cat("  各聚类细胞数:
")
print(table(biofsci$seurat_clusters))

# ==================== 3. 聚类可视化 ====================
cat("

绘制聚类可视化...
")

# UMAP 聚类图
p_umap_clusters <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                           label = TRUE, label.size = 5, repel = TRUE)
ggsave("results/07_UMAP_clusters.pdf", p_umap_clusters, width = 9, height = 7)

# tSNE 聚类图
p_tsne_clusters <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.5,
                           label = TRUE, label.size = 5, repel = TRUE)
ggsave("results/07_tSNE_clusters.pdf", p_tsne_clusters, width = 9, height = 7)

# 按分组拆分的 UMAP 聚类图
p_umap_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                        split.by = "group", label = TRUE, repel = TRUE)
ggsave("results/07_UMAP_clusters_split_by_group.pdf", p_umap_split,
       width = 16, height = 7)

# 按样本拆分的 UMAP 聚类图
p_umap_sample_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                               split.by = "sample_id", ncol = 3,
                               label = TRUE, repel = TRUE)
ggsave("results/07_UMAP_clusters_split_by_sample.pdf", p_umap_sample_split,
       width = 18, height = 12)

# ==================== 4. 差异标记基因分析（All Markers） ====================
cat("

=== 寻找聚类标记基因 ===
")

# 使用 presto 加速 (wilcoxauc 替代 FindAllMarkers)
cluster_markers <- wilcoxauc(biofsci, "seurat_clusters")

# 计算 avg_log2FC
cluster_markers$avg_logFC <- log2(exp(1)) * cluster_markers$logFC

# 筛选显著标记基因
sig_markers <- cluster_markers %>%
  filter(abs(avg_logFC) > FCfilter & padj < adjPFilter)

fwrite(sig_markers, "results/08_clusterMarkers.csv")
cat(sprintf("  显著标记基因总数: %d
", nrow(sig_markers)))
cat("  各聚类显著标记基因数:
")
print(table(sig_markers$group))

# Top 标记基因热图
top_markers <- cluster_markers %>%
  group_by(group) %>%
  slice_max(order_by = auc, n = 5) %>%
  ungroup()

top_genes <- unique(top_markers$feature)

if (length(top_genes) > 50) {
  top_genes <- head(top_genes, 50)
}

# Step3 只缩放了高变基因。为热图单独缩放这批 marker。
heatmap_obj <- biofsci
heatmap_rna <- heatmap_obj[["RNA"]]
LayerData(heatmap_rna, layer = "scale.data") <- NULL
heatmap_obj[["RNA"]] <- heatmap_rna
heatmap_obj <- ScaleData(heatmap_obj, assay = "RNA", features = top_genes, verbose = FALSE)
p_heatmap <- DoHeatmap(subset(heatmap_obj, downsample = min(100, ncol(heatmap_obj))),
                       features = top_genes, size = 3) +
  NoLegend()
ggsave("results/08_TopMarkers_Heatmap.pdf", p_heatmap,
       width = 16, height = max(8, length(top_genes) * 0.25))
rm(heatmap_obj, heatmap_rna)

cat(sprintf("  已保存: results/08_clusterMarkers.csv
"))
cat(sprintf("  已保存: results/08_TopMarkers_Heatmap.pdf (top %d genes)
", length(top_genes)))

# ==================== 5. SingleR 自动注释 ====================
cat("

=== SingleR 自动注释 ===
")

# 加载参考数据集
ref_dir <- Sys.getenv("ITP_REF_DIR", unset = "references")

ref_human_all     <- file.path(ref_dir, "ref_Human_all.RData")
ref_hematopoietic <- file.path(ref_dir, "ref_Hematopoietic.RData")

if (file.exists(ref_human_all) && file.exists(ref_hematopoietic)) {

  load(ref_human_all)     # 加载 ref_Human_all (这个变量名根据实际数据确定)
  load(ref_hematopoietic) # 加载 ref_Hematopoietic

  # 获取 RNA 表达矩阵用于注释
  test_data <- GetAssayData(biofsci, assay = "RNA", layer = "data")

  # SingleR 注释 - Human_all 作为主参考
  cat("  使用 ref_Human_all 进行注释...
")
  singler_all <- SingleR(
    test = test_data,
    ref = ref_Human_all,
    labels = ref_Human_all$label.main,
    clusters = biofsci$seurat_clusters
  )

  # 输出注释结果
  cat("

--- SingleR 自动注释结果 (ref_Human_all) ---
")
  print(singler_all)

  # 整理聚类注释表
  clusterAnn <- data.frame(
    cluster = rownames(singler_all),
    label   = singler_all$labels,
    stringsAsFactors = FALSE
  )
  # 添加评分 (处理可能的不同格式)
  if (!is.null(singler_all$tuning.scores)) {
    if (is.list(singler_all$tuning.scores)) {
      clusterAnn$score <- round(as.numeric(unlist(singler_all$tuning.scores$first)), 3)
    }
  }

  cat("

集群注释表:
")
  print(clusterAnn)

  # 每个细胞的注释
  cellAnn <- singler_all$labels[match(as.character(biofsci$seurat_clusters), rownames(singler_all))]
  biofsci$cell_type <- cellAnn

  # 保存注释表
  fwrite(clusterAnn, "results/09_clusterAnn.csv")
  fwrite(data.frame(
    cell_barcode = colnames(biofsci),
    cell_type    = cellAnn,
    stringsAsFactors = FALSE
  ), "results/09_cellAnn.csv")

  cat("

已保存注释文件:
")
  cat("  results/09_clusterAnn.csv - 聚类→细胞类型映射表
")
  cat("  results/09_cellAnn.csv   - 每个细胞的注释
")

  # 将 SingleR 标签设置为 Idents
  newLabels <- cellAnn
  names(newLabels) <- colnames(biofsci)
  biofsci$cell_type <- cellAnn

  # ==================== 6. 注释可视化 ====================
  cat("

绘制注释可视化...
")

  # UMAP 注释图
  p_cell_umap <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                         group.by = "cell_type", label = TRUE, repel = TRUE,
                         label.size = 3.5)
  ggsave("results/10_cellAnn_UMAP.pdf", p_cell_umap, width = 12, height = 8)

  # tSNE 注释图
  p_cell_tsne <- DimPlot(biofsci, reduction = "tsne", pt.size = 0.5,
                         group.by = "cell_type", label = TRUE, repel = TRUE,
                         label.size = 3.5)
  ggsave("results/10_cellAnn_tSNE.pdf", p_cell_tsne, width = 12, height = 8)

  # 按分组拆分
  p_group_ann <- DimPlot(biofsci, reduction = "umap", pt.size = 0.5,
                         group.by = "cell_type", label = FALSE,
                         split.by = "group")
  ggsave("results/10_cellAnn_UMAP_split_by_group.pdf", p_group_ann,
         width = 16, height = 7)

  # 细胞类型比例 (Before vs After)
  cell_prop <- as.data.frame(prop.table(
    table(biofsci$cell_type, biofsci$group), margin = 2
  ))
  colnames(cell_prop) <- c("CellType", "Group", "Proportion")

  p_prop <- ggplot(cell_prop, aes(x = Group, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "分组", y = "细胞比例", title = "各分组细胞类型比例") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  ggsave("results/10_cell_type_proportion.pdf", p_prop, width = 10, height = 7)

  # 每个患者各分组细胞类型比例
  patient_prop <- as.data.frame(prop.table(
    table(biofsci$cell_type, biofsci$patient, biofsci$group), margin = c(2, 3)
  ))
  colnames(patient_prop) <- c("CellType", "Patient", "Group", "Proportion")

  p_patient_prop <- ggplot(patient_prop,
                           aes(x = Group, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Patient, nrow = 1) +
    theme_minimal() +
    labs(x = "分组", y = "细胞比例", title = "各患者细胞类型比例") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("results/10_cell_type_proportion_by_patient.pdf", p_patient_prop,
         width = 16, height = 6)

  cat("

--- 细胞类型统计 ---
")
  print(table(biofsci$cell_type))
  cat(sprintf("

共注释出 %d 种细胞类型
", length(unique(biofsci$cell_type))))

  # 每聚类对应细胞类型与占比
  cat("

每聚类的细胞类型分布:
")
  for (cl in sort(unique(biofsci$seurat_clusters))) {
    cells_in_cluster <- colnames(biofsci)[biofsci$seurat_clusters == cl]
    if (length(cells_in_cluster) > 0) {
      cl_types <- table(biofsci$cell_type[cells_in_cluster])
      cl_types <- sort(cl_types, decreasing = TRUE)
      top_types <- cl_types[1:min(3, length(cl_types))]
      cat(sprintf("  聚类 %s (%d cells): %s
",
                  cl, length(cells_in_cluster),
                  paste(names(top_types),
                        sprintf("(%d)", top_types), collapse = ", ")))
    }
  }

} else {
  cat("

!!! 警告: 未找到 SingleR 参考数据集 !!!
")
  cat(sprintf("  期望路径: %s
", ref_dir))
  cat("  请确认以下文件存在:
")
  cat(sprintf("    %s
", ref_human_all))
  cat(sprintf("    %s
", ref_hematopoietic))
  cat("
  跳过注释步骤。后续可手动运行注释。
")
}

# ==================== 7. 保存结果 ====================
# 保存 Seurat 对象
saveRDS(biofsci, "biofsci04_clustered.rds")
cat("
已保存: biofsci04_clustered.rds
")

# 保存聚类标记基因
if (exists("sig_markers")) {
  save(biofsci, sig_markers, file = "Seurat04.Rdata")
  cat("已保存: Seurat04.Rdata
")
}

# ==================== 8. 输出注释审核提示 ====================
cat("
========================================
===== 细胞类型注释完成，请审核！ =====
========================================

自动注释结果已保存到:
  - results/09_clusterAnn.csv  (聚类→细胞类型映射)
  - results/09_cellAnn.csv     (每个细胞的注释)
  - results/10_cellAnn_UMAP.pdf (UMAP注释可视化)
  - results/10_cell_type_proportion.pdf (细胞比例)

请检查以下内容:
  1. 每个聚类的注释是否合理 (对照标记基因)
  2. 考虑是否需要合并相似的细胞类型
  3. 考虑是否需要手动调整特定聚类的注释

如需调整某聚类的注释，请编辑 results/09_clusterAnn.csv
中的label列，然后重新运行后续分析。

========== Step4: 聚类与注释 完毕 ==========
")
