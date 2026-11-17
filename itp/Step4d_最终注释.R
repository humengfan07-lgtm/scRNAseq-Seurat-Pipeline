#######################################################
# ITP 外周血 - 最终手动注释 (严格按照用户审核)
# 根据用户审核的cluster→cell type映射表
#######################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(ComplexHeatmap)
library(circlize)

cat("========================================\n")
cat("ITP 外周血 | 最终手动注释 (v2)\n")
cat("========================================\n\n")

biofsci <- readRDS("biofsci04_clustered.rds")
cat(sprintf("加载数据: %d cells\n", ncol(biofsci)))

# ==================== 1. 手动注释映射表 ====================
# 严格按照用户审核结果

clusters <- as.character(biofsci$seurat_clusters)
cluster_ids <- sort(unique(clusters))
annotation <- read.csv("cluster_annotations.csv", colClasses = "character",
                       stringsAsFactors = FALSE)
if (!all(c("cluster", "cell_type") %in% names(annotation))) {
  stop("cluster_annotations.csv 必须包含 cluster 和 cell_type 列。")
}
annotation$cell_type <- trimws(annotation$cell_type)
if (anyDuplicated(annotation$cluster) || !setequal(annotation$cluster, cluster_ids) ||
    anyNA(annotation$cell_type) || any(!nzchar(annotation$cell_type))) {
  stop("请给本批每个聚类填写一条最终注释，聚类编号须与 Step4 完全匹配。")
}
manual_annotation <- setNames(annotation$cell_type, annotation$cluster)

# ==================== 2. 应用到Seurat对象 ====================
new_labels <- manual_annotation[as.character(biofsci$seurat_clusters)]
names(new_labels) <- colnames(biofsci)
biofsci <- AddMetaData(biofsci, metadata = new_labels, col.name = "cell_type_final")

# 设置因子顺序
cell_type_order <- c(
  "Monocytes", "DC", "Neutrophil", "Erythroid", "Platelets",
  "CD4+ T", "CD8+ T", "NK",
  "B cell", "Proliferating cell"
)
biofsci$cell_type_final <- factor(biofsci$cell_type_final, levels = unique(c(cell_type_order, annotation$cell_type)))

cat("\n\n最终细胞类型分布:\n")
print(table(biofsci$cell_type_final))

# ==================== 3. 可视化 ====================
cat("\n绘制可视化...\n")

# 颜色方案
colors_final <- c(
  "Monocytes"         = "#E64B35",
  "DC"                = "#F39B7F",
  "Neutrophil"        = "#D95F02",
  "Erythroid"         = "#CC0033",
  "Platelets"         = "#E41A1C",
  "CD4+ T"            = "#4DBBD5",
  "CD8+ T"            = "#00A087",
  "NK"                = "#3C5488",
  "B cell"            = "#7E6148",
  "Proliferating cell" = "#B09C85"
)
extra_types <- setdiff(annotation$cell_type, names(colors_final))
if (length(extra_types)) {
  colors_final <- c(colors_final, setNames(grDevices::hcl.colors(length(extra_types), "Dark 3"), extra_types))
}


# UMAP 最终注释
p_umap <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                  group.by = "cell_type_final", label = TRUE, repel = TRUE,
                  cols = colors_final, label.size = 3.5) +
  ggtitle("ITP Peripheral Blood - Manual Annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("results/13_final_annotation_UMAP.pdf", p_umap, width = 13, height = 9)
ggsave("results/13_final_annotation_UMAP.png", p_umap, width = 2600, height = 1800,
       dpi = 200, units = "px")

# UMAP 按 Before/After 拆分
p_split <- DimPlot(biofsci, reduction = "umap", pt.size = 0.3,
                   group.by = "cell_type_final", label = FALSE,
                   split.by = "group", cols = colors_final, ncol = 2)
ggsave("results/13_final_annotation_UMAP_by_group.pdf", p_split, width = 16, height = 7)

# UMAP 按患者拆分
p_patient <- DimPlot(biofsci, reduction = "umap", pt.size = 0.2,
                     group.by = "cell_type_final", label = FALSE,
                     split.by = "patient", cols = colors_final, ncol = 3)
ggsave("results/13_final_annotation_UMAP_by_patient.pdf", p_patient, width = 20, height = 7)

# ==================== 4. 细胞比例 ====================
cell_prop <- as.data.frame(prop.table(
  table(biofsci$cell_type_final, biofsci$group), margin = 2
))
colnames(cell_prop) <- c("CellType", "Group", "Proportion")

p_prop <- ggplot(cell_prop, aes(x = Group, y = Proportion * 100, fill = CellType)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = colors_final) +
  theme_minimal(base_size = 12) +
  labs(x = "", y = "Proportion (%)",
       title = "Cell Type Proportions: Before vs After Treatment",
       fill = "Cell Type") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("results/13_cell_proportion_before_after.pdf", p_prop, width = 11, height = 7)

# 每个患者各组的细胞比例
patient_prop <- as.data.frame(prop.table(
  table(biofsci$cell_type_final, biofsci$patient, biofsci$group),
  margin = c(2, 3)
))
colnames(patient_prop) <- c("CellType", "Patient", "Group", "Proportion")

p_patient_prop <- ggplot(patient_prop,
                         aes(x = Group, y = Proportion * 100, fill = CellType)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = colors_final) +
  facet_wrap(~ Patient, nrow = 1) +
  theme_minimal(base_size = 11) +
  labs(x = "", y = "Proportion (%)",
       title = "Cell Type Proportions by Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("results/13_cell_proportion_by_patient.pdf", p_patient_prop, width = 18, height = 7)

# ==================== 5. Before/After 变化统计 ====================
cat("\n\n========== Cell Type Changes: Before vs After ==========\n")

summary_df <- data.frame(
  CellType = character(),
  Before = integer(),
  After = integer(),
  Before_Pct = numeric(),
  After_Pct = numeric(),
  Change = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

total_before <- sum(biofsci$group == "Before")
total_after  <- sum(biofsci$group == "After")

for (ct in levels(biofsci$cell_type_final)) {
  n_before <- sum(biofsci$cell_type_final == ct & biofsci$group == "Before")
  n_after  <- sum(biofsci$cell_type_final == ct & biofsci$group == "After")
  pct_before <- n_before / total_before * 100
  pct_after  <- n_after  / total_after  * 100
  change <- pct_after - pct_before
  direction <- ifelse(change > 0, "Increase", ifelse(change < 0, "Decrease", "Stable"))

  summary_df <- rbind(summary_df, data.frame(
    CellType = ct,
    Before = n_before,
    After = n_after,
    Before_Pct = round(pct_before, 2),
    After_Pct = round(pct_after, 2),
    Change = round(change, 2),
    Direction = direction,
    stringsAsFactors = FALSE
  ))

  arrow <- ifelse(change > 0, "↑", ifelse(change < 0, "↓", "-"))
  cat(sprintf("  %-20s | Before: %5d (%6.2f%%) | After: %5d (%6.2f%%) | %+6.2f%% %s\n",
              ct, n_before, pct_before, n_after, pct_after, change, arrow))
}

fwrite(summary_df, "results/13_cell_proportion_summary.csv")
cat("\nSaved: results/13_cell_proportion_summary.csv\n")

# 瀑布图
p_change <- ggplot(summary_df, aes(x = reorder(CellType, Change), y = Change)) +
  geom_bar(stat = "identity", aes(fill = Direction), width = 0.7) +
  scale_fill_manual(values = c("Increase" = "#E41A1C", "Decrease" = "#377EB8", "Stable" = "#999999")) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(x = "", y = "Change in Proportion (%)",
       title = "Cell Type Changes: Before -> After Treatment",
       subtitle = "ITP Peripheral Blood") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("results/13_cell_change_waterfall.pdf", p_change, width = 10, height = 6)

# ==================== 6. ITP 关键Marker热图 ====================
cat("\nGenerating key marker heatmap...\n")

itp_key_markers <- c(
  # T cells
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B",
  # Treg / exhaustion
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "PDCD1",
  # NK
  "NKG7", "GNLY", "KLRF1", "PRF1",
  # B cells
  "MS4A1", "CD79A", "CD79B",
  # Plasma / Plasma-like
  "SDC1", "MZB1", "IGHG1", "JCHAIN",
  # Monocyte
  "CD14", "FCGR3A", "S100A8", "S100A9", "LYZ",
  # DC
  "FCER1A", "CLEC10A", "CLEC4C",
  # Neutrophil
  "FCGR3B", "CSF3R", "CXCR2", "MMP8", "CAMP", "LTF",
  # Platelet/MK
  "PPBP", "PF4", "ITGA2B", "ITGB3", "GP9",
  # Erythroid
  "HBB", "HBA1", "HBA2", "GYPA",
  # Proliferation
  "MKI67", "TOP2A", "PCNA"
)

avail <- intersect(itp_key_markers, rownames(biofsci))
cat(sprintf("Key markers available: %d/%d\n", length(avail), length(itp_key_markers)))

avg_expr <- AverageExpression(biofsci, features = avail,
                               assays = "RNA", group.by = "cell_type_final")
expr_mat <- as.matrix(avg_expr$RNA)
expr_scaled <- t(scale(t(expr_mat)))
expr_scaled[!is.finite(expr_scaled)] <- 0
expr_scaled <- pmax(pmin(expr_scaled, 2), -2)

ht <- Heatmap(expr_scaled,
              name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_rot = 45,
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(direction = "horizontal"),
              column_title = "ITP Key Marker Expression by Cell Type",
              width = unit(12, "cm"),
              height = unit(16, "cm"))

pdf("results/13_key_marker_heatmap.pdf", width = 12, height = 14)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

png("results/13_key_marker_heatmap.png", width = 1400, height = 1600, res = 150)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

cat("Saved: results/13_key_marker_heatmap.pdf/png\n")

# ==================== 7. 保存 ====================
biofsci$cell_type <- biofsci$cell_type_final

# 清理旧的meta列
for (col in c("cell_type_fine", "cell_type_fine_clean", "cell_type_corrected")) {
  if (col %in% colnames(biofsci@meta.data)) {
    biofsci[[col]] <- NULL
  }
}

saveRDS(biofsci, "biofsci04_clustered.rds")
cat("\nUpdated: biofsci04_clustered.rds\n")

annot_table <- data.frame(
  Cluster = cluster_ids,
  CellType = unname(manual_annotation[cluster_ids]),
  NCells = as.integer(table(clusters)[cluster_ids]),
  stringsAsFactors = FALSE
)
fwrite(annot_table, "results/13_final_annotation_table.csv")
cat("Saved: results/13_final_annotation_table.csv\n")

cat("\n========== Final annotation complete ==========\n")
cat("10 cell types, 38,229 cells\n")