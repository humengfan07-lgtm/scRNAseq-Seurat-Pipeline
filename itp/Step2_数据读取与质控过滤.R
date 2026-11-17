#######################################################
# ITP 单细胞测序分析 - Step2: 数据读取与质控过滤
# 治疗前(Before) vs 治疗后(After) 配对设计
# Seurat v5 兼容 - 使用 AddMetaData 添加元数据
#######################################################

library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)

set.seed(2024)

cat("========================================\n")
cat("ITP | Step2: 数据读取与质控过滤\n")
cat("========================================\n\n")

# ==================== 1. 加载分组信息 ====================
samples <- read.csv("samples.csv", stringsAsFactors = FALSE,
                    colClasses = c(sample_id = "character", patient = "character"))
cat(sprintf("分组表: %d 个样本\n", nrow(samples)))
if (!nrow(samples) || !all(samples$valid)) stop("请先运行 Step1 验证本批样本。")
cat(sprintf("有效样本: %d 个\n", nrow(samples)))

for (i in seq_len(nrow(samples))) {
  cat(sprintf("  %s: %s (%s)\n",
              samples[["sample_id"]][i],
              samples[["patient"]][i],
              samples[["group"]][i]))
}

# ==================== 2. 批量读取10X数据 ====================
cat("\n读取10X数据...\n")
dirs <- setNames(samples[["path"]], samples[["sample_id"]])
counts <- Read10X(data.dir = dirs)
cat(sprintf("  矩阵维度: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# ==================== 3. 创建Seurat对象 ====================
biofsci <- CreateSeuratObject(
  counts = counts,
  project = "ITP",
  assay = "RNA",
  min.cells = 3,
  min.features = 200
)
cat(sprintf("Seurat对象: %d cells, %d genes\n", ncol(biofsci), nrow(biofsci)))

# ==================== 4. 添加元数据 (Seurat v5 兼容: AddMetaData) ====================
# 从细胞名提取样本编号: "ITP_条形码序列"
cell_sample <- sub("_[^_]+$", "", colnames(biofsci))
if (!all(cell_sample %in% samples$sample_id)) stop("细胞条形码无法对应到样本。")

biofsci <- AddMetaData(biofsci, metadata = cell_sample, col.name = "sample_id")

sample_to_patient <- setNames(samples[["patient"]], samples[["sample_id"]])
sample_to_group   <- setNames(samples[["group"]],   samples[["sample_id"]])

pat <- sample_to_patient[biofsci$sample_id]
names(pat) <- colnames(biofsci)
biofsci <- AddMetaData(biofsci, metadata = pat, col.name = "patient")

grp <- sample_to_group[biofsci$sample_id]
names(grp) <- colnames(biofsci)
grp <- factor(grp, levels = c("Before", "After"))
biofsci <- AddMetaData(biofsci, metadata = grp, col.name = "group")

cat("\n元数据统计:\n")
cat(sprintf("  患者数: %d\n", length(unique(biofsci$patient))))
print(table(biofsci$patient, biofsci$group))
cat(sprintf("\n  总细胞数: %d\n", ncol(biofsci)))

# ==================== 5. QC指标计算 ====================
biofsci[["percent.mt"]] <- PercentageFeatureSet(biofsci, pattern = "^MT-")
biofsci[["percent.rb"]] <- PercentageFeatureSet(biofsci, pattern = "^RP[SL]")

cat("\n--- 过滤前统计 ---\n")
cat(sprintf("  nFeature_RNA 中位数: %d\n", round(median(biofsci$nFeature_RNA))))
cat(sprintf("  nCount_RNA 中位数: %d\n", round(median(biofsci$nCount_RNA))))
cat(sprintf("  percent.mt 中位数: %.2f%%\n", median(biofsci$percent.mt)))
cat(sprintf("  percent.rb 中位数: %.2f%%\n", median(biofsci$percent.rb)))

# ==================== 6. 过滤前QC可视化 ====================
dir.create("QC_plots", showWarnings = FALSE)

p1 <- FeatureScatter(biofsci, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "sample_id", pt.size = 0.3)
ggsave("QC_plots/01_nCount_vs_nFeature_before.pdf", p1, width = 12, height = 8)

p2 <- VlnPlot(biofsci, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
              layer = "counts", ncol = 4, pt.size = 0, group.by = "sample_id")
ggsave("QC_plots/02_QC_before_by_sample.pdf", p2, width = 20, height = 7)

p3 <- VlnPlot(biofsci, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
              layer = "counts", ncol = 4, pt.size = 0, group.by = "group")
ggsave("QC_plots/03_QC_before_by_group.pdf", p3, width = 16, height = 7)

# ==================== 7. 细胞过滤 ====================
cat("\n=== 细胞过滤 ===\n")

count_upper <- quantile(biofsci$nCount_RNA, 0.98)
feat_upper  <- quantile(biofsci$nFeature_RNA, 0.98)
mt_upper    <- 20
rb_upper    <- 30

cat(sprintf("  nCount_RNA 上限 (98%%): %.0f\n", count_upper))
cat(sprintf("  nFeature_RNA 上限 (98%%): %.0f\n", feat_upper))
cat(sprintf("  percent.mt 上限: %.0f%%\n", mt_upper))
cat(sprintf("  percent.rb 上限: %.0f%%\n", rb_upper))
cat(sprintf("  nCount_RNA 下限: 1000\n"))
cat(sprintf("  nFeature_RNA 下限: 200\n"))

mask <- (biofsci$nCount_RNA >= 1000) & (biofsci$nCount_RNA <= count_upper) &
        (biofsci$nFeature_RNA >= 200) & (biofsci$nFeature_RNA <= feat_upper) &
        (biofsci$percent.mt <= mt_upper) &
        (biofsci$percent.rb <= rb_upper)

biofsci_filtered <- biofsci[, mask]

cat(sprintf("\n  过滤后细胞数: %d (保留 %.1f%%)\n",
            ncol(biofsci_filtered),
            ncol(biofsci_filtered) / ncol(biofsci) * 100))

# 各样本过滤统计
filter_stats <- data.frame(
  sample  = names(table(biofsci$sample_id)),
  before  = as.vector(table(biofsci$sample_id)),
  after   = as.vector(table(factor(biofsci_filtered$sample_id, levels = names(table(biofsci$sample_id))))),
  stringsAsFactors = FALSE
)
filter_stats$retained <- round(filter_stats$after / filter_stats$before * 100, 1)
cat("\n各样本过滤统计:\n")
print(filter_stats)

# ==================== 8. 过滤后QC可视化 ====================
p4 <- VlnPlot(biofsci_filtered,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
              layer = "counts", ncol = 4, pt.size = 0.1, group.by = "sample_id")
ggsave("QC_plots/04_QC_after_by_sample.pdf", p4, width = 20, height = 7)

p5 <- VlnPlot(biofsci_filtered,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
              layer = "counts", ncol = 4, pt.size = 0.1, group.by = "group")
ggsave("QC_plots/05_QC_after_by_group.pdf", p5, width = 16, height = 7)

p6 <- FeatureScatter(biofsci_filtered,
                     feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                     group.by = "sample_id", pt.size = 0.3)
ggsave("QC_plots/06_nCount_vs_nFeature_after.pdf", p6, width = 12, height = 8)

p7 <- ggplot(filter_stats, aes(x = sample, y = after, fill = sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "Cells", title = "Filtered cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave("QC_plots/07_cell_counts_by_sample.pdf", p7, width = 10, height = 6)

p8 <- ggplot(as.data.frame(table(biofsci_filtered$patient)),
             aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() + labs(x = "", y = "Cells", title = "Cells per patient") +
  theme(legend.position = "none")
# 用 macOS 的 Quartz PDF 设备和苹方字体保留中文横轴标签。
ggsave("QC_plots/08_cell_counts_by_patient.pdf",
       p8 + theme(axis.text.x = element_text(family = "PingFang SC")),
       width = 8, height = 5, device = grDevices::quartz, type = "pdf")

p9 <- ggplot(as.data.frame(table(biofsci_filtered$group)),
             aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() + labs(x = "", y = "Cells", title = "Before vs After") +
  theme(legend.position = "none")
ggsave("QC_plots/09_cell_counts_by_group.pdf", p9, width = 6, height = 5)

# ==================== 9. 保存 ====================
saveRDS(biofsci_filtered, "biofsci02_filtered.rds")
cat("\n已保存: biofsci02_filtered.rds\n")

# ==================== 10. 质控报告 ====================
qc_report <- data.frame(
  指标 = c("过滤前细胞总数","过滤后细胞总数","保留比例",
           "nCount_RNA中位数","nFeature_RNA中位数",
           "percent.mt中位数","percent.rb中位数",
           "Before组细胞数","After组细胞数",
           paste0(sort(unique(samples$patient)), "细胞数")),
  数值 = c(
    ncol(biofsci), ncol(biofsci_filtered),
    sprintf("%.1f%%", ncol(biofsci_filtered)/ncol(biofsci)*100),
    round(median(biofsci_filtered$nCount_RNA)),
    round(median(biofsci_filtered$nFeature_RNA)),
    sprintf("%.2f%%", median(biofsci_filtered$percent.mt)),
    sprintf("%.2f%%", median(biofsci_filtered$percent.rb)),
    sum(biofsci_filtered$group == "Before"),
    sum(biofsci_filtered$group == "After"),
    as.integer(table(factor(biofsci_filtered$patient, levels = sort(unique(samples$patient)))))
  ),
  stringsAsFactors = FALSE
)
write.csv(qc_report, "QC_report.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("\n========== 质控报告 ==========\n")
print(qc_report)
cat("\n========== Step2 完毕 ==========\n")
