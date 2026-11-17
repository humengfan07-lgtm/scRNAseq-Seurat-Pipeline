# 输入：Step4d 保存的最终注释对象。
# 输出：各细胞类型的配对 DESeq2 表、逐患者 Wilcoxon 表和火山图。
project_dir <- getwd()
input_file <- file.path(project_dir, "biofsci04_clustered.rds")
output_root <- file.path(project_dir, "差异富集 result")
cell_folders <- c(Monocytes = "单核细胞", NK = "NK细胞")

make_pseudobulk <- function(counts, metadata, cell_type) {
  # 每种细胞类型内，按实际测序样本分组，同一样本的细胞原始 UMI 逐基因相加。
  cells <- metadata[which(metadata$cell_type_final == cell_type), , drop = FALSE]
  if (anyNA(cells[, c("sample_id", "patient", "group")])) {
    stop(cell_type, " 的样本、患者或时间点信息有缺失。")
  }
  sample_info <- unique(cells[, c("sample_id", "patient", "group")])
  if (anyDuplicated(sample_info$sample_id) || !all(sample_info$group %in% c("Before", "After"))) {
    stop(cell_type, " 的样本编号与患者或时间点的对应关系有误。")
  }
  sample_info$patient <- factor(sample_info$patient)
  sample_info$group <- factor(sample_info$group, levels = c("Before", "After"))
  pairs <- table(sample_info$patient, sample_info$group)
  if (nrow(pairs) < 2L || any(pairs != 1L)) {
    stop(cell_type, " 需要至少两位完整配对患者，每位各有一个 Before 和一个 After 样本。")
  }

  sample_cells <- split(rownames(cells), as.character(cells$sample_id))
  rownames(sample_info) <- as.character(sample_info$sample_id)
  sample_info <- sample_info[names(sample_cells), , drop = FALSE]
  sample_info$n_cells <- lengths(sample_cells)
  pb_counts <- vapply(
    sample_cells,
    function(barcodes) Matrix::rowSums(counts[, barcodes, drop = FALSE]),
    numeric(nrow(counts))
  )
  rownames(pb_counts) <- rownames(counts)
  pb_counts <- pb_counts[rowSums(pb_counts) > 0, , drop = FALSE]
  list(counts = pb_counts, sample_info = sample_info)
}

save_patient_volcano <- function(de, title, output_file) {
  de <- de[is.finite(de$avg_log2FC) & is.finite(de$p_val_adj), , drop = FALSE]
  de$significance <- "未达阈值"
  passes <- de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.25
  de$significance[passes & de$avg_log2FC > 0] <- "治疗后升高"
  de$significance[passes & de$avg_log2FC < 0] <- "治疗后降低"
  de$minus_log10_p <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))

  plot <- ggplot2::ggplot(de, ggplot2::aes(avg_log2FC, minus_log10_p, color = significance)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = c(
      "未达阈值" = "grey70", "治疗后升高" = "#B63A3A", "治疗后降低" = "#3366A6"
    )) +
    ggplot2::labs(title = title, x = "log2 fold change (After vs Before)",
                  y = "-log10(adjusted p-value)", color = NULL) +
    ggplot2::theme_classic(base_size = 12)
  pdf_device <- function(filename, width, height, bg) {
    if (capabilities("aqua")) {
      grDevices::quartz(type = "pdf", file = filename, family = "PingFang SC",
                        width = width, height = height, bg = bg)
    } else {
      grDevices::pdf(file = filename, family = "GB1", width = width, height = height, bg = bg)
    }
  }
  ggplot2::ggsave(output_file, plot, device = pdf_device, width = 7.5, height = 5)
  if (!file.exists(output_file)) stop("火山图未写入：", output_file)
}

run_patient_de <- function(rna, metadata, cell_type, patients, patient_file_labels,
                           table_dir, plot_dir) {
  # 沿用交付脚本的 Wilcoxon 和 min.pct 设置，每次只选同一患者的一种细胞。
  cells <- metadata[which(metadata$cell_type_final == cell_type), , drop = FALSE]
  for (patient in patients) {
    patient_cells <- cells[which(cells$patient == patient), , drop = FALSE]
    before <- rownames(patient_cells)[which(patient_cells$group == "Before")]
    after <- rownames(patient_cells)[which(patient_cells$group == "After")]
    if (length(before) + length(after) < 10L || min(length(before), length(after)) < 3L) {
      message("跳过 ", cell_type, " / ", patient, "：Before ", length(before),
              " 个细胞，After ", length(after), " 个细胞。")
      next
    }

    label <- paste(patient_file_labels[[patient]], cell_type, sep = "_")
    table_file <- file.path(table_dir, paste0("DE_", label, ".csv"))
    plot_file <- file.path(plot_dir, paste0("Volcano_", label, ".pdf"))
    if (file.exists(table_file) && file.exists(plot_file)) next

    if (file.exists(table_file)) {
      de <- read.csv(table_file, stringsAsFactors = FALSE)
    } else {
      de <- Seurat::FindMarkers(
        rna, cells.1 = after, cells.2 = before, slot = "data",
        test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
        densify = FALSE, verbose = FALSE
      )
      if (nrow(de) == 0L) {
        message("跳过 ", cell_type, " / ", patient, "：没有达到最低检测比例的基因。")
        next
      }
      de$gene <- rownames(de)
      de$patient <- patient
      write.csv(de, table_file, row.names = FALSE)
      message("已生成逐患者差异表：", basename(table_file))
    }
    if (!file.exists(plot_file)) {
      title <- paste(patient, cell_type, "（治疗后 vs 治疗前）")
      save_patient_volcano(de, title, plot_file)
      message("已生成火山图：", basename(plot_file))
    }
  }
}

run_step5 <- function(input = input_file, output = output_root) {
  table_dirs <- file.path(output, cell_folders, "05_differential")
  plot_dirs <- file.path(output, cell_folders, "06_visualization")
  paired_files <- file.path(table_dirs, "DE_pseudobulk_DESeq2.csv")
  pending <- !file.exists(paired_files)
  if (any(pending) && !requireNamespace("DESeq2", quietly = TRUE)) {
    stop('请先在 Console 运行 BiocManager::install("DESeq2", update = FALSE, ask = FALSE)。')
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("请先安装 Seurat。")
  biofsci <- readRDS(input)
  rna <- biofsci[["RNA"]]
  metadata <- biofsci[[]]
  rm(biofsci)
  required <- c("cell_type_final", "patient", "group", "sample_id")
  if (!all(required %in% names(metadata))) {
    stop("注释对象缺少这些列：", paste(setdiff(required, names(metadata)), collapse = ", "))
  }
  selected_cells <- metadata[metadata$cell_type_final %in% names(cell_folders), required]
  if (anyNA(selected_cells$patient) || any(!nzchar(trimws(selected_cells$patient))) ||
      !all(selected_cells$group %in% c("Before", "After"))) {
    stop("目标细胞的 patient 必须有效，group 必须是 Before 或 After。")
  }
  patients <- sort(unique(as.character(selected_cells$patient)))
  # 保留元数据中的患者标识，只替换不适合用于文件名的字符。
  labels <- gsub(r"([^\p{L}\p{N}_.-]+)", "_", patients, perl = TRUE)
  if (anyDuplicated(tolower(labels))) {
    stop("患者标识转换成文件名后发生冲突，请检查 patient 列。")
  }
  patient_file_labels <- setNames(labels, patients)

  # Step4d 对象已有合并的 counts/data；新数据若拆成多个同类 layer，先合并。
  for (layer in c("counts", "data")) {
    layer_names <- SeuratObject::Layers(rna, search = paste0("^", layer, "($|[.])"))
    if (length(layer_names) > 1L || (length(layer_names) == 1L && layer_names != layer)) {
      rna <- SeuratObject::JoinLayers(rna, layers = layer, new = layer)
    }
  }
  if (!"data" %in% SeuratObject::Layers(rna)) {
    stop("RNA assay 缺少归一化 data layer，请先完成 NormalizeData，再运行逐患者差异分析。")
  }
  for (i in seq_along(cell_folders)) {
    cell_type <- names(cell_folders)[i]
    dir.create(table_dirs[i], recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dirs[i], recursive = TRUE, showWarnings = FALSE)
    run_patient_de(rna, metadata, cell_type, patients, patient_file_labels,
                   table_dirs[i], plot_dirs[i])

    if (!pending[i]) {
      message("已存在配对差异表，沿用：", paired_files[i])
      next
    }
    counts <- SeuratObject::LayerData(rna, layer = "counts")
    pb <- make_pseudobulk(counts, metadata, cell_type)
    message("正在分析 ", cell_type, "，完整配对患者数：", nlevels(pb$sample_info$patient))
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = pb$counts, colData = pb$sample_info, design = ~ patient + group
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    # 保留全表，供 GSEA 排序。正的 log2FoldChange 表示治疗后表达较高。
    result <- as.data.frame(DESeq2::results(dds, contrast = c("group", "After", "Before")))
    result$gene <- rownames(result)
    write.csv(result, paired_files[i], row.names = FALSE)
    message("已生成配对差异表：", paired_files[i], "（", nrow(result), " 个基因）")
  }
}

run_step5()
