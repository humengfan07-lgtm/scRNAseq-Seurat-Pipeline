# 扫描 Cell Ranger 的 filtered 矩阵目录，生成本地样本表。
prepare_itp_samples <- function(raw_root, metadata_file = NULL,
                                output_file = "samples_input.csv") {
  if (file.exists(output_file)) {
    stop("样本表已存在，请使用现有表或指定另一个输出文件：", output_file)
  }
  raw_root <- path.expand(raw_root)
  if (!dir.exists(raw_root)) stop("找不到原始数据目录：", raw_root)
  folders <- list.dirs(normalizePath(raw_root), full.names = TRUE, recursive = TRUE)
  matrix_names <- c("filtered_cell_gene_matrix", "filtered_feature_bc_matrix")
  folders <- folders[basename(folders) %in% matrix_names]
  required_files <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  valid <- vapply(folders, function(path) {
    all(file.exists(file.path(path, required_files)))
  }, logical(1))
  if (any(!valid)) {
    stop("以下 filtered 矩阵目录缺少文件：\n", paste(folders[!valid], collapse = "\n"))
  }
  if (!length(folders)) stop("未找到包含三件套的 filtered 矩阵目录。")

  # 同时支持 sample/filtered_* 和 sample/outs/filtered_*。
  parent <- dirname(folders)
  ids <- basename(ifelse(basename(parent) == "outs", dirname(parent), parent))
  if (anyDuplicated(ids)) {
    repeated <- unique(ids[duplicated(ids)])
    stop("发现同名样本目录，请检查重复矩阵：", paste(repeated, collapse = ", "))
  }
  samples <- data.frame(sample_id = ids, patient = "", group = "", path = folders)
  samples <- samples[order(samples$sample_id), , drop = FALSE]

  if (!is.null(metadata_file)) {
    metadata <- read.csv(path.expand(metadata_file), colClasses = "character",
                         fileEncoding = "UTF-8", stringsAsFactors = FALSE)
    required <- c("sample_id", "patient", "group")
    if (!all(required %in% names(metadata))) {
      stop("样本记录需要包含 sample_id、patient、group 三列。")
    }
    metadata[required] <- lapply(metadata[required], trimws)
    if (!nrow(metadata) || anyNA(metadata[required]) ||
        any(vapply(metadata[required], function(x) any(!nzchar(x)), logical(1)))) {
      stop("样本记录的 sample_id、patient 和 group 不能为空。")
    }
    if (anyDuplicated(metadata$sample_id) || !all(metadata$group %in% c("Before", "After"))) {
      stop("样本记录的 sample_id 必须唯一，group 只能用 Before 或 After。")
    }
    index <- match(metadata$sample_id, samples$sample_id)
    if (anyNA(index)) {
      stop("没有找到这些样本的矩阵目录：", paste(metadata$sample_id[is.na(index)], collapse = ", "))
    }
    samples <- data.frame(metadata[required], path = samples$path[index])
  }

  write.csv(samples, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  message("已生成：", output_file, "（", nrow(samples), " 个样本）")
  if (is.null(metadata_file)) message("请按采样记录补齐 patient 和 group，再运行 Step1。")
  invisible(samples)
}
