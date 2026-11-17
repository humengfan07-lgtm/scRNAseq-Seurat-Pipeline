# 在当前批次目录运行。每行对应一个测序样本，治疗前后来自真实采样记录。
samples <- read.csv("samples_input.csv", colClasses = "character", fileEncoding = "UTF-8")
required <- c("sample_id", "patient", "group", "path")
if (!all(required %in% names(samples))) {
  stop("样本表需要四列：sample_id, patient, group, path。")
}
samples[required] <- lapply(samples[required], trimws)
if (nrow(samples) < 4L || anyNA(samples[required]) ||
    any(vapply(samples[required], function(x) any(!nzchar(x)), logical(1)))) {
  stop("请填完整样本表，至少包含两位患者的治疗前后配对。")
}
if (anyDuplicated(samples$sample_id) || !all(samples$group %in% c("Before", "After"))) {
  stop("sample_id 必须唯一，group 只能用 Before 或 After。")
}
pairs <- table(samples$patient, factor(samples$group, levels = c("Before", "After")))
if (any(pairs != 1L)) stop("每位患者需要恰好一个 Before 和一个 After 样本。")
files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
samples$path <- path.expand(samples$path)
samples$valid <- vapply(samples$path, function(path) {
  all(file.exists(file.path(path, files)))
}, logical(1))
if (any(!samples$valid)) {
  stop("以下样本的 10x 矩阵目录不完整：", paste(samples$sample_id[!samples$valid], collapse = ", "))
}
samples$path <- normalizePath(samples$path)
if (anyDuplicated(samples$path)) stop("每行应指向不同的样本目录。")
samples$status <- "OK"
write.csv(samples, "samples.csv", row.names = FALSE, fileEncoding = "UTF-8")
for (folder in c("QC_plots", "results", "diff_genes")) dir.create(folder, showWarnings = FALSE)
print(pairs)
message("Step1 完成。Step2 将读取本批目录下的 samples.csv。")
