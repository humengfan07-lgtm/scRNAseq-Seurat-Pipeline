# 导出 Azimuth 人 PBMC 官方 Annotation Markers，保留原始标签和基因列表。
download_azimuth_pbmc_markers <- function(
    output_file = file.path("references", "azimuth_human_pbmc_markers.csv")) {
  if (file.exists(output_file)) stop("导出文件已存在：", output_file)
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop('请先运行 install.packages("xml2")。')
  }

  source_url <- "https://azimuth.hubmapconsortium.org/references/human_pbmc/"
  page <- xml2::read_html(source_url)
  tables <- xml2::xml_find_all(page, "//details//table[contains(@class, 'marker_table')]")
  expected_levels <- c("celltype.l1", "celltype.l2 (default)", "celltype.l3")
  levels <- vapply(tables, function(table) {
    xml2::xml_text(xml2::xml_find_first(table, "./ancestor::details[1]/summary"), trim = TRUE)
  }, character(1))
  if (length(levels) != 3L || !setequal(levels, expected_levels)) {
    stop("Azimuth 官方页面结构可能已更新，请直接查看该页面的 Annotation Markers。")
  }

  results <- Map(function(table, level) {
    headers <- xml2::xml_text(xml2::xml_find_all(table, "./thead//th"), trim = TRUE)
    expected_headers <- c("Label", "Expanded Label", "OBO Ontology ID", "Markers")
    rows <- xml2::xml_find_all(table, ".//tr[td]")
    if (!identical(headers, expected_headers) || !length(rows)) {
      stop("无法读取 Azimuth 官方 marker 表：", level)
    }
    values <- t(vapply(rows, function(row) {
      xml2::xml_text(xml2::xml_find_all(row, "./td"), trim = TRUE)
    }, character(length(expected_headers))))
    data.frame(reference_level = level, cell_label = values[, 1],
               expanded_label = values[, 2], ontology_id = values[, 3],
               markers = values[, 4], source_url = source_url,
               stringsAsFactors = FALSE)
  }, tables, levels)
  markers <- do.call(rbind, results)
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(markers, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  message("已导出 Azimuth 官方 marker 表：", output_file, "（", nrow(markers), " 行）")
  invisible(markers)
}
