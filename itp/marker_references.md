# 外周血单细胞注释的官方 marker 参考

## Azimuth 人 PBMC 图谱

在官方页面展开 **Annotation Markers**。表格按 celltype.l1、celltype.l2 和 celltype.l3 列出不同细分程度的细胞标签、基因组合及对应的细胞本体术语。审核 Step4 时，可先对照 l1 或 l2。需要细分亚群时，再看 l3 的相应条目。 [Azimuth 人 PBMC 官方 marker 表](https://azimuth.hubmapconsortium.org/references/human_pbmc/) · [Hao 等，Cell，2021](https://doi.org/10.1016/j.cell.2021.04.048)

需要在本地搜索或筛选时，可在 RStudio 中运行以下代码。这里的 `code` 变量指向仓库的 `itp` 目录，与运行说明一致。

```r
source(file.path(code, "download_azimuth_pbmc_markers.R"))
download_azimuth_pbmc_markers()
```

CSV 写入当前批次工作目录下的 `references/azimuth_human_pbmc_markers.csv`。它保留官网三个层级的所有行及原始标签和 marker 列，附有来源链接。导出已有文件时会提示换文件名。

## CellMarker 2.0

查询其他细胞类型或组织时，可按人类、组织和细胞类型查看原始数据库条目。给本项目的外周血细胞做注释时，优先检查记录是否来自血液或 PBMC，核对条目的证据来源。 [CellMarker 2.0 数据库](https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/) · [数据库论文](https://doi.org/10.1093/nar/gkac947)

## 接到 Step4 的人工审核

1. 从本批的 `results/08_clusterMarkers.csv` 按聚类查看候选基因。
2. 查官网同一细胞类型的多基因组合，并用本批的表达图核对。
3. 审核后，把最终标签填入本地 `cluster_annotations.csv`，再运行 Step4d。
