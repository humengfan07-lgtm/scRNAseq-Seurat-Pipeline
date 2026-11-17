# ITP 治疗前后配对分析

每批使用一个独立的运行目录。原始矩阵可以留在原来的位置。若将新患者加入原队列，在新目录里用包含原有和新增样本的样本表，从 Step1 重跑。Step5 会沿用当前目录的同名结果，适合继续补跑同一批数据。

## 哪些自动，哪些手动

| 步骤 | 脚本自动做 | 你需要确认或填写 |
| --- | --- | --- |
| Step0 扫描样本路径 | 找到 filtered 矩阵、提取样本编号、填好绝对路径。有样本记录时，按编号导入患者和分组 | 指定原始数据总目录。没有样本记录时，按实际采样填写患者标识和 Before/After |
| Step1 检查分组 | 检查样本编号、治疗前后配对和矩阵文件，保存 `samples.csv` | 核对表中的配对和时间点是否符合采样记录 |
| Step2 质控 | 读入矩阵、按脚本阈值过滤、输出质控图表 | 看各样本的 QC 图和保留细胞数，必要时调整阈值 |
| Step3 整合与降维 | 标准化、PCA、Harmony、UMAP 和 t-SNE | 检查 PCA 和样本在降维图中的分布 |
| Step4 / 4b 聚类与候选注释 | 生成聚类、marker 和可用的参考注释。两个 Step4b 可选 | 对照 marker 和参考标签审核本批各簇 |
| Step4c / 4d 最终注释 | 生成本批聚类审核表，检查完整性并写回 Seurat | 在 `cluster_annotations.csv` 填每个聚类的 `cell_type`。目标类型用 `Monocytes` 和 `NK` |
| Step5 差异分析 | 逐患者 Wilcoxon 和火山图，配对 pseudobulk DESeq2 | 检查输出和同一细胞类型的可用配对情况 |
| Step6 富集分析 | 用 Step5 配对表做 KEGG ORA 和 GSEA，保存图表 | 按研究问题解读通路和方向 |

## 1. 自动填样本路径

RStudio Console 中把下面的路径改成你本机的实际位置。工作目录决定结果保存到哪里。

```r
code <- "/你的代码目录/scRNAseq-Seurat-Pipeline/itp"
setwd("/你的数据目录/ITP_新批次")

source(file.path(code, "Step0_扫描样本路径.R"))
prepare_itp_samples(raw_root = "/你的原始数据总目录")
```

脚本扫描 `filtered_cell_gene_matrix` 和 `filtered_feature_bc_matrix`，支持 `S01/filtered_feature_bc_matrix` 和 `S01/outs/filtered_feature_bc_matrix` 两种结构。自动填 ID 和路径，打开生成的 `samples_input.csv` 补齐患者和治疗前后。

如果已有本地 `sample_metadata.csv`，包含 `sample_id`、`patient`、`group` 三列，可直接用下面的调用代替上面的 `prepare_itp_samples()`。编号需与样本文件夹名一致。已有的 `samples_input.csv` 会保留。

```r
prepare_itp_samples(
  raw_root = "/你的原始数据总目录",
  metadata_file = "/你的本地记录目录/sample_metadata.csv"
)
```

也可以直接复制[虚拟样本表](samples_input.example.csv)手工填表。原始矩阵、真实样本表、注释表、RDS 和结果都保存在本地。需要 SingleR 参考注释时，把 `ref_Human_all.RData` 和 `ref_Hematopoietic.RData` 放在批次目录的 `references/`，或设置 `Sys.setenv(ITP_REF_DIR = "/本地参考目录")`。没有参考文件时，Step4 仍可生成 marker。

## 2. 在同一个批次目录逐步运行

```r
source(file.path(code, "Step1_数据整理与分组.R"))
source(file.path(code, "Step2_数据读取与质控过滤.R"))
source(file.path(code, "Step3_数据整合与降维.R"))
source(file.path(code, "Step4_聚类与注释.R"))
source(file.path(code, "Step4c_校正注释.R"))
```

对照 `results/08_clusterMarkers.csv`、marker 表和参考注释填好 `cluster_annotations.csv`，再继续。每批的聚类编号都要重新审核。

```r
source(file.path(code, "Step4d_最终注释.R"))
source(file.path(code, "Step5_NK与单核细胞配对差异分析.R"))
source(file.path(code, "Step6_NK与单核细胞富集分析.R"))
```

Step5 的 DESeq2 设计式是 `~ patient + group`。差异与富集结果按细胞类型保存在 `差异富集 result/`。CSV 在 `05_differential/`，PDF 在 `06_visualization/`。

## 主要 R 依赖

- CRAN：Seurat、dplyr、magrittr、ggplot2、patchwork、harmony、clustree、data.table、circlize、tidyr、msigdbr、ggridges。
- Bioconductor：SingleR、DESeq2、ComplexHeatmap、clusterProfiler、enrichplot、org.Hs.eg.db、fgsea、DOSE。
- Step4 还依赖 presto。
