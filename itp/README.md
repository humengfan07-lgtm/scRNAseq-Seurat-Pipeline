# ITP 治疗前后配对分析

每次分析使用一个独立的批次目录。同一批患者的治疗前后样本一起进入该目录。如果将新患者加入原队列，用包含原有和新增样本的样本表，从 Step1 在新目录重跑。Step5 按当前目录查找已有结果，同名 CSV 已存在时会沿用，适合继续补跑当前批次。

## 1. 准备本地样本表

在空的批次目录中，将 [样本表模板](samples_input.example.csv) 复制为 \`samples_input.csv\`，一行对应一个测序样本。

| 列 | 填写内容 |
| --- | --- |
| \`sample_id\` | 每个样本一个唯一编号 |
| \`patient\` | 同一患者的治疗前后填写相同标识。逐患者结果会使用该标识命名 |
| \`group\` | 按实验记录填写 \`Before\` 或 \`After\`，每位患者各一个 |
| \`path\` | 包含 \`matrix.mtx.gz\`、\`features.tsv.gz\`、\`barcodes.tsv.gz\` 的目录。可填绝对路径，也可填相对于批次目录的路径 |

原始矩阵、实际样本表、注释表、RDS 和结果始终保存在本地。建议使用患者研究编号。若需要运行 SingleR，把 \`ref_Human_all.RData\` 和 \`ref_Hematopoietic.RData\` 放在批次目录的 \`references/\`，或用 \`Sys.setenv(ITP_REF_DIR = "/你的本地参考目录")\` 指向现有参考文件。没有参考文件时，Step4 仍可生成 marker 结果，人工注释时可以据此审核。

## 2. 在 RStudio 逐步运行

将以下两个路径换成代码和这次批次目录的实际位置。脚本通过 \`source()\` 读取，输出写入 \`setwd()\` 指定的批次目录。

\`\`\`r
code <- "/你的代码目录/scRNAseq-Seurat-Pipeline/itp"
setwd("/你的数据目录/ITP_新批次")

source(file.path(code, "Step1_数据整理与分组.R"))      # 确认样本和配对
source(file.path(code, "Step2_数据读取与质控过滤.R"))  # Read10X 和 QC
source(file.path(code, "Step3_数据整合与降维.R"))      # 标准化、PCA、Harmony、UMAP/t-SNE
source(file.path(code, "Step4_聚类与注释.R"))          # 聚类、marker、参考注释
source(file.path(code, "Step4c_校正注释.R"))           # 生成聚类审核表和 marker 表
\`\`\`

可选的两个 Step4b 脚本提供更细的参考注释。先看 \`results/08_clusterMarkers.csv\` 和可用的参考注释结果，再在本批生成的 \`cluster_annotations.csv\` 填写每个聚类的最终 \`cell_type\`。新批次的聚类编号需重新核对。要分析的两类细胞应使用标签 \`Monocytes\` 和 \`NK\`。

填好注释后，再在同一批次目录运行下面三步。

\`\`\`r
source(file.path(code, "Step4d_最终注释.R"))
source(file.path(code, "Step5_NK与单核细胞配对差异分析.R"))
source(file.path(code, "Step6_NK与单核细胞富集分析.R"))
\`\`\`

Step5 对每位患者分别做细胞层面的 Wilcoxon 比较，并用原始计数按患者和时间点汇总，拟合 \`~ patient + group\` 的 DESeq2 模型。Step6 从本批生成的配对差异表做 KEGG ORA 和 GSEA。结果分细胞类型存放于 \`差异富集 result/\`，CSV 在 \`05_differential/\`，PDF 在 \`06_visualization/\`。

## 主要 R 依赖

- CRAN：Seurat、dplyr、magrittr、ggplot2、patchwork、harmony、clustree、data.table、circlize、tidyr、msigdbr、ggridges。
- Bioconductor：SingleR、DESeq2、ComplexHeatmap、clusterProfiler、enrichplot、org.Hs.eg.db、fgsea、DOSE。
- Step4 还依赖 presto。参考注释和富集数据库使用本机可用版本。
