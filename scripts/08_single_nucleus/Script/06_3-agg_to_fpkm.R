# 加载必要的包
# dplyr: 数据操作
# Seurat: 单细胞数据分析
# patchwork: 用于将多个图形组合成一个布局
# Biobase: 生物信息学分析工具
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
library(tidyverse)

# 读取注释好的 h5ad 文件转换后的数据
# 读取细胞信息（cellinfo）、基因信息（geneinfo）和表达矩阵（counts）
cellinfo <- read.csv("./02_processData/zz_01-cellinfo.csv")  # 细胞元数据
geneinfo <- read.csv("./02_processData/zz_02-geneinfo.csv", row.names = 1)  # 基因元数据
counts <- Matrix::readMM(file = "./02_processData/zz_03-sparse_matrix_counts.mtx")  # 稀疏表达矩阵

#------------- Step 1. 读取 barcode -------------
cell00dpi1 <- read_tsv("00dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell00dpi2 <- read_tsv("00dpi_2_barcodes.tsv", col_names = FALSE)$X1
cell02dpi1 <- read_tsv("02dpi_1_barcodes.tsv", col_names = FALSE)$X1
cell02dpi2 <- read_tsv("02dpi_2_barcodes.tsv", col_names = FALSE)$X1

# 合并为两个集合
set00dpi <- unique(c(cell00dpi1, cell00dpi2))
set02dpi <- unique(c(cell02dpi1, cell02dpi2))

#------------- Step 2. 清洗 cellinfo$X（去掉第二个 "-" 及后内容）-------------
# 原格式如：AAACCCACAACAGCCC-1-0
cellinfo$X_clean <- sub("^([^-]+-[^-]+).*", "\\1", cellinfo$X)

#------------- Step 3. 判断属于哪一类并加后缀 -------------
cellinfo$X_final <- case_when(
  cellinfo$X_clean %in% set00dpi ~ paste0(cellinfo$X_clean, "-0"),
  cellinfo$X_clean %in% set02dpi ~ paste0(cellinfo$X_clean, "-1"),
  TRUE ~ cellinfo$X_clean    # 不在任何列表就保持原样（一般不会出现）
)

#------------- Step 4. 处理重复并设置行名 -------------
cellinfo$X_final_unique <- make.unique(cellinfo$X_final, sep = "_")
rownames(cellinfo) <- cellinfo$X_final_unique
cellinfo <- cellinfo %>% select(-X, -X_clean, -X_final, -X_final_unique)

# 为稀疏矩阵添加行名（基因）和列名（细胞）
rownames(counts) <- rownames(geneinfo)
colnames(counts) <- rownames(cellinfo)

# 创建 Seurat 对象
# 使用表达矩阵 counts 和细胞元数据 cellinfo
All.merge <- CreateSeuratObject(counts = counts, meta.data = cellinfo)

# 将 Seurat 对象保存为 pbmc（更常用的命名）
pbmc <- All.merge

### 将对象拆分为不同时期处理，比较不同时期的表达差异
# 假设 pbmc 是一个 Seurat 对象，提取元数据
metadata <- pbmc@meta.data

# 提取细胞名称，并获取最后一位数字
cell_names <- rownames(metadata)
last_digits <- as.numeric(substr(cell_names, nchar(cell_names), nchar(cell_names)))

# 动态调整 cell_type 的值 0-8 9-17 18-26 27-35 36-44
# 动态调整 cell_type 的值 0-9 10-19 20-29 30-39 40-49
metadata$adjusted_celltype <- metadata$celltype + 10 * last_digits

# 将修改后的元数据更新到 Seurat 对象
pbmc@meta.data <- metadata


# 计算不同聚类的平均表达量
# group.by = "celltype" 按照列 celltype 分组
# 通过AggregateExpression生成伪样本计数矩阵
pseudobulk <- Seurat::AggregateExpression(
  pbmc,
  group.by = "adjusted_celltype",  # 分组字段
  slot = "counts",                     # 使用原始counts
  fun = "sum"                          # 求和聚合
)
count_matrix <- pseudobulk$RNA          # 提取计数矩阵

# 从伪样本列名解析元数据（假设列名为"sample_celltype"）
sample_info <- data.frame(
  celltype = colnames(count_matrix),  # 直接使用列名
  row.names = colnames(count_matrix)
)

# 从GTF文件加载基因长度（需替换为实际数据）
library(GenomicFeatures)
txdb <- txdbmaker::makeTxDbFromGFF("{PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf")
gene_lengths <- transcriptLengths(txdb, with.cds_len=TRUE) %>%
  group_by(gene_id) %>%
  summarise(length = max(cds_len))  # 取最长转录本长度

# 确保基因名称与计数矩阵一致
rownames(gene_lengths) <- gene_lengths$gene_id
gene_lengths <- gene_lengths[rownames(count_matrix), ]

library(DESeq2)

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ celltype  # 根据分析目标设定design
)

# 添加基因长度信息
mcols(dds)$basepairs <- gene_lengths$length

# 方法1：直接计算（未标准化）
fpkm_matrix <- fpkm(dds, robust = TRUE)  # robust=TRUE对长度做中位数校正

# 方法2：先运行标准化（推荐）
dds <- estimateSizeFactors(dds)
fpkm_matrix <- fpkm(dds, robust = TRUE)

# 保存为CSV
write.csv(fpkm_matrix, "./02_processData/00_02_pseudobulk_fpkm.csv")


