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

# 数据归一化
# NormalizeData 对基因表达数据进行归一化处理
pbmc <- NormalizeData(pbmc)

# 筛选高变基因
# FindVariableFeatures 用于识别在不同细胞间具有高度可变性的基因
pbmc <- FindVariableFeatures(pbmc)

# 计算不同聚类的平均表达量
# group.by = "celltype" 按照列 celltype 分组
dat <- AverageExpression(pbmc, group.by = "celltype")$RNA

# 将计算的平均表达量保存为 CSV 文件
write.csv(dat, file = "./02_processData/zz_05-AverageExpression_celltype.csv", row.names = TRUE)

### 计算两个cluster的差异表达情况
#Idents(pbmc) <- pbmc@meta.data$celltype
#cluster5.markers <- FindMarkers(pbmc, ident.1 = c(5,6), ident.2 = c(0,3))

### 将对象拆分为不同时期处理，比较不同时期的表达差异
# 假设 pbmc 是一个 Seurat 对象，提取元数据
metadata <- pbmc@meta.data

# 提取细胞名称，并获取最后一位数字
cell_names <- rownames(metadata)
last_digits <- as.numeric(substr(cell_names, nchar(cell_names), nchar(cell_names)))

# 动态调整 cell_type 的值 0-9 10-19 20-29 30-39 40-49
metadata$adjusted_celltype <- metadata$celltype + 10 * last_digits

# 将修改后的元数据更新到 Seurat 对象
pbmc@meta.data <- metadata

# 查看调整后的结果
#head(metadata[, c("leiden_0.6", "adjusted_leiden_0.6")])

# 计算基因表达
dat <- AverageExpression(pbmc, group.by = "adjusted_celltype")$RNA
write.csv(dat, file = "./02_processData/zz_06-AverageExpression_celltype_difftime.csv", row.names = TRUE)

saveRDS(pbmc, file = "./02_processData/zz_00_pea.rds")


################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$adjusted_celltype
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_07_celltype_difftime_markers_gene.csv", row.names = FALSE)

### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$celltype
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_07_celltype_markers_gene.csv", row.names = FALSE)

################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
Idents(pbmc) <- pbmc@meta.data$leiden_0.4
##所有cluster
all_markers <- FindAllMarkers(
  object = pbmc, 
  min.pct = 0.25, 
  logfc.threshold = 0.25, 
  only.pos = TRUE
)
#写出marker gene
write.csv(all_markers, file = "./02_processData/zz_08_cluster_markers_gene.csv", row.names = FALSE)