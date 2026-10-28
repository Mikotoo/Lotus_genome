# 加载必要的包
# dplyr: 数据操作
# Seurat: 单细胞数据分析
# patchwork: 用于将多个图形组合成一个布局
# Biobase: 生物信息学分析工具
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)

# 读取注释好的 h5ad 文件转换后的数据
# 读取细胞信息（cellinfo）、基因信息（geneinfo）和表达矩阵（counts）
cellinfo <- read.csv("04_Anno/cluster_01-cellinfo.csv", row.names = 1)  # 细胞元数据
geneinfo <- read.csv("04_Anno/cluster_02-geneinfo.csv", row.names = 1)  # 基因元数据
counts <- Matrix::readMM(file = "04_Anno/cluster_03-sparse_matrix_counts.mtx")  # 稀疏表达矩阵

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
write.csv(dat, file = "04_Anno/cluster_05-AverageExpression_celltype.csv", row.names = TRUE)

### 计算两个cluster的差异表达情况
#Idents(pbmc) <- pbmc@meta.data$celltype
#cluster5.markers <- FindMarkers(pbmc, ident.1 = c(5,6), ident.2 = c(0,3))

saveRDS(pbmc, file = "04_Anno/cluster_00_pea.rds")

################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
#读取注释好的 h5ad 文件转换后的数据
seurat_obj <- readRDS("04_Anno/cluster_00_pea.rds")
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
write.csv(all_markers, file = "04_Anno/cluster_06_cluster_markers_gene.csv", row.names = FALSE)

all_markers <- FindAllMarkers(
  object = pbmc
)
all_markers <- all_markers[all_markers$p_val_adj<0.05&abs(all_markers$avg_log2FC)>1,]
#写出marker gene
write.csv(all_markers, file = "04_Anno/cluster_06_cluster_markers_gene1.csv", row.names = FALSE)