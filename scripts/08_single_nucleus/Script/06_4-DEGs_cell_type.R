library(dplyr)
library(Seurat)
library(patchwork)
library(Biobase)
library(ggplot2)

#读取注释好的 h5ad 文件转换后的数据
pbmc <- readRDS("./02_processData/zz_00_pea.rds")
### 计算cluster的marker gene
#Idents(pbmc) <- pbmc@meta.data$adjusted_celltype
##所有cluster
#all_markers <- FindAllMarkers(
#  object = pbmc
#)
#写出marker gene
#write.csv(all_markers, file = "zz_07_celltype_difftime_markers_gene.csv", row.names = FALSE)

### 计算cluster的marker gene
# 定义 cluster 组别
cluster_groups <- list(
  Xylem = c(1, 11),
  Phloem = c(2, 12),
  Pericycle = c(3, 13),
  Vascular = c(4, 14),
  Cortex = c(5, 15),
  Meristem = c(6, 16),
  Early = c(7, 17)
)

# 设置细胞标识
Idents(pbmc) <- pbmc@meta.data$adjusted_celltype

# 初始化结果列表
marker_results <- list()

# 计算组内差异基因，以dpi01为背景
group_names <- names(cluster_groups)
for (i in 1:(length(group_names))) {
    group_name <- names(cluster_groups)[i]
    clusters <- cluster_groups[[group_name]]
    markers <- FindMarkers(
      pbmc,
      ident.1 = clusters[1],
      ident.2 = clusters[2]
      )
    markers <- markers %>%
      mutate(Significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ))

      # 保存结果
    comparison_name <- paste(group_name, "02dpivs00dpi", sep = "_")
    marker_results[[comparison_name]] <- markers

    ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-Value"
      )

    ggsave(
      filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
      width = 6, height = 5)

    # 输出结果到文件
    write.csv(markers, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

    # 打印提示信息
    cat("已完成比较：", comparison_name, "\n")
}

# 初始化结果列表
marker_results <- list()

# 遍历所有组合
group_names <- names(cluster_groups)
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    # 获取当前组合的组名和 cluster 编号
    group1 <- group_names[i]
    group2 <- group_names[j]
    clusters1 <- cluster_groups[[group1]]
    clusters2 <- cluster_groups[[group2]]

    # 计算差异基因
    markers <- FindMarkers(
      pbmc,
      ident.1 = clusters1,
      ident.2 = clusters2 
    )
    markers <- markers %>%
      mutate(Significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ))

    # 保存结果
    comparison_name <- paste(group1, "vs", group2, sep = "_")
    marker_results[[comparison_name]] <- markers

    # 绘制火山图
    ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-Value"
      )

    ggsave(
    filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
    width = 6, height = 5)

    # 输出结果到文件
    write.csv(markers, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

    # 打印提示信息
    cat("已完成比较：", comparison_name, "\n")
  }
}

###############################################
### 额外分析：Early symbiotic signaling cells vs 非 Early symbiotic signaling cells
###############################################

cat("开始 Early symbiotic signaling cells vs Others 分析...\n")

# 定义 Early symbiotic signaling cells 的 clusters
Early_clusters <- cluster_groups$Early  # 即 c(7, 17)

# 所有 clusters
all_clusters <- unique(unlist(cluster_groups))

# 除 Early symbiotic signaling cells 外的其他 clusters
other_clusters <- setdiff(all_clusters, Early_clusters)

# 执行差异分析（Early symbiotic signaling cells 为 ident.1，上调表示 Early symbiotic signaling cells 相对 others 上调）
markers_inf <- FindMarkers(
  pbmc,
  ident.1 = Early_clusters,
  ident.2 = other_clusters
)

markers_inf <- markers_inf %>%
  mutate(Significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1  ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# 保存结果
comparison_name <- "Early_symbiotic_signaling_cells_vs_Others"
marker_results[[comparison_name]] <- markers_inf

# 绘制火山图
p_inf <- ggplot(markers_inf, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Early symbiotic signaling cells vs Others",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  )

ggsave(
  filename = paste0("./03_figures/", comparison_name, "_volcano_plots.pdf"),
  width = 6, height = 5
)

# 输出到 CSV
write.csv(markers_inf, file = paste0("./02_processData/", comparison_name, "_markers.csv"))

cat("已完成比较：", comparison_name, "\n")
