#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# =============================
# 1. 读取命令行参数
# =============================
if len(sys.argv) < 2:
    print("请提供 n 值，例如：python 03_1-cluster_statitcs.py 4000")
    sys.exit(1)

n = sys.argv[1]

# =============================
# 2. 读取 .h5ad 数据
# =============================
h5_file = f"02_processData/data_clusterd_{n}.h5ad"
adata = anndata.read_h5ad(h5_file)

df = adata.obs.copy()

# =============================
# 3. 提取 cluster & Sample
# =============================
df["cluster"] = df["leiden_0.4"].astype(str)
df["sample"] = df["Sample"].astype(str)

# =============================
# 4. cluster 指定顺序
# =============================
cluster_order = ['9', '3', '4', '7', '2', '6', '8', '0', '1', '5']
df["cluster"] = pd.Categorical(df["cluster"], categories=cluster_order, ordered=True)

# =============================
# 5. 每个 cluster 的总细胞数
# =============================
cluster_counts = df.groupby("cluster").size().reset_index(name="total_cells")

# =============================
# 6. 各 Sample 在每 cluster 的真实细胞数（raw counts）
# =============================
cluster_sample_count = df.groupby(["cluster", "sample"]).size().reset_index(name="count")

cluster_sample_count = cluster_sample_count.merge(cluster_counts, on="cluster", how="left")

# =============================
# 7. pivot 创建 raw count 堆叠矩阵
# =============================
plot_df = cluster_sample_count.pivot_table(
    index="cluster", columns="sample", values="count", fill_value=0
)

plot_df = plot_df.reindex(cluster_order)

# =============================
# 8. 绘图（横向 raw count stacked bar）
# =============================
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

# 自定义颜色
palette = ["#1F77B4", "#FF7F0E"]
samples = plot_df.columns.tolist()
colors = (palette * (len(samples) // len(palette) + 1))[:len(samples)]

# 横向堆叠条形图（使用真实数量）
plot_df.plot(
    kind="barh",
    stacked=True,
    ax=ax,
    color=colors,
    width=0.7
)

# =============================
# 9. 添加 cluster 总细胞数注释（右侧）
# =============================
for idx, row in cluster_counts.iterrows():
    cluster = row["cluster"]
    total_cells = row["total_cells"]
    y_pos = cluster_order.index(cluster)

    ax.text(
        plot_df.sum(axis=1).max() * 1.02,   # 在右侧一点点
        y_pos,
        str(total_cells),
        va="center",
        ha="left",
        fontsize=9,
        fontweight="bold"
    )

# =============================
# 10. 删除外边框，仅保留坐标轴
# =============================
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.tick_params(axis='both', which='both', length=4)

# =============================
# 图形标签
# =============================
ax.set_title(f"Cluster Composition (raw cell counts, n_top_genes = {n})")
ax.set_xlabel("Number of Cells")
ax.set_ylabel("Cluster")
ax.legend(samples, title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()

# =============================
# 11. 保存 PDF
# =============================
os.makedirs("03_figures", exist_ok=True)
out_file = f"03_figures/cluster_stacked_raw_{n}.pdf"
plt.savefig(out_file, bbox_inches="tight")

print("完成绘图:", out_file)
