#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os

# --------- 目录路径（你的 Arial 字体存放目录） ------------
font_dir = "{HOME_HX}/.fonts/Arial"
for f in os.listdir(font_dir):
    if f.endswith(".ttf") or f.endswith(".TTF"):
        font_manager.fontManager.addfont(os.path.join(font_dir, f))

# --------- 强制 matplotlib 使用 Arial 字体 ------------
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42   # 可编辑文本
matplotlib.rcParams['ps.fonttype'] = 42

# ------------ 保证AI中可编辑的Arial字体 ------------
matplotlib.rcParams['pdf.fonttype'] = 42   # 不要转轮廓
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
# ------------------------------------------------------

import os
import pandas as pd
import anndata
import scanpy as sc
import numpy as np

# ==========================
# 0. Read h5ad
# ==========================
data = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")

# ==========================
# 1. Cluster order
# ==========================
# cluster_order = ['9', '3', '4', '7', '2', '6', '8', '0', '1', '5']
cluster_order = ['5', '1', '0', '8', '6', '2', '7', '4', '3', '9']

data.obs["leiden_0.4"] = pd.Categorical(
    data.obs["leiden_0.4"],
    categories=cluster_order,
    ordered=True
)

# ==========================
# 2. Read marker gene CSV
# ==========================
csv_path = "./04_Anno/cluster_06_cluster_markers_gene.csv"
df = pd.read_csv(csv_path)

df["cluster"] = df["cluster"].astype(str)

# ---- 统一 gene 名格式：替换 - 为 _ ----
df["gene"] = df["gene"].str.replace("-", "_")

# ==========================
# 3. special cluster 8 with corrected gene names
# ==========================
special_cluster8_genes_raw = [
    "Lotus_GifuT2T_04G2692200",
    "Lotus-GifuT2T-01G0002200",
    "Lotus-GifuT2T-02G1393700"
]

# 替换为与矩阵一致的格式
special_cluster8_genes = [g.replace("-", "_") for g in special_cluster8_genes_raw]

cluster_top_genes = {}

# ==========================
# 4. Process each cluster
# ==========================
for cl in cluster_order:
    sub = df[df["cluster"] == cl]

    # Cluster8 special handling
    if cl == "8":
        cluster_top_genes[cl] = special_cluster8_genes
        continue

    # Apply filtering rules
    sub_filtered = sub[
        (sub["avg_log2FC"] > 0.25) &
        (sub["pct.1"] > 0.25) &
        (sub["p_val_adj"] < 0.05)
    ]

    if sub_filtered.shape[0] == 0:
        print(f"⚠️ Cluster {cl} has no genes passing filter")
        continue

    # Sort by p_val_adj → avg_log2FC
    sub_filtered = sub_filtered.sort_values(
        by=["p_val_adj", "avg_log2FC"],
        ascending=[True, False]
    )

    # Select top3
    cluster_top_genes[cl] = list(sub_filtered["gene"].head(3))

# print final used genes
print("\n最终用于绘图的基因：")
for cl in cluster_order:
    if cl in cluster_top_genes:
        print(f"Cluster {cl}: {cluster_top_genes[cl]}")
    else:
        print(f"Cluster {cl}: 无基因符合条件")

# ==============================================================
# 5. DOTPLOT ALL CLUSTERS
# ==============================================================

all_genes = []
for cl in cluster_order:
    if cl in cluster_top_genes:
        all_genes.extend(cluster_top_genes[cl])

# remove genes not in var
all_genes = [g for g in all_genes if g in data.var.index]

# 如果仍然没有基因，直接退出防止报错
if len(all_genes) == 0:
    raise ValueError("❌ dotplot 无基因可绘制，检查基因名格式是否都正确匹配 data.var.index")

# 自动根据基因数调整画布宽度
n_genes = len(all_genes)
fig_width = max(12, n_genes * 0.8)

print(f"绘制 dotplot，基因数 {n_genes}，figsize=({fig_width}, 10)")

plt.figure(figsize=(fig_width, 10), dpi=300)

sc.pl.dotplot(
    data,
    var_names=all_genes,
    groupby="leiden_0.4",
    standard_scale="var",
    dendrogram=False,
    show=False
)

plt.savefig("./04_Anno/dotplot_all_clusters_top3.pdf", dpi=300)
plt.close()

print("\nDOTPLOT Done!")

# ==============================================================
# 6. UMAP for each cluster
# ==============================================================

out_root = "./04_Anno/UMAP_top/"
os.makedirs(out_root, exist_ok=True)

red_gray_cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
    "redgray",
    ["lightgray", "red"]
)

for cl in cluster_order:
    if cl not in cluster_top_genes:
        continue

    cl_dir = os.path.join(out_root, f"cluster_{cl}")
    os.makedirs(cl_dir, exist_ok=True)

    for gene in cluster_top_genes[cl]:
        if gene not in data.var.index:
            print(f"⚠️ {gene} not found in data")
            continue

        sc.pl.umap(
            data,
            color=gene,
            show=False,
            color_map=red_gray_cmap,
            s=10
        )

        plt.title(f"Cluster {cl} - {gene}")
        plt.savefig(os.path.join(cl_dir, f"{gene}_UMAP.pdf"), dpi=300)
        plt.close()

print("\nUMAP All Done!")
