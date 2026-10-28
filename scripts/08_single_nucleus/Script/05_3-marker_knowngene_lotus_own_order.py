import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scripts
import matplotlib.pyplot as plt

data_clustered = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data=data_clustered


### marker基因注释
marker_genes = {"":[
"Lotus_GifuT2T_06G3179800",
"Lotus_GifuT2T_01G0061000",
"Lotus_GifuT2T_01G0087900",
"Lotus_GifuT2T_01G0128800"]}

# check if the markers are in the data
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in data.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found


sc.tl.dendrogram(data,'leiden_0.4',use_rep="X_scVI")

print(marker_genes_in_data)

# 指定聚类顺序，按实际需求调整列表顺序
cluster_order = ['0', '1', '5', '8', '6', '2', '7', '3', '4', '9'] # 替换为你的实际顺序

# 将 leiden_0.4 列转换为有序分类
data.obs['leiden_0.4'] = pd.Categorical(
    data.obs['leiden_0.4'],
    categories=cluster_order,  # 强制指定顺序
    ordered=True
)

categories = data.obs['Sample'].cat.categories
for sample in categories:
    selected_cells = data.obs['Sample'] == sample
    sample_data = data[selected_cells].copy()

    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in sample_data.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found

    sc.tl.dendrogram(sample_data,'leiden_0.4',use_rep="X_scVI")

    fig, axes=plt.subplots(2,1,figsize=(16,12),dpi=300)
    sc.pl.dotplot(
        sample_data,
        ax=axes[0],
        groupby="leiden_0.4",
        var_names=marker_genes_in_data,
        dendrogram=False,
        standard_scale="var",
    )
    axes[1].axis("off")

    plt.savefig(f"03_figures/knowngene_dot_{sample}_1_own_0.4.pdf")

# ================================
#  绘制 marker 基因的 UMAP 分布图
# ================================
# 四个 marker 基因列表
genes_to_plot = [
    "Lotus_GifuT2T_06G3179800",
    "Lotus_GifuT2T_01G0061000",
    "Lotus_GifuT2T_01G0087900",
    "Lotus_GifuT2T_01G0128800"
]

# 逐个 sample 绘制 UMAP
for sample in categories:
    selected_cells = data.obs['Sample'] == sample
    sample_data = data[selected_cells].copy()

    # 筛选当前 sample 中确实存在的基因
    valid_genes = [g for g in genes_to_plot if g in sample_data.var.index]

    if len(valid_genes) == 0:
        print(f"⚠️ Sample {sample} 中没有任何 marker 基因")
        continue

    # ---- 创建一个 2×2 的 panel 图 ----
    n = len(valid_genes)
    ncols = 2
    nrows = int(np.ceil(n / ncols))

    fig, axs = plt.subplots(nrows, ncols, figsize=(10, 10), dpi=300)
    axs = axs.flatten()

    for i, gene in enumerate(valid_genes):
        sc.pl.umap(
            sample_data,
            color=gene,
            ax=axs[i],
            show=False,
            color_map="viridis",  # 可换 plasma, magma
            s=10
        )
        axs[i].set_title(gene, fontsize=10)

    # 隐藏多余子图
    for j in range(i+1, len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    plt.savefig(f"03_figures/knowngene_umap_{sample}_4genes.pdf")
    plt.close()

print("UMAP 分布图已全部绘制完成！")
