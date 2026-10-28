import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
import scripts
import matplotlib.pyplot as plt

data_clustered = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data=data_clustered

marker_genes = {
    "Cortex": [
        "Lotus_GifuT2T_01G0461100",
        "Lotus_GifuT2T_02G1487100",
    ],
    "Epidermis": [
        "Lotus_GifuT2T_05G2864900"
    ],
    "Meristem": [
        "Lotus_GifuT2T_03G1828600",
        "Lotus_GifuT2T_04G2253600",
        "Lotus_GifuT2T_04G2253800"
    ],
    "Pericycle": [
        "Lotus_GifuT2T_02G1146000",
        "Lotus_GifuT2T_04G2574000"
    ],
    "Phloem": [
        "Lotus_GifuT2T_03G2130600",
        "Lotus_GifuT2T_05G2958500",
    ],
    "Vascular": [
        "Lotus_GifuT2T_01G0541600", #
        "Lotus_GifuT2T_01G0093700",
        "Lotus_GifuT2T_01G0703300",
        "Lotus_GifuT2T_01G0626500"
    ],
    "Xylem": [
        "Lotus_GifuT2T_01G0273300",
        "Lotus_GifuT2T_04G2644400",
    ],
    "Infected cell": [
        "Lotus_GifuT2T_04G2692200",
        "Lotus_GifuT2T_01G0784100",
        "Lotus_GifuT2T_03G1515300",
        "Lotus_GifuT2T_05G3052700",
    ],
    "Nodule primordium": [
        "Lotus_GifuT2T_01G0025100",
        "Lotus_GifuT2T_04G2544700"
    ],
    "LjLb": [
        "Lotus_GifuT2T_05G2767100",
        "Lotus_GifuT2T_05G2767300",
        "Lotus_GifuT2T_05G2792500"
    ]
}

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
# cluster_order = ['3', '2', '6','12','10','7','8','0','1','4','5','9','11','13','14']  # 替换为你的实际顺序
# cluster_order = ['8', '6', '1', '0', '5', '4', '3', '2', '7', '9']
cluster_order = ['0', '1', '5', '8', '6', '2', '7', '3', '4', '9']

# 将 leiden_0.4 列转换为有序分类
data.obs['leiden_0.4'] = pd.Categorical(
    data.obs['leiden_0.4'],
    categories=cluster_order,  # 强制指定顺序
    ordered=True
)

# 创建画布
fig, axes = plt.subplots(2, 1, figsize=(16, 12), dpi=300)

# 绘制点图（无需 groups_order 参数）
sc.pl.dotplot(
    data,
    ax=axes[0],
    groupby="leiden_0.4",
    var_names=marker_genes_in_data,
    dendrogram=False,
    standard_scale="var"
)

axes[1].axis("off")

plt.savefig("03_figures/markerAnno_filter_dot_v5_own_0.4.pdf")
