#!/usr/bin/env python3

import anndata
import pandas as pd
import matplotlib.pyplot as plt

# ===============================
# 1. 读取 h5ad
# ===============================
adata = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
adata.obs_names_make_unique()

# ===============================
# 2. 注释 celltype
# ===============================
cluster2annotation = {
    "0": "5",
    "1": "5",
    "2": "3",
    "3": "4",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "2",
    "8": "7",
    "9": "1"
}

adata.obs["celltype"] = (
    adata.obs["leiden_0.4"]
    .astype(str)
    .map(cluster2annotation)
    .astype("category")
)

# ===============================
# 3. 读取 barcode
# ===============================
cell02dpi1 = pd.read_csv(
    "02dpi_1_barcodes.tsv",
    header=None,
    sep="\t"
)[0].astype(str).tolist()

cell02dpi2 = pd.read_csv(
    "02dpi_2_barcodes.tsv",
    header=None,
    sep="\t"
)[0].astype(str).tolist()

# ===============================
# 4. 处理 obs_names
# adata.obs_names: AAAC...-1-0
# barcode 文件: AAAC...-1
# 去掉最后的 -0 / -1 / -2 ...
# ===============================
obs = adata.obs.copy()
obs["cell_id"] = obs.index
obs["barcode"] = obs["cell_id"].str.replace(r"-[0-9]+$", "", regex=True)

obs["group"] = "other"
obs.loc[obs["barcode"].isin(cell02dpi1), "group"] = "cell02dpi1"
obs.loc[obs["barcode"].isin(cell02dpi2), "group"] = "cell02dpi2"

# ===============================
# 5. 加入 UMAP 坐标
# ===============================
obs["UMAP_1"] = adata.obsm["X_umap"][:, 0]
obs["UMAP_2"] = adata.obsm["X_umap"][:, 1]

# ===============================
# 6. 统计 celltype == 7 数量
# ===============================
n_02dpi1_total = (obs["group"] == "cell02dpi1").sum()
n_02dpi2_total = (obs["group"] == "cell02dpi2").sum()

n_02dpi1_ct7 = ((obs["group"] == "cell02dpi1") & (obs["celltype"].astype(str) == "7")).sum()
n_02dpi2_ct7 = ((obs["group"] == "cell02dpi2") & (obs["celltype"].astype(str) == "7")).sum()

print(f"cell02dpi1 total cells: {n_02dpi1_total}")
print(f"cell02dpi1 celltype==7 cells: {n_02dpi1_ct7}")

print(f"cell02dpi2 total cells: {n_02dpi2_total}")
print(f"cell02dpi2 celltype==7 cells: {n_02dpi2_ct7}")

# ===============================
# 7. 绘图函数
# ===============================
def plot_highlight_umap(
    obs,
    group_name,
    outfile_prefix
):
    target = (obs["group"] == group_name) & (obs["celltype"].astype(str) == "7")
    other = ~target

    plt.figure(figsize=(5, 4.5))

    plt.scatter(
        obs.loc[other, "UMAP_1"],
        obs.loc[other, "UMAP_2"],
        s=2,
        c="lightgrey",
        alpha=0.5,
        linewidths=0
    )

    plt.scatter(
        obs.loc[target, "UMAP_1"],
        obs.loc[target, "UMAP_2"],
        s=6,
        c="#D7B76A",
        alpha=0.9,
        linewidths=0
    )

    plt.title(f"{group_name}: celltype 7", fontsize=14, fontweight="bold")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()

    plt.savefig(f"{outfile_prefix}.pdf")
    plt.savefig(f"{outfile_prefix}.png", dpi=300)
    plt.close()

# ===============================
# 8. 分别绘制两个 UMAP
# ===============================
plot_highlight_umap(
    obs,
    group_name="cell02dpi1",
    outfile_prefix="cell02dpi1_celltype7_umap"
)

plot_highlight_umap(
    obs,
    group_name="cell02dpi2",
    outfile_prefix="cell02dpi2_celltype7_umap"
)

print("Done.")