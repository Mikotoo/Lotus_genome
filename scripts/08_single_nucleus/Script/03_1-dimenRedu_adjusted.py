import sys
sys.path.append('{CONDA_SINGLECELL}/lib/python3.9/site-packages/Scripts')
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
import Scripts
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
clustree = importr('clustree')
ggplot2 = importr('ggplot2')

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    print(f"==== Running n_top_genes = {n} ====")

    in_file = f"02_processData/data_integrated_{n}.h5ad"

    data_integrated = anndata.read_h5ad(in_file)
    data = data_integrated

    palette = ["#1F77B4","#FF7F0E"]


    categories = data.obs['Sample'].cat.categories
    color_map = {category: color for category, color in zip(categories, palette)}

    ###降维
    ## UMAP
    # 通过降低n_neighbors和min_dist让各类之间的边界更清晰
    sc.pp.neighbors(data, n_neighbors=10, use_rep='X_scVI')
    sc.tl.umap(data, min_dist=0.3)
    data.obsm["scvi_umap"]=data.obsm["X_umap"]

    fig,ax = plt.subplots(figsize=(10,6),dpi=300)
    sc.pl.embedding(data,basis='scvi_umap',color="Sample",palette=color_map,show=False, ax=ax)
    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
    out_file = f"03_figures/dimenRedu_{n}_adjusted.png"
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()

    ##聚类
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

    for res in resolutions:
        sc.tl.leiden(data, resolution=res, key_added=f'leiden_{res}',random_state=0)


    fig,axes=plt.subplots(3,2,figsize=(30,15))
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.9,right=0.8,hspace=0.2,wspace=0.4)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0,0],color="leiden_0.2")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0,1],color="leiden_0.4")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1,0],color="leiden_0.6")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1,1],color="leiden_0.8")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[2,0],color="leiden_1.0")
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[2,1],color="leiden_1.2")
    out_file = f"03_figures/cluster_{n}_adjusted.png"
    plt.savefig(out_file)
    plt.close()

    cluster_df = pd.DataFrame({f'leiden_{res}': data.obs[f'leiden_{res}'] for res in resolutions})
    r_cluster_df = pandas2ri.py2rpy(cluster_df)
    ro.globalenv['cluster_df'] = r_cluster_df
    ro.globalenv['n'] = n
    r_code = """
    library(clustree)
    library(ggplot2)
    p <- clustree(cluster_df, prefix = "leiden_")
    out_file = paste0("03_figures/clustree_plot_", n, ".pdf")
    ggsave(out_file, plot = p, device = "pdf", width=16, height=8)
    """
    ro.r(r_code)

    # final
    fig,axes=plt.subplots(2,1,figsize=(8,10),dpi=300)
    plt.subplots_adjust(left=0.1,bottom=0.1,top=0.9,right=0.8,hspace=0.2,wspace=0.25)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[0],color="Sample",palette=color_map,show=False)
    sc.pl.embedding(data,basis='scvi_umap',ax=axes[1],color="leiden_0.4",legend_loc="on data",show=False)
    axes[0].set_title("Samples")
    axes[1].set_title("Clusters")
    out_file = f"03_figures/final_cluster_{n}_adjusted.png"
    plt.savefig(out_file)
    plt.close()

    data.write_h5ad(f"02_processData/data_clusterd_{n}_adjusted.h5ad")

    # cluster_statitcs
    clusters = data.obs["leiden_0.4"]
    samples = data.obs['Sample']
    cluster_counts = clusters.value_counts().sort_index()
    cluster_sample_counts=pd.crosstab(clusters,samples)
    cluster_sample_props = cluster_sample_counts.div(cluster_sample_counts.sum(axis=1), axis=0)


    fig, axs = plt.subplots(nrows=2, figsize=(12, 16), dpi=300)
    cluster_counts.plot(kind='bar', ax=axs[0], color='grey')
    axs[0].set_title('Total Number of Cells in Each Cluster')
    axs[0].set_xlabel('Cluster')
    axs[0].set_ylabel('Number of Cells')

    cluster_sample_props.plot(kind='bar', stacked=True, color=[color_map[color] for color in categories], ax=axs[1])
    axs[1].set_title('Proportion of Samples in Each Cluster')
    axs[1].set_xlabel('Cluster')
    axs[1].set_ylabel('Proportion')

    handles, labels = axs[1].get_legend_handles_labels()
    axs[1].legend(handles, labels, title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')

    out_file = f"03_figures/cluster_statitcs_{n}_adjusted.png"
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()