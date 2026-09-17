import anndata
import os
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.io as sio

os.makedirs("04_Anno", exist_ok=True)

data_cluster = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
selected_data = data_cluster
selected_data.obs_names_make_unique()
#############cell type
cluster2annotation = {
    '0': '0',
    '1': '1',
    '2': '2',
    '3': '3',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '7',
    '8': '8',
    '9': '9'
}

selected_data.obs['celltype'] = selected_data.obs['leiden_0.4'].map(cluster2annotation).astype('category')
#############################################

cellinfo=selected_data.obs
geneinfo=selected_data.var
mtx2=selected_data.layers['counts'].T
mtx1=selected_data.X.T

cellinfo.to_csv("04_Anno/cluster_01-cellinfo.csv")
geneinfo.to_csv("04_Anno/cluster_02-geneinfo.csv")
sio.mmwrite("04_Anno/cluster_03-sparse_matrix_counts.mtx",mtx2)
sio.mmwrite("04_Anno/cluster_04-sparse_matrix.mtx",mtx1)

