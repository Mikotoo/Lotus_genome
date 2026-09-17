import anndata
import os
import matplotlib.pyplot as plt
import scanpy as sc
import scipy.io as sio

data_cluster = anndata.read_h5ad("./02_processData/data_clusterd_all_gene_4000.h5ad")
data_cluster.obs_names_make_unique()
selected_data = data_cluster

#############cell type
#Xylem:1   Phloem:2   Pericycle:3   Vascular:4   Cortex:5   Meristem:6    Early symbiotic signaling cells:7 
cluster2annotation = {
    '0': '5',
    '1': '5',
    '2': '3',
    '3': '4',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '2',
    '8': '7',
    '9': '1'
}

selected_data.obs['celltype'] = selected_data.obs['leiden_0.4'].map(cluster2annotation).astype('category')

#############################################

cellinfo=selected_data.obs
geneinfo=selected_data.var
mtx2=selected_data.layers['counts'].T
mtx1=selected_data.X.T

cellinfo.to_csv("./02_processData/zz_01-cellinfo.csv")
geneinfo.to_csv("./02_processData/zz_02-geneinfo.csv")
sio.mmwrite("./02_processData/zz_03-sparse_matrix_counts.mtx",mtx2)
sio.mmwrite("./02_processData/zz_04-sparse_matrix.mtx",mtx1)

