import anndata

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    data_concatenated = anndata.read_h5ad("./02_processData/data_all.h5ad")

    print(f"==== Running n_top_genes = {n} ====")

    in_file = f"02_processData/data_clusterd_{n}.h5ad"
    data_cluster = anndata.read_h5ad(in_file)

    # Copy 'obs' and 'var' if necessary
    data_concatenated.obs = data_cluster.obs
    data_concatenated.uns = data_cluster.uns
    data_concatenated.obsm = data_cluster.obsm
    #data_concatenated.layers = data_cluster.layers
    data_concatenated.obsp = data_cluster.obsp

    data_concatenated.write_h5ad(f"02_processData/data_clusterd_all_gene_{n}.h5ad")