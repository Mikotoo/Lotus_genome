import sys
import anndata
import scanpy as sc
import scvi

# 读取数据
data_scaled = anndata.read_h5ad("./02_processData/data_all.h5ad")

# 循环运行不同的 n_top_genes
for n in [1000, 2000, 3000, 4000, 5000]:

    print(f"==== Running n_top_genes = {n} ====")

    # 每次都复制数据，避免被上一次处理修改
    data = data_scaled.copy()

    # 选择高变基因
    sc.experimental.pp.highly_variable_genes(
        data, 
        flavor="pearson_residuals",
        layer="counts",
        batch_key="Sample",
        n_top_genes=n,
        subset=True,
        inplace=True,
    )

    # 设置 scVI anndata
    scvi.model.SCVI.setup_anndata(
        data,
        layer="counts",
        batch_key="Sample",
        categorical_covariate_keys=["Sample"],
        continuous_covariate_keys=["total_counts"]
    )

    # 固定随机种子
    scvi.settings.seed = 7

    # 初始化模型
    model = scvi.model.SCVI(data)

    # 训练模型
    model.train()

    # 保存 latent space
    data.obsm['X_scVI'] = model.get_latent_representation()

    # 保存 normalized expression
    data.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)

    # 输出文件名
    out_file = f"02_processData/data_integrated_{n}.h5ad"
    data.write_h5ad(out_file)

    print(f"==== Saved: {out_file} ====\n")
