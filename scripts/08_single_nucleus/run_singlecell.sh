#!/bin/bash
#CSUB -J singlecell
#CSUB -q c01
#CSUB -o singlecell.out
#CSUB -e singlecell.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# Steps 01-03 and 06_2-06_3 of the single-nucleus pipeline, which were
# previously one six-line wrapper each. The longer steps keep their own
# scripts: 00_cellranger.sh, 00_mkref.sh, 05_cluster_DEGs_finder.sh,
# 06_1.marker_gene_finder.sh, 06_1.Individual_clusters_marker_gene_finder.sh
# and 07_batch_multiple_scTenifoldKnk.sh.
#
# Run from 08_single_nucleus/:
#   bash run_singlecell.sh

# 加载 conda
source ~/anaconda3/etc/profile.d/conda.sh

# 激活环境
conda activate singlecell

# ── 01 QC ────────────────────────────────────────────────────────────────
python ./Script/01_annaDataQC.py

# ── 02 integration ───────────────────────────────────────────────────────
python ./Script/02_integrate.py

# ── 03 dimensionality reduction and clustering ───────────────────────────
python ./Script/03_1-dimenRedu_adjusted.py

# ── 06_2 marker annotation ───────────────────────────────────────────────
# 已有的文献的marker gene效果不好，故改用先发掘各簇DEG，在比对到大豆marker gene的思路。
# python ./Script/05_3-marker_knowngene_lotus_own_order.py
python ./Script/05_3-markerAnno_Ps_v5_own.py

# ── 06_3 cell-type expression and DEGs ───────────────────────────────────
# The three calls below feed this one and are commented out in the original run.
# python ./Script/06_1-Seurat_Exp_h5ad_to_rds.py
# Rscript ./Script/06_2-Exp_h5ad_to_rds_cell_type.R
# Rscript ./Script/06_3-agg_to_fpkm.R
Rscript ./Script/06_4-DEGs_cell_type.R
