#!/bin/bash
#CSUB -J clusterDEG
#CSUB -q c01
#CSUB -o clusterDEG.out
#CSUB -e clusterDEG.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 加载 conda
source ~/anaconda3/etc/profile.d/conda.sh

# 激活环境
conda activate singlecell

# 运行脚本
python ./Script/05_1-Seurat_Exp_h5ad_to_rds.py
Rscript ./Script/05_2-Exp_h5ad_to_rds_cluster.R

# 对每簇的marker genes的功能进行富集
#========================
# 0. 基本路径设置
#========================
output=./04_Anno/cluster_makers_genes_enrichment
mkdir -p "$output"

script={SOFTWARE_HX}/Script/RNA-seq
anno={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_Gifu.emapper.annotations
file=./04_Anno/cluster_06_cluster_markers_gene1.csv


#========================
# 1. 获取所有 cluster 编号
#========================
clusters=$(tail -n +2 "$file" \
  | awk -F',' '{gsub(/"/,""); print $6}' \
  | sort -u)

#========================
# 2. 逐 cluster 筛选基因 + 富集
#========================
for cl in $clusters; do
  echo "处理 cluster ${cl}..."

  # 固定文件名，不使用 mktemp，避免报错
  tmpfile="${output}/cluster_${cl}.genelist"

  # 筛选基因并替换 "-" 为 "_"
  awk -F',' -v c="$cl" '
    BEGIN { OFS="\t" }
    NR>1 {
      gsub(/"/, "", $0)
      if ($6 == c && $2 > 1 && $5 < 0.05) {
        gsub("-", "_", $7)
        print $7
      }
    }
  ' "$file" > "$tmpfile"


  # 若无基因，则跳过
  if [ ! -s "$tmpfile" ]; then
    echo "  cluster ${cl} 没有满足条件的基因，跳过富集分析"
    rm "$tmpfile"
    continue
  fi

  # cluster 单独输出目录
  outdir="${output}/cluster_${cl}"
  mkdir -p "$outdir"

  echo "  正在进行富集分析..."
  Rscript "${script}/Enrichment.R" "$anno" "$tmpfile" "$outdir"

done

echo "全部 cluster 处理完毕！"