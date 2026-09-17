#!/bin/bash
set -euo pipefail

cul="$(pwd)"
runlist="$cul/scTenifoldKnk_run.list"
genelist=$1
celltype=$2

# 清空 run list
: > "$runlist"

while IFS=$'\t ' read -r gene group || [ -n "$gene" ]; do
  # 去掉 Windows CR
  gene="${gene//$'\r'/}"
  group="${group//$'\r'/}"

  # 跳过非法行
  [ -z "$gene" ] && continue
  [ -z "$group" ] && continue

  # 构建目录：06_batch_multiple_scTenifoldKnk / group / gene
  base_dir="$cul/06_batch_multiple_scTenifoldKnk/$group/$gene"
  mkdir -p "$base_dir"

  job="$base_dir/${gene}_scTenifoldKnk.sh"

  cat > "$job" <<EOF
#!/bin/bash
#CSUB -J ${gene}_Knk
#CSUB -q c02
#CSUB -o $base_dir/${gene}_Knk.out
#CSUB -e $base_dir/${gene}_Knk.error
#CSUB -n 16
#CSUB -R span[hosts=1]

set -eo pipefail   # conda 前不要 -u

source ~/anaconda3/etc/profile.d/conda.sh
conda activate singlecell

set -euo pipefail  # 激活后再启用 -u

# 进入 gene 专属目录
cd "$base_dir"

# 软链接数据（避免重复拷贝）
ln -sf $cul/02_processData/zz_00_pea.rds
ln -sf $cul/All_species_merged_by_GifuT2T.SNF_Symbol.tsv

# 执行 scTenifoldKnk
Rscript $cul/Script/07_batch_multiple_scTenifoldKnk.R "$gene" "$celltype"
EOF

  chmod +x "$job"
  echo "$job" >> "$runlist"

done < "$genelist"

echo "Generated run list: $runlist"
