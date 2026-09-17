#!/usr/bin/env bash
set -euo pipefail

BLAST="$1"      # trf-blast 输出（至少包含：cluster_id chr start end ...）
FASTA="$2"      # 参考基因组 fasta
OUTDIR="${3:-supercluster_0.7_100kb}"

# ===================== 参数（都从这里控制，python 不写死） =====================
WIN=100000          # bin 大小（100kb）
JAC=0.7             # jaccard 阈值
MERGE_D=50          # merge 间隔（用于计算 super union 区间 & jaccard）
DENS_CUT=0.5        # density 阈值（>0.5 认为高密度）
BIN_KEEP=5          # 一个 super 中 density>DENS_CUT 的 bin 数量必须 > BIN_KEEP 才保留

mkdir -p "$OUTDIR"/{merged,superbeds,superunits,density,tmp}

# 取 blast 中某个 cluster 的所有命中区间（不 merge），输出 BED
blast_units_for_cluster() {
  local cl="$1"
  awk -v C="$cl" 'BEGIN{OFS="\t"}
    $1==C{
      chr=$2; s=$3; e=$4;
      if(s>e){t=s;s=e;e=t}
      s0=s-1; if(s0<0)s0=0;
      print chr,s0,e
    }' "$BLAST"
}

echo "[1/7] Prepare genome sizes & windows (${WIN}bp)"
samtools faidx "$FASTA"
cut -f1,2 "${FASTA}.fai" > "$OUTDIR/genome.sizes"
bedtools makewindows -g "$OUTDIR/genome.sizes" -w "$WIN" \
  | sort -k1,1 -k2,2n \
  > "$OUTDIR/windows.${WIN}.bed"

echo "[2/7] Split BLAST by cluster and merge hits (only for Jaccard/super union)"
cut -f1 "$BLAST" | sort -u > "$OUTDIR/clusters.list"

while read -r cl; do
  blast_units_for_cluster "$cl" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/merged/${cl}.bed"
done < "$OUTDIR/clusters.list"

echo "[3/7] Compute pairwise Jaccard and keep >= ${JAC}"
: > "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"
mapfile -t CLS < "$OUTDIR/clusters.list"
N=${#CLS[@]}

for ((i=0;i<N;i++)); do
  for ((j=i+1;j<N;j++)); do
    A="${CLS[i]}"
    B="${CLS[j]}"
    J=$(bedtools jaccard -a "$OUTDIR/merged/${A}.bed" -b "$OUTDIR/merged/${B}.bed" \
        | awk 'NR==2{print $3}')
    awk -v a="$A" -v b="$B" -v j="$J" -v cut="$JAC" 'BEGIN{
      if((j+0) >= cut) print a"\t"b"\t"j
    }' >> "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"
  done
done

echo "[4/7] Union-Find to build superclusters"
python build_supercluster_unionfind.py \
  --outdir "$OUTDIR" \
  --jac "$JAC" \
  --pairs "$OUTDIR/pairs.jaccard_ge${JAC}.tsv"

echo "[5/7] Build supercluster union BEDs (merge across member clusters)"
cut -f2 "$OUTDIR/cluster2super.tsv" | sort -u > "$OUTDIR/supers.list"

while read -r sup; do
  awk -v S="$sup" '$2==S{print $1}' "$OUTDIR/cluster2super.tsv" \
  | while read -r cl; do
      cat "$OUTDIR/merged/${cl}.bed"
    done \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/superbeds/${sup}.bed"
done < "$OUTDIR/supers.list"

echo "[6/7] Compute density per bin (col5 = density)"
# 输出列：chr start end group densitybin_len
calc_one () {
  local NAME="$1"
  local BED="$2"
  bedtools coverage -a "$OUTDIR/windows.${WIN}.bed" -b "$BED" \
    | awk -v G="$NAME" 'BEGIN{OFS="\t"}
        {
          # bedtools coverage 默认输出：
          # $1 chr, $2 start, $3 end,
          # $4 numFeaturesInA, $5 lenA,
          # $6 basesCovered, $7 fractionCovered
          density=$7
          covered_bp=$6
          bin_len=$5
          print $1,$2,$3,G,density,covered_bp,bin_len
        }'
}


# ALL：所有 super union 再合并
cat "$OUTDIR"/superbeds/super_*.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d "$MERGE_D" \
  > "$OUTDIR/superbeds/ALL.bed"

calc_one "ALL" "$OUTDIR/superbeds/ALL.bed" > "$OUTDIR/density/density.ALL.tsv"

: > "$OUTDIR/density/density.supers.tsv"
while read -r sup; do
  calc_one "$sup" "$OUTDIR/superbeds/${sup}.bed" >> "$OUTDIR/density/density.supers.tsv"
done < "$OUTDIR/supers.list"

cat "$OUTDIR/density/density.ALL.tsv" "$OUTDIR/density/density.supers.tsv" \
  > "$OUTDIR/density/density.long.tsv"

echo "[7/7] Filter supers; export unit beds from BLAST"
# 统计每个 super：density>DENS_CUT 的 bin 数
awk -v cut="$DENS_CUT" '
  $4!="ALL" && $5>cut {cnt[$4]++}
  END{
    for(g in cnt){
      if(cnt[g] > '"$BIN_KEEP"') print g"\t"cnt[g]
    }
  }' "$OUTDIR/density/density.supers.tsv" \
  | sort -k1,1 \
  > "$OUTDIR/kept_supers.tsv"

cut -f1 "$OUTDIR/kept_supers.tsv" > "$OUTDIR/kept_supers.list" || true

KEPT_N=$(wc -l < "$OUTDIR/kept_supers.list" 2>/dev/null || echo 0)
echo "[INFO] Kept supers: ${KEPT_N} (criterion: bins(density>${DENS_CUT}) > ${BIN_KEEP})"

if [[ "${KEPT_N}" -eq 0 ]]; then
  echo "[WARN] No super passed filter."
  echo "Done."
  exit 0
fi

# 导出 unit bed：从 BLAST 中提取（不使用 merge 后结果）
while read -r sup; do
  out="$OUTDIR/superunits/${sup}.units.bed"
  : > "$out"
  awk -v S="$sup" '$2==S{print $1}' "$OUTDIR/cluster2super.tsv" \
    | while read -r cl; do
        blast_units_for_cluster "$cl" >> "$out"
      done
  sort -k1,1 -k2,2n "$out" -o "$out"
done < "$OUTDIR/kept_supers.list"

echo "Done."
echo "  Mapping:       $OUTDIR/cluster2super.tsv"
echo "  Density long:  $OUTDIR/density/density.long.tsv"
echo "  Kept list:     $OUTDIR/kept_supers.list"
echo "  Units bed:     $OUTDIR/superunits/*.units.bed"
