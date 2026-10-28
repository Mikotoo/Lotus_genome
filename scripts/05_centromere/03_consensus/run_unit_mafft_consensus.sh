#!/usr/bin/env bash
set -euo pipefail

# ========= 需要你改的 =========
GIFU_FA="Lotus_GifuT2T_v1.0.fasta"
MG20_FA="Lotus_MG20T2T_v1.0.fasta"

UNITS=(unit_161 unit_172 unit_186 unit_330)
GIFU_BED_PREFIX="gifu_"
MG20_BED_PREFIX="mg20_"

THREADS=32
OUTDIR="separate_then_merge_consensus"

CONS_PY="get_consensus.py"   # 你给的共识脚本文件名
THRESH=0.5
USE_IUPAC=""   # 不用IUPAC就置空：USE_IUPAC=""
INCLUDE_GAP="--include-gap"            # 要算gap就：INCLUDE_GAP="--include-gap"
# =================================

mkdir -p "$OUTDIR"/{01_seqs,02_msa,03_consensus,03_dist,04_merge2cons_msa,05_final}

command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found"; exit 1; }
command -v mafft    >/dev/null 2>&1 || { echo "[ERROR] mafft not found"; exit 1; }
[[ -f "$CONS_PY" ]] || { echo "[ERROR] missing $CONS_PY"; exit 1; }

extract_and_consensus () {
  local genome="$1"
  local bed="$2"
  local tag="$3"     # Gifu or MG20
  local unit="$4"

  local seqfa="$OUTDIR/01_seqs/${tag}_${unit}.fa"
  local alnfa="$OUTDIR/02_msa/${tag}_${unit}.aln.fa"
  local consfa="$OUTDIR/03_consensus/${tag}_${unit}.cons.fa"
  local disttsv="$OUTDIR/03_dist/${tag}_${unit}.dist.tsv"

  if [[ ! -s "$bed" ]]; then
    echo "[WARN] $tag $unit bed missing/empty: $bed"
    return 1
  fi

  # 提取序列（按strand统一方向）
  bedtools getfasta -fi "$genome" -bed "$bed" -fo "$seqfa" -name -s

  # header 加前缀避免冲突（可选但建议）
  awk -v p="${tag}|" '/^>/{sub(/^>/,">"p);} {print}' "$seqfa" > "$seqfa.tmp" && mv "$seqfa.tmp" "$seqfa"

  local nseq
  nseq=$(grep -c "^>" "$seqfa" || true)
  if [[ "$nseq" -lt 2 ]]; then
    echo "[WARN] $tag $unit sequences <2 ($nseq), skip"
    return 1
  fi

  # MAFFT
  mafft --thread "$THREADS" --auto "$seqfa" > "$alnfa"

  # 共识 + 距离
  python "$CONS_PY" \
    -i "$alnfa" \
    -o "$consfa" \
    -d "$disttsv" \
    -t "$THRESH" \
    $USE_IUPAC \
    $INCLUDE_GAP \
    --name "${tag}_${unit}_consensus"

  echo "[OK] $tag $unit done: nseq=$nseq"
  return 0
}

for u in "${UNITS[@]}"; do
  echo "==== $u ===="

  gbed="${GIFU_BED_PREFIX}${u}.bed"
  mbed="${MG20_BED_PREFIX}${u}.bed"

  g_ok=0
  m_ok=0
  extract_and_consensus "$GIFU_FA" "$gbed" "Gifu" "$u" && g_ok=1 || true
  extract_and_consensus "$MG20_FA" "$mbed" "MG20" "$u" && m_ok=1 || true

  # 两个基因组都成功时：把两个共识再对齐，得到跨基因组比较/最终合并共识
  if [[ "$g_ok" -eq 1 && "$m_ok" -eq 1 ]]; then
    two_cons="$OUTDIR/04_merge2cons_msa/${u}.Gifu_MG20.cons.fa"
    cat "$OUTDIR/03_consensus/Gifu_${u}.cons.fa" \
        "$OUTDIR/03_consensus/MG20_${u}.cons.fa" > "$two_cons"

    two_aln="$OUTDIR/04_merge2cons_msa/${u}.Gifu_MG20.cons.aln.fa"
    mafft --thread "$THREADS" --auto "$two_cons" > "$two_aln"

    # 可选：再生成“最终跨基因组共识”（其实就是两条共识的共识）
    final_cons="$OUTDIR/05_final/${u}.final.cons.fa"
    final_dist="$OUTDIR/05_final/${u}.final.dist.tsv"
    python "$CONS_PY" \
      -i "$two_aln" \
      -o "$final_cons" \
      -d "$final_dist" \
      -t 0.5 \
      $USE_IUPAC \
      $INCLUDE_GAP \
      --name "${u}_Final_GifuMG20"

    echo "[OK] $u merged consensus done"
  else
    echo "[WARN] $u cannot merge 2-consensus (need both genomes succeed)"
  fi
done

echo "[DONE] outputs in: $OUTDIR"
