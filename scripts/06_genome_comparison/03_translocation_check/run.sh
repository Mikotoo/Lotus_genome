#!/bin/bash
#CSUB -J run
#CSUB -q c01
#CSUB -o run.out
#CSUB -e run.error
#CSUB -n 88
#CSUB -R span[hosts=1]

# Two entry points share this directory:
#   (a) batch structural typing of the Chr1<->Chr2 translocation across samples
#       (the `python batch_translocation_typing.py` block below);
#   (b) the per-pair counting commands used for the per-pair plots, recorded
#       at the end of the file. They were run interactively rather than from
#       this script, which is why they are commented out.

echo "[*] STEP 1: 从 trans.conf 提取 BAM 列表 ..."
cut -f2,3 trans.conf | tr '\t' '\n' | sort -u > bam.list

echo "[*] STEP 2: 自动检查/补齐 BAM 索引 ..."
while read bam; do
    if [ ! -f "$bam" ]; then
        echo "[WARN] 找不到 BAM: $bam"
        continue
    fi

    bai1="${bam}.bai"
    bai2="${bam%.bam}.bai"

    if [ -f "$bai1" ] || [ -f "$bai2" ]; then
        echo "[OK] 已存在索引: $bam"
    else
        echo "[INDEX] 开始建立索引: $bam"
        samtools index "$bam"
    fi
done < bam.list


python batch_translocation_typing.py \
  --conf trans.conf \
  --region-gifu Chr1:22829936-22849936 \
  --region-mg   Chr2:99260160-99280160 \
  --find-script find_cross_mapped_reads.py \
  --out-summary lotus_translocation_typing.tsv


# ── Per-pair counting for the per-pair plots (run interactively) ─────────
#samtools index gifu2gifu.bam
#samtools index gifu2mg.bam
#samtools index mg2gifu.bam
#samtools index mg2mg.bam

#python count_pairs_chr_chr2.py \
#  --bam gifu2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > gifu2gifu_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam mg2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > mg2gifu_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam gifu2mg.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > gifu2mg_chr1_chr2_pairs.txt

#python count_pairs_chr_chr2.py \
#  --bam mg2mg.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  > mg2mg_chr1_chr2_pairs.txt

#python count_pairs_region1_chr1tail_chr2.py \
#  --bam gifu2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  --bp 22839936 \
#  > gifu2gifu_region_pairs.txt

#python count_pairs_region1_chr1tail_chr2.py \
#  --bam mg2gifu.bam \
#  --chr1 Chr1 \
#  --chr2 Chr2 \
#  --bp 22839936 \
#  > mg2gifu_region_pairs.txt

#python find_cross_mapped_reads.py \
#  --bam-gifu gifu2gifu.bam \
#  --bam-mg   gifu2mg.bam \
#  --region-gifu Chr1:22829936-22849936 \
#  --region-mg   Chr2:99260160-99280160 \
#  --out-prefix gifu_cross

#python find_cross_mapped_reads.py \
#  --bam-gifu mg2gifu.bam \
#  --bam-mg   mg2mg.bam \
#  --region-gifu Chr1:22829936-22849936 \
#  --region-mg   Chr2:99260160-99280160 \
#  --out-prefix mg_cross
