#!/bin/bash
#CSUB -J hisat_07
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# hisat2 比对 + stringtie 定量 (参考: 01_ChiDou_genome/08_RNAseq/02.ChiDou_tissue/align.sh)

export PATH={CONDA_PY3}/bin:$PATH

INDEX="{HOME_ZC}/Project/Ref/Medicago_v4/Medic_v4"
REF="{HOME_ZC}/Project/Ref/Medicago_v4/Medic_truncatula.gtf"
THREADS=32

mkdir -p align count
: > prepDE.list

while read -r sample fq1 fq2; do
    [[ -z "${sample}" || "${sample}" == "#"* ]] && continue

    echo "[START] ${sample}"

    # ---- hisat2 比对 ----
    if [[ -n "${fq2}" ]]; then
        # 双端
        hisat2 -p "${THREADS}" -x "${INDEX}" -1 "${fq1}" -2 "${fq2}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1 且正确配对 YT:Z:CP)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" | grep "YT:Z:CP" > "${sample}_unique.sam"
    else
        # 单端
        hisat2 -p "${THREADS}" -x "${INDEX}" -U "${fq1}" \
            -S "${sample}.sam" --summary-file "${sample}_hisat.log" --new-summary \
            || { echo "ERROR: hisat2 failed for ${sample}"; continue; }
        # unique 比对 (NH:i:1)
        grep -E "NH:i:1([[:space:]]|$)" "${sample}.sam" > "${sample}_unique.sam"
    fi

    # ---- unique -> 排序 bam ----
    samtools view -H "${sample}.sam" > "${sample}_head.sam"
    cat "${sample}_unique.sam" >> "${sample}_head.sam"
    samtools sort -@ "${THREADS}" -o "${sample}.bam" "${sample}_head.sam"
    rm -f "${sample}.sam" "${sample}_unique.sam" "${sample}_head.sam"
    samtools index "${sample}.bam"

    # ---- stringtie 定量 ----
    mkdir -p "align/${sample}"
    stringtie -p "${THREADS}" -G "${REF}" -e -B \
        -o "align/${sample}/${sample}.gtf" \
        -A "align/${sample}/${sample}.tsv" \
        "${sample}.bam" \
        || { echo "ERROR: stringtie failed for ${sample}"; continue; }

    printf '%s\t%s\n' "${sample}" "align/${sample}/${sample}.gtf" >> prepDE.list
    echo "[DONE] ${sample}"
done < samples.conf

# ---- prepDE 汇总 count 矩阵 ----
python {CONDA_PY3}/bin/prepDE.py \
    -i prepDE.list \
    -g count/prepDE_counts.csv \
    -t count/prepDE_transcript.csv

echo "[ALL DONE] count matrix written to count/"
