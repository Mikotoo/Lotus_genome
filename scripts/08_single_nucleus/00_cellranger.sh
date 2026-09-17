#!/bin/bash
#CSUB -J cellranger
#CSUB -q c01
#CSUB -o cellranger.out
#CSUB -e cellranger.error
#CSUB -n 64
#CSUB -R span[hosts=1]

REF={PROJ_LOTUS_HX}/07_singleCell/00_ref/fasta/genome.fa
DATADIR={HOME_HX}/Single_cell_data
SAMPLELIST=sample.list

# module load cellranger/7.2.0   # 如果集群环境要求，请加载软件

while read ProjectID CNName RunID SampleName
do
    echo ">>> Processing $SampleName ($RunID)"

    R1=${DATADIR}/${ProjectID}/${RunID}/${RunID}_S1_L001_R1_001.fastq.gz
    R2=${DATADIR}/${ProjectID}/${RunID}/${RunID}_S1_L001_R2_001.fastq.gz

    if [[ ! -f $R1 ]]; then
        echo "❌ Missing R1 file: $R1"
        continue
    fi

    outdir=./01_cellranger_out/${SampleName}
    mkdir -p ${outdir}

    cellranger count \
        --id=${SampleName} \
        --transcriptome={PROJ_LOTUS_HX}/07_singleCell/00_ref \
        --fastqs=$(dirname $R1) \
        --sample=${RunID} \
        --create-bam=true \
        --localcores=32 \
        --localmem=120 \
        --chemistry=SC3Pv3 \
        --output-dir=${outdir}

done < ${SAMPLELIST}
