#!/bin/bash
#CSUB -J bwa
#CSUB -q c01
#CSUB -o bwa.out
#CSUB -e bwa.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#samtools view -b -f 4 01.align/GifuT2T.sorted.bam > Gifu_unmapped.bam

#samtools fastq -f 4 -1 Gifu_unmapped_1.fq -2 Gifu_unmapped_2.fq -0 /dev/null -s /dev/null 01.align/GifuT2T.sorted.bam
#samtools fastq -f 4 -1 MG20_unmapped_1.fq -2 MG20_unmapped_2.fq -0 /dev/null -s /dev/null 01.align/MG20T2T.sorted.bam

#cat Gifu_unmapped_1.fq Gifu_unmapped_2.fq \
# | awk 'NR%4==1{ id=$1; sub(/^@/,"",id); print id }' \
# | sort -u > Gifu_unmapped.ids


#cat MG20_unmapped_1.fq MG20_unmapped_2.fq \
# | awk 'NR%4==1{ id=$1; sub(/^@/,"",id); print id }' \
# | sort -u > MG20_unmapped.ids

#awk '{print $1"/1"; print $1"/2"}' MG20_unmapped.ids \
#  | sort -u > MG20_unmapped.literal.ids

#awk '{print $1"/1"; print $1"/2"}' Gifu_unmapped.ids \
#  | sort -u > Gifu_unmapped.literal.ids




filterbyname.sh \
  in={PROJ_04_LOTUS_GENOME}/01_data/Gifu/WGS/Gifu_Clean_R1.fq in2={PROJ_04_LOTUS_GENOME}/01_data/Gifu/WGS/Gifu_Clean_R2.fq \
  out=Gifu_Clean_R1.fq out2=Gifu_Clean_R2.fq \
  names=Gifu_unmapped.literal.ids include=f ow=t

filterbyname.sh \
  in={PROJ_04_LOTUS_GENOME}/01_data/MG20/WGS/MG20_Clean_R1.fq in2={PROJ_04_LOTUS_GENOME}/01_data/MG20/WGS/MG20_Clean_R2.fq \
  out=MG20_Clean_R1.fq out2=MG20_Clean_R2.fq \
  names=MG20_unmapped.literal.ids include=f ow=t