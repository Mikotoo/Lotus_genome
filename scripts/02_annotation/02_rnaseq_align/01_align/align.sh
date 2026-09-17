#!/bin/bash
#CSUB -J hisat
#CSUB -q c01
#CSUB -o hisat.out
#CSUB -e hisat.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#cd {PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/00.index
#hisat2-build -p 64 Gifu_v1.0.fasta Gifu
#index={PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq/00.index/Gifu
dir={PROJ_LOTUS_ZC}/Gifu/09.annotation/04.rnasq
ref1={PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Lotus_GifuT2T.gtf
ref2={PROJ_LOTUS_ZC}/Gifu/09.annotation/06.braker/Gifu_filtered.gtf
cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/align.conf | while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

cd $dir
#hisat2 -p 64 -x $index -1 $fq1 -2 $fq2 -S $sample.sam --summary-file $sample_hisat.log --new-summary
#grep "NH:i:1" $sample.sam | grep "YT:Z:CP" > ${sample}_unique.sam
#samtools view -H $sample.sam > ${sample}_head.sam
#cat ${sample}_unique.sam >> ${sample}_head.sam
#samtools sort -@ 64 -o $sample.bam ${sample}_head.sam
#rm *sam
#samtools index $sample.bam
mkdir 02.stringtie_geta
mkdir 02.stringtie_geta/$sample && cd 02.stringtie_geta/$sample
stringtie -p 64 -G $ref1 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam

cd $dir
mkdir 03.stringtie_braker
mkdir 03.stringtie_braker/$sample && cd 03.stringtie_braker/$sample
stringtie -p 64 -G $ref2 -e -B -o $sample.gtf -A $sample.tsv $dir/01.align/$sample.bam

cd $dir
mkdir 02.tpm_geta
sed '1d' 02.stringtie_geta/$sample/$sample.tsv | awk '{print $1"\t"$9}'|sort > 02.tpm_geta/$sample.tpm

cd $dir
mkdir 03.tpm_braker
sed '1d' 03.stringtie_braker/$sample/$sample.tsv | awk '{print $1"\t"$9}'|sort > 03.tpm_braker/$sample.tpm

done

cd $dir/02.tpm_geta
ls *.tpm | sort > tpm.files
mapfile -t files < tpm.files
cp "${files[0]}" __merged.tmp
for f in "${files[@]:1}"; do
  join -t $'\t' -a1 -a2 -e 0 -o auto __merged.tmp "$f" > __merged.next
  mv __merged.next __merged.tmp
done
header="gene_id"
for f in "${files[@]}"; do header+="\t${f%.tpm}"; done
{ echo -e "$header"; cat __merged.tmp; } > all_samples.tpm.tsv

rm -f __merged.tmp tpm.files


cd $dir/03.tpm_braker
ls *.tpm | sort > tpm.files
mapfile -t files < tpm.files
cp "${files[0]}" __merged.tmp
for f in "${files[@]:1}"; do
  join -t $'\t' -a1 -a2 -e 0 -o auto __merged.tmp "$f" > __merged.next
  mv __merged.next __merged.tmp
done
header="gene_id"
for f in "${files[@]}"; do header+="\t${f%.tpm}"; done
{ echo -e "$header"; cat __merged.tmp; } > all_samples.tpm.tsv

rm -f __merged.tmp tpm.files