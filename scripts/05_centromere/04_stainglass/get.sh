cat mg20_cent.fa.fai|awk '{print $1}'|while read id;
do 
cd {PROJ_04_LOTUS_GENOME}/08.cent/stainGlass
mkdir ${id} && cd ${id}
echo ${id} > ${id}
perl ~/script/getseq.pl ${id} ../mg20_cent.fa > ${id}.fa
samtools faidx ${id}.fa
bedtools makewindows -g ${id}.fa.fai -w 2000 > ${id}_2k.bed
bedtools getfasta -fi ${id}.fa -bed ${id}_2k.bed > ${id}_2k.fasta
done