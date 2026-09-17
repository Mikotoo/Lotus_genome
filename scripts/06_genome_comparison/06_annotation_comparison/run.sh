#awk 'BEGIN{OFS="\t"}
#{
#  split($1,a,":");
#  split(a[2],b,"\\.\\.");
#  chr=a[1];
#  start=b[1]-1;
#  end=b[2];
#  if(start<0) start=0;
#  print chr,start,end,$2,$3
#}' {PROJ_04_LOTUS_GENOME}/03_annotation/MG20/03.edta/MG20_LTR.insertionTime.txt \
#> LTR_insertion.bed
#cat LTR_insertion.bed |awk '$2<$3{print $0}' > tmp
#cat LTR_insertion.bed |awk '$2>$3{print $1"\t"$3"\t"$2"\t"$4"\t"$5}' > tmp1
#cat tmp tmp1 |sort -k 1,1 -k 2n > LTR_insertion.bed


#bedtools intersect \
#  -a LTR_insertion.bed \
#  -b MG20_repeat.bed \
#  -wa \
#> LTR_in_special_regions.bed


bedtools intersect -v \
  -a LTR_insertion.bed \
  -b LTR_in_special_regions.bed \
> LTR_outside.bed
