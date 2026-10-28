#!/bin/bash
#CSUB -J trf
#CSUB -q c01
#CSUB -o trf.out
#CSUB -e trf.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#trf Gifu_v1.0.fasta 2 7 7 80 10 50 500 -d -h

#awk -f trf_dat2bed.awk Gifu_v1.0.fasta.2.7.7.80.10.50.500.dat > Gifu_TRF_all.bed

#awk '($7 >= 10) && ($3-$2 >= 1000)' Gifu_TRF_all.bed > Gifu_TRF_satellite.bed

#sort -k1,1 -k2,2n -k7,7n Gifu_TRF_satellite.bed \
#| awk '!seen[$1":"$2":"$3]++' \
#> Gifu_TRF_minUnit.bed


# 提取 unique unit
#cut -f9 Gifu_TRF_minUnit.bed | sort -u > units.unique.txt

# 转 fasta
#awk '{print ">"$1"\n"$1}' units.unique.txt > units.fa

# cd-hit-est 聚类
#cd-hit-est -i units.fa -o units_clustered -c 0.85 -n 3 -d 0 -T 8 -M 16000

# 生成 unit -> cluster 映射
#awk '
#BEGIN{ OFS="\t" }
#/^>Cluster/{
#    cluster = $2;
#    next;
#}
#{
#    match($0, />[^ ]+/);
#    if (RSTART > 0) {
#        seq = substr($0, RSTART+1, RLENGTH-1);
#        sub(/\.\.\.$/, "", seq);
#        print seq, cluster;
#    }
#}' units_clustered.clstr > unit2cluster.txt

# 加回 TRF 记录
#awk '
#BEGIN{ OFS="\t" }
#NR==FNR { unit2cl[$1]=$2; next }
#{
#    unit = $9;
#    cl = unit2cl[unit];
#    if (cl=="") cl="NA";
#    print $0, cl;
#}
#' unit2cluster.txt Gifu_TRF_minUnit.bed \
#  > Gifu_TRF_satellite_unitCluster.bed



#awk '
#BEGIN{ OFS="\t" }
#{
#    chr  = $1
#    beg  = $2
#    end  = $3
#    cl   = $NF      # 最后一列视为 cluster
#    if (cl == "NA") next    # 不要 NA 的 cluster（如果有）

#    len = end - beg
#    key = chr "\t" cl
#    a[key] += len
#}
#END{
#    for (k in a) {
#        print k, a[k]   # chr cluster total_len
#    }
#}' Gifu_TRF_satellite_unitCluster.bed \
#  > Gifu_cluster_len_per_chr.txt

#sort -k1,1 -k 3,3nr Gifu_cluster_len_per_chr.txt > Gifu_cluster_len_per_chr.sorted.txt


#TOPN=10
#awk -v N=$TOPN '
#BEGIN{ OFS="\t" }
#{
#    chr = $1
#    if (chr != last_chr) {
#        last_chr = chr
#        count = 0
#    }
#    if (count < N) {
#        print
#        count++
#    }
#}' Gifu_cluster_len_per_chr.sorted.txt \
#  > Gifu_top${TOPN}_clusters_per_chr.txt


#cut -f2 Gifu_top${TOPN}_clusters_per_chr.txt \
#  | sort -u \
#  > Gifu_candidate_clusters.txt


#awk '
#BEGIN{ OFS="\t" }
#{
#    chr  = $1
#    beg  = $2
#    end  = $3
#    cl   = $10
#    unit = $9
#    len  = end - beg
#    print chr, beg, end, cl, unit, len
#}' Gifu_TRF_satellite_unitCluster.bed \
#  | sort -k4,4 -k6,6nr \
#  > Gifu_unit_by_cluster_len.txt


#awk '
#BEGIN{ OFS="\t" }
#NR==FNR{
#    cand[$1] = 1;   # 候选 cluster 集合
#    next;
#}
#{
#    cl = $4
#    if (!(cl in cand)) next
#    if (!(cl in seen)) {
#        seen[cl] = 1
#        print
#    }
#}' Gifu_candidate_clusters.txt \
#   Gifu_unit_by_cluster_len.txt \
# > Gifu_cluster_representatives.txt

awk '
BEGIN{ OFS="" }
{
    cl   = $4
    unit = $5
    print ">cluster_", cl
    print unit
}' Gifu_cluster_representatives.txt \
  > blast/Gifu_cluster_reps.fa
