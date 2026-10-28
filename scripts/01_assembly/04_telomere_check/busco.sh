#4.1 busco初步检查组装完整性

busco -i contigs.fa \
      -l {SOFTWARE}/busco_downloads/lineages/embryophyta_odb10 \
      -o Gifu_lotus \
      -m genome \
      --cpu 64