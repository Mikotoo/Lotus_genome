#3.1 排除长度较长（10Mb）的染色体序列，提取其余序列（可能是质体），比对到nt库（需要从ncbi下载到本地），输出fmt5格式

blastn -query contigs.fa -out nt.blast.xml -db {SOFTWARE_ZC}/nt_lib/nt -outfmt 5 -evalue 1e-5 -num_threads 88

python 03.2.get_gi.py
python 03.3.get_name.py