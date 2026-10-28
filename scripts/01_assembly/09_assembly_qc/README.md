# 09 组装质量评估 / Assembly QC

## 中文
**作用**：用回贴深度、BUSCO 与 Merqury 三块独立评估最终组装的完整度与碱基准确性。
**入口**：`01_mapping/{01_ont,02_hifi,03_wgs}/`、`02_busco/`、`03_merqury/{01_hifi,02_wgs,03_pan_wgs}/` 下的脚本
**输入**：组装 `{PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta`；reads ONT `Gifu-jing-ye.pass.ul.fa`、HiFi `Gifu_hifi.fastq.gz`、WGS `Gifu_1-R{1,2}.fastq.gz`；窗口 `Gifu_{100k,10k}.bed`；HiFi k-mer 库 `ccs.merylDB`；BUSCO lineage `embryophyta_odb10`。
**输出**：`gifu_{ont,hifi,wgs}_{100k,10k}_coverage.txt`；BUSCO `Gifu_lotus/`；Merqury `Gifu_hifimer`、`Gifu_wgsmer`、`Gifu_<sample>_wgsmer`；pan-WGS `{Gifu,MG20}_Clean_R{1,2}.fq` 与 `MG20_subsample_{summary,scientific_name}.txt`。
**运行**：进入相应子目录后 `bash <script>`
**工具**：minimap2 `-ax map-ont`/`-ax map-pb`、bwa `mem`；samtools `view -b -q 30`、`sort`、`flagstat`；bedtools `coverage -mean`（100 kb 与 10 kb 窗口）；BBMap `filterbyname.sh`；BLAST+ `blastn -outfmt 5`；BUSCO `-m genome --offline --cpu 64`；Merqury 与 meryl `count k=21`。
**说明**：`01_mapping/01_ont`、`02_hifi`、`03_wgs` 各自完成比对、`-q 30` 过滤、`samtools index` 与 `bedtools coverage`；输入的 BAM 由同一目录脚本产出。

## English
**Purpose**: Assess the finished assembly in three independent blocks: read-remapping depth, BUSCO gene-space completeness and Merqury k-mer QV.
**Entry point**: scripts under `01_mapping/{01_ont,02_hifi,03_wgs}/`, `02_busco/` and `03_merqury/{01_hifi,02_wgs,03_pan_wgs}/`
**Inputs**: assembly `{PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta`; reads ONT `Gifu-jing-ye.pass.ul.fa`, HiFi `Gifu_hifi.fastq.gz`, WGS `Gifu_1-R{1,2}.fastq.gz`; windows `Gifu_{100k,10k}.bed`; HiFi k-mer database `ccs.merylDB`; BUSCO lineage `embryophyta_odb10`.
**Outputs**: `gifu_{ont,hifi,wgs}_{100k,10k}_coverage.txt`; BUSCO `Gifu_lotus/`; Merqury `Gifu_hifimer`, `Gifu_wgsmer`, `Gifu_<sample>_wgsmer`; pan-WGS `{Gifu,MG20}_Clean_R{1,2}.fq` and `MG20_subsample_{summary,scientific_name}.txt`.
**Run**: `cd` into the relevant sub-directory and run `bash <script>`
**Tools**: minimap2 `-ax map-ont`/`-ax map-pb`, bwa `mem`; samtools `view -b -q 30`, `sort`, `flagstat`; bedtools `coverage -mean` (100 kb and 10 kb windows); BBMap `filterbyname.sh`; BLAST+ `blastn -outfmt 5`; BUSCO `-m genome --offline --cpu 64`; Merqury with meryl `count k=21`.
**Notes**: `01_mapping/01_ont`, `02_hifi` and `03_wgs` each cover alignment, `-q 30` filtering, `samtools index` and `bedtools coverage`; the input BAMs are produced by the scripts in the same directory.
