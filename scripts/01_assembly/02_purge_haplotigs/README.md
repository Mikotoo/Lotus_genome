# 02 去单倍型冗余 / Purge haplotigs

## 中文
**作用**：把 HiFi reads 回贴到 hifiasm primary contig 统计读深与覆盖度，并去除单倍型冗余（haplotig）contig。
**入口**：`02.1.minimap.sh`（比对）→ `02.2.purge.sh`（去冗余）
**输入**：`Gifu.hifiasm.fasta`（步骤 01 输出）；HiFi reads。
**输出**：`aligned.bam`、`aligned.bam.200.gencov`、`coverage_stats.csv`；去冗余后的 contig FASTA 与 haplotig FASTA（`purge_haplotigs purge` 的默认名）。
**运行**：`bash 02.1.minimap.sh`，然后 `bash 02.2.purge.sh`
**工具**：minimap2 `-ax map-pb`；samtools `view -hF 256`、`sort -@ 20 -m 1G -T tmp.ali`；purge_haplotigs `readhist -t 88`、`contigcov -l 15 -m 100 -h 185`、`purge -t 88 -a 50`。

## English
**Purpose**: Align HiFi reads back to the hifiasm primary contigs to measure read depth and coverage, then purge haplotig (duplicated) contigs.
**Entry point**: `02.1.minimap.sh` (alignment) → `02.2.purge.sh` (purging)
**Inputs**: `Gifu.hifiasm.fasta` (step 01 output); HiFi reads.
**Outputs**: `aligned.bam`, `aligned.bam.200.gencov`, `coverage_stats.csv`; the purged contig FASTA and the haplotig FASTA (the `purge_haplotigs purge` default names).
**Run**: `bash 02.1.minimap.sh`, then `bash 02.2.purge.sh`
**Tools**: minimap2 `-ax map-pb`; samtools `view -hF 256`, `sort -@ 20 -m 1G -T tmp.ali`; purge_haplotigs `readhist -t 88`, `contigcov -l 15 -m 100 -h 185`, `purge -t 88 -a 50`.
