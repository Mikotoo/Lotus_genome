# 01 hifiasm 组装 / hifiasm assembly

## 中文
**作用**：用 hifiasm 对 ONT ultra-long reads 与 PacBio HiFi reads 进行初步 contig 组装，并把 primary contig GFA 转成 FASTA 并建索引。
**入口**：`01.hifiasm.sh`
**输入**：ONT UL reads; HiFi reads
**输出**：`Gifu.sam.bp.p_ctg.gfa`、`Gifu.hifiasm.fasta`、`Gifu.hifiasm.fasta.fai`
**运行**：`bash 01.hifiasm.sh`
**工具**：hifiasm `-t 64 -o Gifu.sam --telo-m CCCTAAA --ul`（UL reads 在前、HiFi FASTQ 在后）；gfatools `gfa2fa`；samtools `faidx`

## English
**Purpose**: Assemble the Gifu accession from ONT ultra-long and PacBio HiFi reads with hifiasm, then convert the primary contig GFA to FASTA and index it.
**Entry point**: `01.hifiasm.sh`
**Inputs**: ONT UL reads; HiFi reads
**Outputs**: `Gifu.sam.bp.p_ctg.gfa`, `Gifu.hifiasm.fasta`, `Gifu.hifiasm.fasta.fai`.
**Run**: `bash 01.hifiasm.sh`
**Tools**: hifiasm `-t 64 -o Gifu.sam --telo-m CCCTAAA --ul` (UL reads first, then the HiFi FASTQ); gfatools `gfa2fa`; samtools `faidx`.
