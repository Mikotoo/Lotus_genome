# 08 抛光 / Polishing

## 中文
**作用**：对补洞后的组装做长读长与 k-mer 迭代抛光：Winnowmap2 比对、falconc 过滤、Racon 共识、Merfin/bcftools 校正。
**入口**：`bash.sh`（调用第三方 `polish.sh`；`build_merylDB.sh` 构建 k-mer 库）
**输入**：`ontPolish.iter_2.consensus.fasta`（ONT 抛光后的草图）；HiFi reads `{PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz`；读 k-mer 库 `ccs.merylDB`。
**输出**：`ccs.merylDB`；每轮 `<prefix>.iter_N.*`（`winnowmap.bam`、`falconc.sam`、`racon.fasta`、merfin 校正 VCF、`consensus.fasta`）与 `.memtime` 日志。
**运行**：`bash bash.sh`
**工具**：`polish.sh 88 2 ontPolish.iter_2.consensus.fasta <hifi.fastq.gz> ccs.merylDB hifiPolish`；meryl `count k=21`；winnowmap `-k 21 -ax map-pb`；falconc `bam-filter-clipped -F 0x104`；racon `-L -S`；merfin `-polish -peak 106.7`；bcftools `consensus -H 1`。
**说明**：`bash.sh` 先跑前缀 `ontPolish` 的 ONT 抛光，再以该输出为草图跑 HiFi 抛光；`polish.sh` 为 T2T-Polish 的公开版本（public domain），来源与许可见 `third_party/THIRD_PARTY_NOTICES.md`。

## English
**Purpose**: Iteratively polish the gap-closed assembly with long reads and k-mers: Winnowmap2 mapping, falconc filtering, Racon consensus and Merfin/bcftools correction.
**Entry point**: `bash.sh` (calls the third-party `polish.sh`; `build_merylDB.sh` builds the k-mer database)
**Inputs**: `ontPolish.iter_2.consensus.fasta` (the draft after ONT polishing); HiFi reads `{PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz`; read k-mer database `ccs.merylDB`.
**Outputs**: `ccs.merylDB`; per iteration `<prefix>.iter_N.*` (`winnowmap.bam`, `falconc.sam`, `racon.fasta`, the Merfin correction VCF, `consensus.fasta`) and `.memtime` logs.
**Run**: `bash bash.sh`
**Tools**: `polish.sh 88 2 ontPolish.iter_2.consensus.fasta <hifi.fastq.gz> ccs.merylDB hifiPolish`; meryl `count k=21`; winnowmap `-k 21 -ax map-pb`; falconc `bam-filter-clipped -F 0x104`; racon `-L -S`; merfin `-polish -peak 106.7`; bcftools `consensus -H 1`.
**Notes**: `bash.sh` first runs the ONT polishing with prefix `ontPolish`, then uses its output as the draft for the HiFi run; `polish.sh` is the public T2T-Polish script (public domain), see `third_party/THIRD_PARTY_NOTICES.md` for its source and licence.
