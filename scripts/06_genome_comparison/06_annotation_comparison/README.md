# 06_annotation_comparison / 注释比较

## 中文

**作用**：比较 accession 之间的重复/转座子注释：由 EDTA 的 LTR 插入时间表生成 LTR 插入 BED，再筛出落在特殊区域（着丝粒、rDNA、端粒）之外的插入。
**入口**：`run.sh`
**输入**：`{PROJ_04_LOTUS_GENOME}/03_annotation/MG20/03.edta/MG20_LTR.insertionTime.txt`（EDTA 输出）与 `MG20_repeat.bed`（特殊区域 BED）。
**输出**：`LTR_insertion.bed`、`LTR_in_special_regions.bed`，以及 `LTR_outside.bed`（在 `LTR_in_special_regions.bed` 中无重叠的 LTR 插入）。
**运行**：`bash run.sh`（依次执行 awk 转换、`bedtools intersect -wa` 与 `bedtools intersect -v`）。
**工具**：awk（把 `MG20_LTR.insertionTime.txt` 转成 BED，并修正 start > end 的记录）；bedtools 2.31.1（`intersect -wa`、`intersect -v`）。

## English

**Purpose**: Compare the repeat/transposon annotation between the accessions: an LTR insertion BED is derived from the EDTA insertion-time table and the insertions outside the special regions (centromere, rDNA, telomere) are selected.
**Entry point**: `run.sh`
**Inputs**: `{PROJ_04_LOTUS_GENOME}/03_annotation/MG20/03.edta/MG20_LTR.insertionTime.txt` (EDTA output) and `MG20_repeat.bed` (special-region BED).
**Outputs**: `LTR_insertion.bed`, `LTR_in_special_regions.bed` and `LTR_outside.bed` (LTR insertions with no overlap in `LTR_in_special_regions.bed`).
**Run**: `bash run.sh` (awk conversion, then `bedtools intersect -wa` and `bedtools intersect -v`).
**Tools**: awk (converting `MG20_LTR.insertionTime.txt` into BED, including a swap fix for records where start > end); bedtools 2.31.1 (`intersect -wa`, `intersect -v`).
