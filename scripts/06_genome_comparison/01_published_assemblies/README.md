# 01_published_assemblies / 已发表组装比较

## 中文

**作用**：用 MUMmer 点阵图比较已发表的 Gifu 与 MG20 组装；gap 内容比较见 `gap_analysis/`。
**入口**：`mummer.sh`
**输入**：`ref.filter.delta`；`-R {PROJ_LOTUS_ZC}/Gifu/05_3.mummer/Gifu_ref.fa`、`-Q {PROJ_LOTUS_ZC}/MG20/05_3.mummer/MG20_ref.fa`。
**输出**：`ref.filter.delta`、`ref.filter.delta.coords` 与前缀 `ref.filter` 的 `mummerplot` PNG。
**运行**：`bash mummer.sh`。
**工具**：MUMmer4 4.0.0rc1 — `nucmer --prefix ref`、`delta-filter -i 89 -l 1000 -1`、`show-coords -THrd`、`mummerplot --layout --png --large`；CSUB（`-q c01`、`-n 64`、`span[hosts=1]`）。

## English

**Purpose**: Dot-plot comparison of the published Gifu and MG20 assemblies; the gap content comparison is in `gap_analysis/`.
**Entry point**: `mummer.sh`
**Inputs**: `ref.filter.delta`; `-R {PROJ_LOTUS_ZC}/Gifu/05_3.mummer/Gifu_ref.fa`, `-Q {PROJ_LOTUS_ZC}/MG20/05_3.mummer/MG20_ref.fa`.
**Outputs**: `ref.filter.delta`, `ref.filter.delta.coords` and the `mummerplot` PNG for prefix `ref.filter`.
**Run**: `bash mummer.sh`.
**Tools**: MUMmer4 4.0.0rc1 — `nucmer --prefix ref`, `delta-filter -i 89 -l 1000 -1`, `show-coords -THrd`, `mummerplot --layout --png --large`; CSUB (`-q c01`, `-n 64`, `span[hosts=1]`).
