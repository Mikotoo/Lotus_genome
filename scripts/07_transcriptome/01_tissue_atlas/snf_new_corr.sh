#!/bin/bash
set -euo pipefail

if [ $# -lt 3 ]; then
  echo "Usage: bash $0 expr.tsv genes.tsv outprefix"
  exit 1
fi

EXPR="$1"
GENES="$2"
OUT="$3"

ODIR="${OUT}.snf_new_corr"
mkdir -p "$ODIR"

# =========================
# 1) Python: corr + top5
# =========================
python3 - "$EXPR" "$GENES" "$OUT" <<'PY'
import sys, os
import pandas as pd
import numpy as np
from scipy.stats import t as tdist

expr, genes, outp = sys.argv[1:4]
odir = f"{outp}.snf_new_corr"
os.makedirs(odir, exist_ok=True)

df = pd.read_csv(expr, sep="\t", header=0)
df = df.rename(columns={df.columns[0]:"GeneID"})
df = df.drop_duplicates("GeneID").set_index("GeneID")

meta = pd.read_csv(genes, sep="\t", header=None)
if meta.shape[1] < 2:
    raise SystemExit("gene.list 至少2列：GeneID 与 Category")
meta.columns = ["GeneID","Category"] + [f"c{i}" for i in range(3, meta.shape[1]+1)]
meta = meta.drop_duplicates("GeneID").set_index("GeneID")

snf_ids = meta.index[meta["Category"].astype(str).str.contains("SNF", case=False, na=False)].tolist()
new_ids = meta.index[meta["Category"].astype(str).str.contains("New", case=False, na=False)].tolist()

snf_ids = [g for g in snf_ids if g in df.index]
new_ids = [g for g in new_ids if g in df.index]

if len(snf_ids) == 0 or len(new_ids) == 0:
    raise SystemExit(f"在表达矩阵中匹配到 SNF={len(snf_ids)}, New={len(new_ids)}；检查 GeneID 是否一致")

X = df.loc[snf_ids].to_numpy(float)
Y = df.loc[new_ids].to_numpy(float)
n = X.shape[1]
dfree = n - 2

def zscore(A):
    mu = A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = np.nan
    return (A - mu) / sd

Xz = zscore(X)
Yz = zscore(Y)

R = (Xz @ Yz.T) / n
Rp = np.clip(R, -0.999999999, 0.999999999)
tval = Rp * np.sqrt(dfree / (1.0 - Rp**2))
P = 2.0 * tdist.sf(np.abs(tval), df=dfree)

out_long = pd.DataFrame({
    "SNF_gene": np.repeat(snf_ids, len(new_ids)),
    "New_gene": np.tile(new_ids, len(snf_ids)),
    "r": R.reshape(-1),
    "p": P.reshape(-1),
}).dropna(subset=["r"])

odir = f"{outp}.snf_new_corr"
out_long.to_csv(f"{odir}/{outp}.SNF_vs_New.all_pairs.tsv", sep="\t", index=False)

top_list = []
for i, snf in enumerate(snf_ids):
    rr = R[i, :]
    pp = P[i, :]
    idx = np.argsort(-np.abs(rr))[:5]
    for j in idx:
        if np.isnan(rr[j]): 
            continue
        top_list.append([snf, new_ids[j], rr[j], pp[j], abs(rr[j])])

top = pd.DataFrame(top_list, columns=["SNF_gene","New_gene","r","p","abs_r"])
top = top.sort_values(["SNF_gene","abs_r"], ascending=[True, False])
top.to_csv(f"{odir}/{outp}.SNF_top5_newgenes.tsv", sep="\t", index=False)

uniq_new = top["New_gene"].unique().tolist()
mat = pd.DataFrame(index=snf_ids, columns=uniq_new, dtype=float)
for _, row in top.iterrows():
    mat.loc[row["SNF_gene"], row["New_gene"]] = row["r"]
mat.to_csv(f"{odir}/{outp}.SNF_top5_matrix.tsv", sep="\t", na_rep="NA")

with open(f"{odir}/{outp}.summary.txt","w") as f:
    f.write(f"Samples: {n}\nSNF genes used: {len(snf_ids)}\nNew genes used: {len(new_ids)}\nAll pairs: {out_long.shape[0]}\nTop5 rows: {top.shape[0]}\n")

print("[OK] written:")
print(f"  {odir}/{outp}.SNF_vs_New.all_pairs.tsv")
print(f"  {odir}/{outp}.SNF_top5_newgenes.tsv")
print(f"  {odir}/{outp}.SNF_top5_matrix.tsv")
print(f"  {odir}/{outp}.summary.txt")
PY

# =========================
# 2) R: plot (write temp R script)
# =========================
RSCRIPT="${ODIR}/plot_top5.R"
cat > "$RSCRIPT" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
outp <- args[1]
odir <- paste0(outp, ".snf_new_corr")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

top_fn <- file.path(odir, paste0(outp, ".SNF_top5_newgenes.tsv"))
mat_fn <- file.path(odir, paste0(outp, ".SNF_top5_matrix.tsv"))

if (!file.exists(top_fn)) stop(paste("Missing:", top_fn))
if (!file.exists(mat_fn)) stop(paste("Missing:", mat_fn))

top <- read_tsv(top_fn, show_col_types = FALSE)
mat <- read_tsv(mat_fn, show_col_types = FALSE)

# Heatmap data
mat_long <- mat %>%
  pivot_longer(-1, names_to="New_gene", values_to="r") %>%
  rename(SNF_gene = 1) %>%
  filter(!is.na(r))

p1 <- ggplot(mat_long, aes(x=New_gene, y=SNF_gene, fill=r)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_gradient2(low="#2b6cb0", mid="white", high="#c53030", midpoint=0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    panel.grid = element_blank()
  ) +
  labs(x="Top New genes (union)", y="SNF genes", fill="Pearson r")

h <- max(3, 0.35*length(unique(mat_long$SNF_gene)))
ggsave(file.path(odir, paste0(outp, ".SNF_top5_heatmap.pdf")), p1, width=10, height=h)
ggsave(file.path(odir, paste0(outp, ".SNF_top5_heatmap.png")), p1, width=10, height=h, dpi=300)

# Dotplot
top2 <- top %>%
  group_by(SNF_gene) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p2 <- ggplot(top2, aes(x=rank, y=SNF_gene, size=abs_r, color=r)) +
  geom_point() +
  scale_color_gradient2(low="#2b6cb0", mid="grey70", high="#c53030", midpoint=0) +
  scale_x_continuous(breaks=1:5, labels=paste0("Top",1:5)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x=NULL, y="SNF genes", size="|r|", color="r")

h2 <- max(3, 0.35*length(unique(top2$SNF_gene)))
ggsave(file.path(odir, paste0(outp, ".SNF_top5_dotplot.pdf")), p2, width=7, height=h2)
ggsave(file.path(odir, paste0(outp, ".SNF_top5_dotplot.png")), p2, width=7, height=h2, dpi=300)
RS

Rscript "$RSCRIPT" "$OUT"

echo "[DONE] outputs in: $ODIR/"
