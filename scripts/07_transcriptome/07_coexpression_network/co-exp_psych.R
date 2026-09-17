library(tidyverse)
library(data.table)
library(psych)
library(igraph)

# DOI: 10.1126/sciadv.adt1113

# ===============================
# 1. 参数 & 数据读入
# ===============================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法：Rscript co_exp_SNF_all.R go_id_list.txt output_dir", call. = FALSE)
}

id_file <- args[1]
out_dir <- args[2]

genes <- fread(id_file, header = FALSE)$V1
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 表达矩阵
Gifu <- read.table("../../../Gifu_all_samples.tpm.tsv",
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
MG20 <- read.table("../../../MG20_all_samples.tpm.tsv",
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

datExpr0 <- cbind(Gifu, MG20)

# 仅保留根 / 接菌根 / 根瘤
noduleTraitList <- c("root_hpi48","root_uni","nodule")
keep_cols <- grep(paste(noduleTraitList, collapse = "|"),
                  colnames(datExpr0), value = TRUE)

datExpr1 <- datExpr0[rownames(datExpr0) %in% genes, keep_cols]

# ===============================
# 2. 基因分类信息
# ===============================

new <- fread("all_new_genes.list", col.names = "GeneID")
SNF <- fread("SNF_gene.list", col.names = "GeneID")

Ortho <- fread("All_species_merged_by_GifuT2T.tsv")

SNF1 <- Ortho %>%
  filter(MG20Old %in% SNF$GeneID) %>%
  pull(GifuT2T) %>%
  na.omit() %>%
  unique()

genes_all <- rownames(datExpr1)
SNF1_genes <- intersect(SNF1, genes_all)
new_genes  <- intersect(new$GeneID, genes_all)

# ===============================
# 3. 表达过滤（保持你原逻辑）
# ===============================

stages <- c("root_uni","root_hpi48",
            "nodule_dpi10","nodule_dpi21","nodule_dpi45")

pass_expression_filter <- function(gene, expr_mat, stages, cutoff = 0.5) {
  for (kw in stages) {
    cols <- grep(kw, colnames(expr_mat), value = TRUE)
    if (length(cols) == 0) next
    if (min(expr_mat[gene, cols]) > cutoff) return(TRUE)
  }
  FALSE
}

SNF1_genes <- SNF1_genes[
  sapply(SNF1_genes, pass_expression_filter,
         expr_mat = datExpr1, stages = stages)
]

genes_all <- genes_all[
  sapply(genes_all, pass_expression_filter,
         expr_mat = datExpr1, stages = stages)
]

# ===============================
# 4. 相关性计算（SNF × 所有基因）
# ===============================

expr_SNF <- t(datExpr1[SNF1_genes, ])
expr_all <- t(datExpr1[genes_all, ])

cor_res <- corr.test(expr_SNF, expr_all,
                     method = "pearson", adjust = "none")

edges_all <- expand.grid(
  gene1 = SNF1_genes,
  gene2 = genes_all,
  stringsAsFactors = FALSE
) %>%
  mutate(
    cor = as.vector(cor_res$r),
    p   = as.vector(cor_res$p)
  ) %>%
  filter(gene1 != gene2)

# ===============================
# 5. Bonferroni + cor > 0.7
# ===============================

alpha <- 0.05
p_cutoff <- alpha / nrow(edges_all)

edges_sig <- edges_all %>%
  filter(p < p_cutoff, cor > 0.7)

# ===============================
# 6. 输出 ①：SNF × 所有基因 网络
# ===============================

dir_all <- file.path(out_dir, "SNF_vs_ALL")
dir.create(dir_all, showWarnings = FALSE)

net_all <- graph_from_data_frame(edges_sig, directed = FALSE)

write.table(edges_sig,
            file.path(dir_all, "edges_bonferroni.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(GeneID = V(net_all)$name),
            file.path(dir_all, "nodes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(net_all,
        file.path(dir_all, "network_igraph.rds"))

# ===============================
# 7. 输出 ②：SNF–SNF & SNF–新基因 子网络
# ===============================

edges_sub <- edges_sig %>%
  filter(
    gene1 %in% SNF1_genes &
    (gene2 %in% SNF1_genes | gene2 %in% new_genes)
  )

dir_sub <- file.path(out_dir, "SNF_vs_New")
dir.create(dir_sub, showWarnings = FALSE)

net_sub <- graph_from_data_frame(edges_sub, directed = FALSE)

write.table(edges_sub,
            file.path(dir_sub, "edges_bonferroni.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(GeneID = V(net_sub)$name),
            file.path(dir_sub, "nodes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(net_sub,
        file.path(dir_sub, "network_igraph.rds"))

# ===============================
# 8. 网络组成统计
# ===============================

nodes_all <- V(net_all)$name

n_SNF <- sum(nodes_all %in% SNF1_genes)
n_new <- sum(nodes_all %in% new_genes)
n_other <- length(nodes_all) - n_SNF - n_new

stat <- data.frame(
  Category = c("SNF", "New_gene", "Other_coexpressed"),
  Count    = c(n_SNF, n_new, n_other)
)

write.table(stat,
            file.path(out_dir, "network_gene_composition.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Done.")
