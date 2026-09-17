#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(ggrepel)
  library(remotes)
  library(scTenifoldKnk)
  library(glue)
  library(gghalves)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 0) 参数
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法：Rscript 10.scTenifoldKnk_replicates.R gKO", call. = FALSE)

gKO <- args[1]
celltype <- args[2]

# 运行次数
N_RUN <- 100

# 子集数量
K_SET <- c(10, 30, 50, 70, 100)

# hit阈值（用于次数统计）
LOG2FC_CUT <- 1
PADJ_CUT <- 0.05

# 细胞选择（与你当前脚本一致）
N_CELLS_USE <- 2000

# HVG数
N_HVG <- 1000

# ============================================================
# 1) 输出目录
# ============================================================
dir.create("01_runs", showWarnings = FALSE)
dir.create("02_summaries", showWarnings = FALSE)
dir.create("03_plots", showWarnings = FALSE)

# ============================================================
# 2) 数据准备：提取count矩阵 + SNF重命名
# ============================================================
obj5 <- readRDS("zz_00_pea.rds")
SNF  <- read.table("All_species_merged_by_GifuT2T.SNF_Symbol.tsv", header = TRUE)

SNF_fix <- SNF %>%
  mutate(GifuT2T_fix = gsub("_", "-", GifuT2T))

Idents(obj5) <- "celltype"
message("celltype table:")
print(table(obj5$celltype))

# 提取该细胞类型的特异表达基因
deg_ct <- FindMarkers(
  object = obj5,
  ident.1 = celltype,
  ident.2 = NULL
)

deg_ct_sig <- deg_ct %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 1
  )

genes_ct_all <- rownames(deg_ct_sig)

# 早期信号响应细胞（与你当前一致：idents="7"）
sc.s <- subset(obj5, idents = celltype)

count <- GetAssayData(sc.s, slot = "counts")

# HVGs
hvg <- VariableFeatures(sc.s)
if (length(hvg) < N_HVG) {
  warning(glue("HVG数量不足({length(hvg)}), 将使用全部HVG。"))
  hvg_use <- hvg
} else {
  hvg_use <- hvg[1:N_HVG]
}

genes_use <- unique(c(hvg_use, gKO, genes_ct_all,SNF_fix$GifuT2T_fix))
cells_use <- colnames(count)[1:min(N_CELLS_USE, ncol(count))]

count_use <- count[genes_use, cells_use]

# --- 将SNF的基因号替换为GeneSymbol（处理同名重复加后缀） ---
genes_count <- rownames(count_use)

snf_map <- SNF_fix %>%
  filter(GifuT2T_fix %in% genes_count) %>%
  select(GeneSymbol, GifuT2T_fix) %>%
  group_by(GeneSymbol) %>%
  mutate(
    index = row_number(),
    NewGene = if (n() > 1) paste0(GeneSymbol, "-", index) else GeneSymbol
  ) %>%
  ungroup()

rename_table <- data.frame(OldGene = genes_count, stringsAsFactors = FALSE) %>%
  left_join(snf_map %>% select(GifuT2T_fix, NewGene),
            by = c("OldGene" = "GifuT2T_fix")) %>%
  mutate(FinalGene = ifelse(is.na(NewGene), OldGene, NewGene))

count_use_renamed <- count_use
rownames(count_use_renamed) <- rename_table$FinalGene

write.table(
  rename_table,
  file = "02_summaries/SNF_gene_rename_map.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 虚拟敲除的基因名称要和矩阵保持一致
if (gKO %in% rename_table$OldGene){
  gKO_raw <- gKO
  gKO <- rename_table[ rename_table$OldGene %in% gKO_raw, "FinalGene"]
}

# ============================================================
# 3) scTenifoldKnk重复运行：100次并保留结果
# ============================================================
run_one <- function(seed_i, count_mat, gKO) {
  set.seed(seed_i)

  res <- scTenifoldKnk(
    countMatrix = count_mat,
    qc = TRUE,
    gKO = gKO,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = min(500, ncol(count_mat)),
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-5,
    td_nDecimal = 3,
    ma_nDim = 2,
    nCores = parallel::detectCores()
  )

  df <- res[["diffRegulation"]] %>% as.data.frame()
  # 统一列名：确保 gene / FC / p.adj 存在
  if (!("gene" %in% colnames(df))) stop("diffRegulation中缺少 gene 列")
  if (!("FC" %in% colnames(df))) stop("diffRegulation中缺少 FC 列")
  if (!("p.adj" %in% colnames(df))) stop("diffRegulation中缺少 p.adj 列")

  df <- df %>%
    mutate(
      log2FC = log2(FC),
      logp.adj = -log10(pmax(p.adj, 1e-300))
    ) %>%
    select(gene, distance, Z, FC, p.value, p.adj, log2FC, logp.adj, everything())

  df
}

# 主循环（带失败保护）
all_runs <- vector("list", N_RUN)
run_status <- tibble(
  run = 1:N_RUN,
  seed = 1000 + 1:N_RUN,
  ok = FALSE,
  msg = NA_character_
)

message(glue("Start {N_RUN} runs ..."))
for (i in 1:N_RUN) {
  seed_i <- 1000 + i
  message(glue("[Run {i}/{N_RUN}] seed={seed_i}"))

  out_rds <- file.path("01_runs", glue("run_{sprintf('%03d', i)}.rds"))
  out_csv <- file.path("01_runs", glue("run_{sprintf('%03d', i)}.csv"))

  df_i <- tryCatch(
    {
      df <- run_one(seed_i, count_use_renamed, gKO)

      saveRDS(df, out_rds)
      write.csv(df, out_csv, row.names = FALSE)

      run_status$ok[i] <- TRUE
      run_status$msg[i] <- "OK"
      df
    },
    error = function(e) {
      run_status$ok[i] <- FALSE
      run_status$msg[i] <- as.character(e$message)
      warning(glue("[Run {i}] FAILED: {e$message}"))
      NULL
    }
  )

  all_runs[[i]] <- df_i
}

all_runs <- lapply(
  all_runs,
  function(df) {
    if (is.null(df)) return(NULL)
    df %>% dplyr::filter(gene != gKO)
  }
)

write.csv(run_status, "02_summaries/run_status.csv", row.names = FALSE)
saveRDS(all_runs, "02_summaries/all_runs_list.rds")

ok_idx <- which(run_status$ok)
if (length(ok_idx) < max(K_SET)) {
  warning(glue(
    "成功运行次数只有 {length(ok_idx)}，不足以覆盖最大K={max(K_SET)}。将只对可用范围内的K做分析。"
  ))
}
# 按顺序取成功的前若干次
all_runs_ok <- all_runs[ok_idx]
N_OK <- length(all_runs_ok)

K_SET2 <- K_SET[K_SET <= N_OK]
if (length(K_SET2) == 0) stop("没有足够的成功运行结果用于后续分析。")

# ============================================================
# 4) (2) 相关性：各K下，K次结果两两FC Spearman相关 + 箱线图 + 均值拟合
# ============================================================
# 将每次结果转为 named vector: gene -> FC
to_fc_vec <- function(df) {
  v <- df$FC
  names(v) <- df$gene
  v
}

# 计算两次run的spearman相关（按gene对齐）
spearman_fc <- function(v1, v2) {
  g <- intersect(names(v1), names(v2))
  if (length(g) < 10) return(NA_real_)
  suppressWarnings(cor(v1[g], v2[g], method = "spearman", use = "pairwise.complete.obs"))
}

cor_records <- list()
mean_records <- list()

for (K in K_SET2) {
  runsK <- all_runs_ok[1:K]
  vecs <- lapply(runsK, to_fc_vec)

  cors <- c()
  if (K >= 2) {
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        cors <- c(cors, spearman_fc(vecs[[i]], vecs[[j]]))
      }
    }
  }

  cor_df <- tibble(
    K = K,
    pair_id = seq_along(cors),
    rho = cors
  ) %>% filter(!is.na(rho))

  cor_records[[as.character(K)]] <- cor_df

  mean_records[[as.character(K)]] <- tibble(
    K = K,
    mean_rho = mean(cor_df$rho, na.rm = TRUE)
  )
}

cor_all <- bind_rows(cor_records)
mean_all <- bind_rows(mean_records)

write.csv(cor_all, "02_summaries/spearman_pairwise_all.csv", row.names = FALSE)
write.csv(mean_all, "02_summaries/spearman_mean_byK.csv", row.names = FALSE)

cor_all <- cor_all %>%
  mutate(K = factor(K, levels = c(10, 30, 50, 70, 100)))

mean_all <- mean_all %>%
  mutate(K = factor(K, levels = c(10, 30, 50, 70, 100)))

mean_all <- mean_all %>%
  mutate(K_index = as.numeric(K))


## -----------------------------
## 1. 数据准备（再次确认）
## -----------------------------
cor_all <- cor_all %>%
  mutate(K = factor(K, levels = c("10", "30", "50", "70", "100")))

mean_all <- mean_all %>%
  mutate(
    K = factor(K, levels = c("10", "30", "50", "70", "100")),
    K_index = as.numeric(K)   # 仅用于趋势线定位
  )

## -----------------------------
## 2. 作图（CNS / PNAS 风格）
## -----------------------------
p_cor <- ggplot() +

  ## 箱线图（主体，黑色）
  geom_boxplot(
    data = cor_all,
    aes(x = K, y = rho),
    width = 0.55,
    outlier.size = 0.6,
    outlier.alpha = 0.25,
    linewidth = 0.4,
    color = "black",
    fill = NA
  ) +

  ## 均值点（小而克制）
  stat_summary(
    data = cor_all,
    aes(x = K, y = rho),
    fun = mean,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +

  ## 均值趋势线（强调色，非回归）
  geom_line(
    data = mean_all,
    aes(x = K_index, y = mean_rho, group = 1),
    linewidth = 0.9,
    color = "#B2182B"   # CNS 常用深红
  ) +

  ## 均值锚点
  geom_point(
    data = mean_all,
    aes(x = K_index, y = mean_rho),
    size = 2.2,
    color = "#B2182B"
  ) +

  ## X 轴：离散标签
  scale_x_discrete(
    limits = c("10", "30", "50", "70", "100"),
    name = "Number of scTenifoldKnk runs (K)"
  ) +

  ## Y 轴
  scale_y_continuous(
    name = "Spearman correlation of FC (pairwise)",
    expand = expansion(mult = c(0.05, 0.1))
  ) +

  ## 主题（PNAS / CNS 通用）
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(
      hjust = 0.5,
      face = "plain",
      size = 11
    )
  ) +

  labs(
    title = glue("Reproducibility of inferred FC across repeated scTenifoldKnk runs (gKO = {gKO})")
  )

ggsave("03_plots/01_FC_spearman_boxplot_with_mean_fit.pdf", p_cor, width = 8, height = 5)

# ============================================================
# 5) (3) 各K下：每个基因 log2FC>1 且 p.adj<0.05 的次数
# ============================================================
# 以及 p.adj<0.05 的次数（给第4步上色用）
count_hits_for_K <- function(runs_list, K, log2fc_cut, padj_cut) {
  runsK <- runs_list[1:K]

  # gene全集（通常一致，但稳妥起见取union）
  gene_all <- sort(unique(unlist(lapply(runsK, function(df) df$gene))))

  hit_counts <- setNames(rep(0L, length(gene_all)), gene_all)
  sigp_counts <- setNames(rep(0L, length(gene_all)), gene_all)

  for (df in runsK) {
    df2 <- df %>% select(gene, log2FC, p.adj)

    hit_gene <- df2 %>% filter(log2FC < log2fc_cut, p.adj < padj_cut) %>% pull(gene)
    sig_gene <- df2 %>% filter(p.adj < padj_cut) %>% pull(gene)

    hit_counts[hit_gene] <- hit_counts[hit_gene] + 1L
    sigp_counts[sig_gene] <- sigp_counts[sig_gene] + 1L
  }

  tibble(
    gene = gene_all,
    K = K,
    hit_n = as.integer(hit_counts[gene_all]),
    sigp_n = as.integer(sigp_counts[gene_all])
  )
}

hit_tables <- lapply(K_SET2, function(K) {
  count_hits_for_K(all_runs_ok, K, LOG2FC_CUT, PADJ_CUT)
})
hit_all <- bind_rows(hit_tables)

write.csv(hit_all, "02_summaries/gene_hit_counts_byK.csv", row.names = FALSE)

# 同时输出宽表（便于查看）
hit_wide <- hit_all %>%
  select(gene, K, hit_n) %>%
  pivot_wider(names_from = K, values_from = hit_n, names_prefix = "K_")

sigp_wide <- hit_all %>%
  select(gene, K, sigp_n) %>%
  pivot_wider(names_from = K, values_from = sigp_n, names_prefix = "K_")

write.csv(hit_wide, "02_summaries/gene_hit_counts_byK_wide.csv", row.names = FALSE)
write.csv(sigp_wide, "02_summaries/gene_sigp_counts_byK_wide.csv", row.names = FALSE)

# ============================================================
# 6) (4) 各K下：Top20基因云雨图（雨=散点，云=小提琴）
#     - 横轴：FC
#     - 纵轴：gene
#     - 基因排序：按FC均值排序
#     - 颜色：由“p.adj<0.05的次数(sigp_n)”控制（渐变）
# ============================================================
get_fc_long_for_K <- function(runs_list, K) {
  runsK <- runs_list[1:K]
  fc_long <- bind_rows(lapply(seq_along(runsK), function(i) {
    df <- runsK[[i]]
    df %>%
      transmute(
        run = i,
        gene = gene,
        FC = FC,
        log2FC = log2FC,
        p.adj = p.adj
      )
  }))
  fc_long
}

make_raincloud_plot <- function(K, fc_long, hit_tbl_K, out_pdf) {
  # Top20：按 hit_n 降序
  top_genes <- hit_tbl_K %>%
    arrange(desc(hit_n), desc(sigp_n)) %>%
    slice_head(n = 20) %>%
    pull(gene)

  dat <- fc_long %>% filter(gene %in% top_genes)

  # gene排序：按FC均值（对这K次）
  gene_order <- dat %>%
    group_by(gene) %>%
    summarise(mean_FC = mean(FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_FC) %>%
    pull(gene)

  # 合并颜色信息：sigp_n（p.adj<0.05次数）
  sigp_map <- hit_tbl_K %>% select(gene, sigp_n, hit_n)

  dat2 <- dat %>%
    left_join(sigp_map, by = "gene") %>%
    mutate(gene = factor(gene, levels = gene_order))

  # 云雨图：violin + boxplot + jitter
  p <- ggplot(dat2, aes(x = gene, y = log2FC)) +

    ## 上半部小提琴（云）
    geom_half_violin(
      aes(fill = sigp_n),
      side = "r",
      trim = TRUE,
      width = 0.9,
      alpha = 0.85,
      color = NA
    ) +

    ## 箱线（雨的“核心”）
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.35
    ) +

    ## 散点（雨）
    geom_jitter(
      aes(color = sigp_n),
      width = 0.08,
      height = 0,
      size = 1.2,
      alpha = 0.6
    ) +

    coord_flip() +

    ## CNS 风格连续渐变（小提琴 & 点共用）
    scale_fill_gradient(
      low = "grey85",
      high = "#B2182B"   # Nature / CNS 常用暗红
    ) +
    scale_color_gradient(
      low = "grey60",
      high = "#B2182B"
    ) +

    labs(
      title = glue(
        "Top 20 genes ranked by hit frequency (log2FC > {LOG2FC_CUT}, adj. P < {PADJ_CUT})\n",
        "K = {K}, gKO = {gKO}"
      ),
      x = "Gene (ordered by mean log2FC)",
      y = "log2(FC)",
      fill  = "Count (adj. P < 0.05)",
      color = "Count (adj. P < 0.05)"
    ) +

    theme_classic(base_size = 11) +
    theme(
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(
        hjust = 0,
        size = 11,
        face = "plain"
      )
    )

  ggsave(out_pdf, p, width = 8.5, height = 6)

  p
}

for (K in K_SET2) {
  message(glue("Make raincloud plot for K={K}"))

  fc_long_K <- get_fc_long_for_K(all_runs_ok, K)

  hit_tbl_K <- hit_all %>% filter(K == !!K)

  out_pdf <- file.path("03_plots", glue("02_raincloud_Top20_byHit_K{K}.pdf"))
  make_raincloud_plot(K, fc_long_K, hit_tbl_K, out_pdf)
}

# ============================================================
# 7) 额外输出：每个K下 Top20基因列表（方便检查）
# ============================================================
top20_list <- bind_rows(lapply(K_SET2, function(K) {
  (hit_all %>% filter(K == !!K) %>%
     arrange(desc(hit_n), desc(sigp_n)) %>%
     slice_head(n = 20) %>%
     mutate(rank = row_number()) %>%
     select(K, rank, gene, hit_n, sigp_n))
}))

write.csv(top20_list, "02_summaries/top20_genes_byHit_eachK.csv", row.names = FALSE)

message("All done.")
