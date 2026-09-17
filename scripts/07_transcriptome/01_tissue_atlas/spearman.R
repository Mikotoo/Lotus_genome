#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

options(stringsAsFactors = FALSE)

# ===============================
# 1. 输入文件
# ===============================
gifu_file <- "Gifu_all_samples_mean.tsv"
mg20_file <- "MG20_all_samples_mean.tsv"
out_prefix <- "Gifu_MG20"

# ===============================
# 2. 读入数据
# ===============================
gifu <- read_tsv(gifu_file, show_col_types = FALSE)
mg20 <- read_tsv(mg20_file, show_col_types = FALSE)

# ===============================
# 3. 转为 long format
#    并去掉 Gifu_root_dpi02
# ===============================
gifu_long <- gifu %>%
  pivot_longer(
    -gene_id,
    names_to = "sample",
    values_to = "expr_gifu"
  ) %>%
  filter(sample != "Gifu_root_dpi02") %>%
  mutate(
    sample = str_remove(sample, "^Gifu_")
  )

mg20_long <- mg20 %>%
  pivot_longer(
    -gene_id,
    names_to = "sample",
    values_to = "expr_mg20"
  ) %>%
  mutate(
    sample = str_remove(sample, "^MG20_")
  )

# ===============================
# 4. 合并 Gifu & MG20
# ===============================
expr_merge_raw <- gifu_long %>%
  inner_join(mg20_long, by = c("gene_id", "sample"))

# ===============================
# 5. 不添加器官分类：
#    按每个 sample 计算表达相似性
# ===============================
sample_corr <- expr_merge_raw %>%
  group_by(sample) %>%
  filter(
    sum(!is.na(expr_gifu) & !is.na(expr_mg20)) > 1,
    sd(expr_gifu, na.rm = TRUE) > 0,
    sd(expr_mg20, na.rm = TRUE) > 0
  ) %>%
  summarise(
    spearman_r = cor(
      expr_gifu,
      expr_mg20,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    n_gene = n_distinct(gene_id),
    .groups = "drop"
  )

write_tsv(
  sample_corr,
  paste0(out_prefix, ".sample_expression_spearman.tsv")
)

# sample 排序
sample_corr <- sample_corr %>%
  mutate(
    sample = factor(sample, levels = sample)
  )

p_sample_corr <- ggplot(sample_corr, aes(x = sample, y = spearman_r)) +
  geom_col(fill = "#59A14F", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", spearman_r)),
    vjust = -0.5,
    size = 3.2
  ) +
  ylim(0, 1) +
  labs(
    x = NULL,
    y = "Spearman correlation (Gifu vs MG20)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggsave(
  paste0(out_prefix, ".sample_expression_correlation.pdf"),
  p_sample_corr,
  width = 9,
  height = 4.5
)

# ===============================
# 6. 器官分类函数
# ===============================
make_sample_info <- function(samples) {
  tibble(
    sample = samples,
    organ = case_when(
      str_detect(sample, "nodule_dpi10|nodule_dpi21|nodule_dpi45") ~ "nodule",
      str_detect(sample, "root_dpi10|root_dpi21|root_dpi45|root_hpi48") ~ "root_infect",
      str_detect(sample, "root_uni") ~ "root",
      str_detect(sample, "leaf") ~ "leaf",
      str_detect(sample, "stem") ~ "stem",
      str_detect(sample, "flower") ~ "flower",
      str_detect(sample, "pod") ~ "pod",
      TRUE ~ NA_character_
    )
  )
}

sample_info <- make_sample_info(
  intersect(gifu_long$sample, mg20_long$sample)
)

# ===============================
# 7. 添加器官分类后：
#    按 organ 计算表达相似性
# ===============================
expr_merge <- expr_merge_raw %>%
  left_join(sample_info, by = "sample") %>%
  filter(!is.na(organ))

organ_corr <- expr_merge %>%
  group_by(organ) %>%
  filter(
    sum(!is.na(expr_gifu) & !is.na(expr_mg20)) > 1,
    sd(expr_gifu, na.rm = TRUE) > 0,
    sd(expr_mg20, na.rm = TRUE) > 0
  ) %>%
  summarise(
    spearman_r = cor(
      expr_gifu,
      expr_mg20,
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    n_gene = n_distinct(gene_id),
    .groups = "drop"
  )

write_tsv(
  organ_corr,
  paste0(out_prefix, ".organ_expression_spearman.tsv")
)

p_corr <- ggplot(organ_corr, aes(x = organ, y = spearman_r)) +
  geom_col(fill = "#4E79A7", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", spearman_r)),
    vjust = -0.5,
    size = 3.5
  ) +
  ylim(0, 1) +
  labs(
    x = NULL,
    y = "Spearman correlation (Gifu vs MG20)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  paste0(out_prefix, ".organ_expression_correlation.pdf"),
  p_corr,
  width = 7,
  height = 4
)