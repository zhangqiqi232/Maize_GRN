setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析")
#筛选关键TFs
#Bundle sheath
library(dplyr)
regulation_count_Bun<- KBoost_Bun %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
Bun_average_targets <- mean(regulation_count_Bun$n_targets)
# 使用泊松分布计算每个转录因子的p值
Bun_ppois <- regulation_count_Bun %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = Bun_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
Bun_key_TFs <- Bun_ppois %>%
  filter(p_value < significance_level)
# 找出调控边最多的转录因子
Bun_max_tf <- Bun_key_TFs[which.max(Bun_key_TFs$n_targets), ]
# 找出调控边最少的转录因子
Bun_min_tf <- Bun_key_TFs[which.min(Bun_key_TFs$n_targets), ]
write.xlsx(Bun_key_TFs, "Bun_key_TFs.xlsx")

#Mesophyll
regulation_count_Mesophyll<- KBoost_Mesophyll %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
Mesophyll_average_targets <- mean(regulation_count_Mesophyll$n_targets)
# 使用泊松分布计算每个转录因子的p值
Mesophyll_ppois <- regulation_count_Mesophyll %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = Mesophyll_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
Mesophyll_key_TFs <- Mesophyll_ppois %>%
  filter(p_value < significance_level)
# 找出调控边最多的转录因子
Mesophyll_max_tf <- Mesophyll_key_TFs[which.max(Mesophyll_key_TFs$n_targets), ]
# 找出调控边最少的转录因子
Mesophyll_min_tf <- Mesophyll_key_TFs[which.min(Mesophyll_key_TFs$n_targets), ]
write.xlsx(Mesophyll_key_TFs, "Mesophyll_key_TFs.xlsx")

#Guard
regulation_count_Guard<- KBoost_Guard %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
Guard_average_targets <- mean(regulation_count_Guard$n_targets)
# 使用泊松分布计算每个转录因子的p值
Guard_ppois <- regulation_count_Guard %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = Guard_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
Guard_key_TFs <- Guard_ppois %>%
  filter(p_value < significance_level)
# 找出调控边最多的转录因子
Guard_max_tf <- Guard_key_TFs[which.max(Guard_key_TFs$n_targets), ]
# 找出调控边最少的转录因子
Guard_min_tf <- Guard_key_TFs[which.min(Guard_key_TFs$n_targets), ]
write.xlsx(Guard_key_TFs, "Guard_key_TFs.xlsx")

#Pavement
regulation_count_Pavement<- KBoost_Pavement  %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
Pavement_average_targets <- mean(regulation_count_Pavement$n_targets)
# 使用泊松分布计算每个转录因子的p值
Pavement_ppois <- regulation_count_Pavement %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = Pavement_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
Pavement_key_TFs <- Pavement_ppois %>%
  filter(p_value < significance_level)
# 找出调控边最多的转录因子
Pavement_max_tf <- Pavement_key_TFs[which.max(Pavement_key_TFs$n_targets), ]
# 找出调控边最少的转录因子
Pavement_min_tf <- Pavement_key_TFs[which.min(Pavement_key_TFs$n_targets), ]
write.xlsx(Pavement_key_TFs, "Pavement_key_TFs.xlsx")

#Subsidiary
regulation_count_Subsidiary<- KBoost_Subsidiary  %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
Subsidiary_average_targets <- mean(regulation_count_Subsidiary$n_targets)
# 使用泊松分布计算每个转录因子的p值
Subsidiary_ppois <- regulation_count_Subsidiary %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = Subsidiary_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
Subsidiary_key_TFs <- Subsidiary_ppois %>%
  filter(p_value < significance_level)
# 找出调控边最多的转录因子
Subsidiary_max_tf <- Subsidiary_key_TFs[which.max(Subsidiary_key_TFs$n_targets), ]
# 找出调控边最少的转录因子
Subsidiary_min_tf <- Subsidiary_key_TFs[which.min(Subsidiary_key_TFs$n_targets), ]
write.xlsx(Subsidiary_key_TFs, "Subsidiary_key_TFs.xlsx")

#使用miRSM包计算五个细胞类型的相似性
#获取五个细胞类型共同的关键转录因子
# 提取每个数据框中的TF列
Bun_key_TFs <- Bun_key_TFs[[1]]
Guard_key_TFs <- Guard_key_TFs[[1]]
Mesophyll_key_TFs <- Mesophyll_key_TFs[[1]]
Pavement_key_TFs <- Pavement_key_TFs[[1]]
Subsidiary_key_TFs <- Subsidiary_key_TFs[[1]]
# 放到一个列表里
key_TFs_list <- list(
  Bundle_sheath = Bun_key_TFs,
  Guard = Guard_key_TFs,
  Mesophyll = Mesophyll_key_TFs,
  Pavement = Pavement_key_TFs,
  Subsidiary = Subsidiary_key_TFs
)
library(miRSM)
n <- length(key_TFs_list)
similarity_matrix <- matrix(NA, n, n)

rownames(similarity_matrix) <- names(key_TFs_list)
colnames(similarity_matrix) <- names(key_TFs_list)

for (i in 1:n) {
  for (j in 1:n) {
    similarity_matrix[i, j] <- module_group_sim(
      list(key_TFs_list[[i]]),  # 包成 list
      list(key_TFs_list[[j]]),  # 包成 list
      sim.method = "Simpson"
    )
  }
}
similarity_matrix
library(ggcorrplot)
library(ggcorrplot)

sim1_plot <- ggcorrplot(
  similarity_matrix,
  type = "lower",
  lab = TRUE,
  lab_size = 4
) +
  scale_fill_gradientn(
    colours = c("#f7e2e6", "#e44345"),  # 从浅粉到深红
    limits = c(0.25, 1)
  ) +
  coord_fixed(ratio = 1) +
  ggtitle("Similarity of Cell-type Specific TF Modules") +
  theme_minimal() +
  theme(  # 指定字体为 Arial
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

#计算细胞类型关键TFs的度中心性数
# 定义一个函数计算每个关键TF调控的靶基因数量
calc_outdegree <- function(grn_df, key_tfs) {
  grn_df %>%
    group_by(TF) %>%
    summarise(out_degree = n()) %>%   # 计算每个TF的靶基因数
    filter(TF %in% key_tfs)
}
# 分别计算五个细胞类型的关键TF调控数
deg_Bun <- calc_outdegree(KBoost_Bun, Bun_key_TFs)
deg_Guard <- calc_outdegree(KBoost_Guard, Guard_key_TFs)
deg_Mesophyll <- calc_outdegree(KBoost_Mesophyll, Mesophyll_key_TFs)
deg_Pavement <- calc_outdegree(KBoost_Pavement, Pavement_key_TFs)
deg_Subsidiary <- calc_outdegree(KBoost_Subsidiary, Subsidiary_key_TFs)
# 合并到列表
deg_list <- list(
  Bundle_sheath = deg_Bun,
  Guard = deg_Guard,
  Mesophyll = deg_Mesophyll,
  Pavement = deg_Pavement,
  Subsidiary = deg_Subsidiary
)
# 获取所有关键TF的并集
all_tfs <- unique(unlist(lapply(deg_list, function(x) x$TF)))

# 转换为矩阵
deg_matrix <- sapply(deg_list, function(df) {
  v <- setNames(df$out_degree, df$TF)
  v[all_tfs]
})
rownames(deg_matrix) <- all_tfs
deg_matrix[is.na(deg_matrix)] <- 0  # 缺失填充为0
library(pheatmap)
library(pheatmap)
library(RColorBrewer)

# 自定义颜色梯度（无白色）
color_palette <- colorRampPalette(c("#4c7cb8", "#c2dfed", "#FFFFBF", "#F46D43", "#A50026"))(100)

# 绘制热图
pheatmap(
  deg_matrix,
  color = color_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  border_color = NA,
  show_rownames = FALSE,     # 不显示TF名称
  show_colnames = TRUE,      # 显示细胞类型
  main = "Degree Centrality of Key TFs across Five Cell Types",
  fontsize_col = 12
)
