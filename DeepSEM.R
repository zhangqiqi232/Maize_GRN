setwd("D:/玉米多模态数据的基因调控网络研究/DeepSEM")
#生成先验数据
#Bundle sheath
# 假设 TF_Bun 是一个字符串向量
# 假设 expr_matrix 是一个表达矩阵，行是基因，列是样本
# 获取基因名（靶基因）
target_genes <- rownames(expr_Bun)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
pirort_Bun <- expand.grid(
  TF = tf_Bun,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(pirort_Bun) <- c("Gene1", "Gene2")
write.csv(pirort_Bun, file = "pirort_Bun.csv", row.names = FALSE)

#Mesophyll
expr_Mesophyll_t <- t(expr_Mesophyll)
write.csv(expr_Mesophyll_t, file = "expr_Mesophyll_t.csv")
target_genes <- rownames(expr_Mesophyll)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
pirort_Mesophyll <- expand.grid(
  TF = tf_Mesophyll,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(pirort_Mesophyll) <- c("Gene1", "Gene2")
write.csv(pirort_Mesophyll, file = "pirort_Mesophyll.csv", row.names = FALSE)

#Guard
expr_Guard_t <- t(expr_Guard)
write.csv(expr_Guard_t, file = "expr_Guard_t.csv")
target_genes <- rownames(expr_Guard)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
prior_Guard <- expand.grid(
  TF = tf_Guard,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(prior_Guard) <- c("Gene1", "Gene2")
write.csv(prior_Guard, file = "prior_Guard.csv", row.names = FALSE)

#Pavement
expr_Pavement_t <- t(expr_Pavement)
write.csv(expr_Pavement_t, file = "expr_Pavement_t.csv")
target_genes <- rownames(expr_Pavement)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
prior_Pavement <- expand.grid(
  TF = tf_Pavement,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(prior_Pavement) <- c("Gene1", "Gene2")
write.csv(prior_Pavement, file = "prior_Pavement.csv", row.names = FALSE)

#Subsidiary
expr_Subsidiary_t <- t(expr_Subsidiary)
write.csv(expr_Subsidiary_t, file = "expr_Subsidiary_t.csv")
target_genes <- rownames(expr_Subsidiary)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
prior_Subsidiary <- expand.grid(
  TF = tf_Subsidiary,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(prior_Subsidiary) <- c("Gene1", "Gene2")
write.csv(prior_Subsidiary, file = "prior_Subsidiary.csv", row.names = FALSE)


#bulk
expr_leaf_t <- t(leaf_matrix)
write.csv(expr_leaf_t, file = "expr_leaf_t.csv")
target_genes <- rownames(leaf_matrix)
# 创建转录因子和靶基因的所有组合（笛卡尔积）
prior_leaf <- expand.grid(
  TF = tf_leaf,
  Target = target_genes,
  stringsAsFactors = FALSE
)
colnames(prior_leaf) <- c("Gene1", "Gene2")
write.csv(prior_leaf, file = "prior_leaf.csv", row.names = FALSE)

#读取DeepSEM推断的GRN筛选Z分值大于0.5的关系
#Bundle sheath
library(vroom)
DeepSEM_Bun<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/GRN_Bundle.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
DeepSEM_Bun$EdgeWeight.zscore <- transform_zscore(DeepSEM_Bun$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_Bun <- biadjacency_matrix(DeepSEM_Bun$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_Bun <- ifelse(DeepSEM_biadjacency_Bun == TRUE, 1, 0)
DeepSEM_Bun <- DeepSEM_Bun[apply(DeepSEM_biadjacency_Bun, 1, function(x) any(x == 1)), ]
save(DeepSEM_Bun,file='DeepSEM_Bun.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepSEM_Bun <- DeepSEM_Bun[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(397重叠调控对)
DeepSEM_overlap_networks_Bun <- merge(DeepSEM_Bun, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(99个TFs)
DeepSEM_unique_TF_Bun <- length(unique(DeepSEM_Bun$TF))
#统计不重复的 Target数量(18443个Targets)
DeepSEM_unique_Target_Bun <- length(unique(DeepSEM_Bun$Target))
#使用超几何分布来验证p值
N <- 1825857 # 背景基因总数(TFxTarget)
K <- 116712   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 397 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 1.914881e-36 


#Mesophyll
#读取tsv文件
library(vroom)
DeepSEM_Mesophyll<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/Mesophyll.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
DeepSEM_Mesophyll$EdgeWeight.zscore <- transform_zscore(DeepSEM_Mesophyll$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_Mesophyll <- biadjacency_matrix(DeepSEM_Mesophyll$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_Mesophyll <- ifelse(DeepSEM_biadjacency_Mesophyll == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepSEM_Mesophyll <- DeepSEM_Mesophyll[apply(DeepSEM_biadjacency_Mesophyll, 1, function(x) any(x == 1)), ]
save(DeepSEM_Mesophyll,file='DeepSEM_Mesophyll.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepSEM_Mesophyll <- DeepSEM_Mesophyll[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(217重叠调控对)
DeepSEM_overlap_networks_Mesophyll <- merge(DeepSEM_Mesophyll, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(121个TFs)
DeepSEM_unique_TF_Mesophyll <- length(unique(DeepSEM_Mesophyll$TF))
#统计不重复的 Target数量(24951个Targets)
DeepSEM_unique_Target_Mesophyll <- length(unique(DeepSEM_Mesophyll$Target))
#使用超几何分布来验证p值
N <- 3019071 # 背景基因总数(TFxTarget)
K <- 184717   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 217 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.05101136 

#Guard
#读取tsv文件
library(vroom)
DeepSEM_Guard<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/Guard.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
DeepSEM_Guard$EdgeWeight.zscore <- transform_zscore(DeepSEM_Guard$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_Guard <- biadjacency_matrix(DeepSEM_Guard$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_Guard <- ifelse(DeepSEM_biadjacency_Guard == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepSEM_Guard <- DeepSEM_Guard[apply(DeepSEM_biadjacency_Guard, 1, function(x) any(x == 1)), ]
save(DeepSEM_Guard,file='DeepSEM_Guard.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepSEM_Guard <- DeepSEM_Guard[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(319重叠调控对)
DeepSEM_overlap_networks_Guard <- merge(DeepSEM_Guard, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(107个TFs)
DeepSEM_unique_TF_Guard <- length(unique(DeepSEM_Guard$TF))
#统计不重复的 Target数量(21607个Targets)
DeepSEM_unique_Target_Guard <- length(unique(DeepSEM_Guard$Target))
#使用超几何分布来验证p值
N <- 2311949 # 背景基因总数(TFxTarget)
K <- 160607   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 319 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 5.265017e-11 

#Pavement
#读取tsv文件
library(vroom)
DeepSEM_Pavement<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/Pavement.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
DeepSEM_Pavement$EdgeWeight.zscore <- transform_zscore(DeepSEM_Pavement$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_Pavement <- biadjacency_matrix(DeepSEM_Pavement$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_Pavement <- ifelse(DeepSEM_biadjacency_Pavement == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepSEM_Pavement <- DeepSEM_Pavement[apply(DeepSEM_biadjacency_Pavement, 1, function(x) any(x == 1)), ]
save(DeepSEM_Pavement,file='DeepSEM_Pavement.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepSEM_Pavement <- DeepSEM_Pavement[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(241重叠调控对)
DeepSEM_overlap_networks_Pavement <- merge(DeepSEM_Pavement, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(112个TFs)
DeepSEM_unique_TF_Pavement <- length(unique(DeepSEM_Pavement$TF))
#统计不重复的 Target数量(22973个Targets)
DeepSEM_unique_Target_Pavement <- length(unique(DeepSEM_Pavement$Target))
#使用超几何分布来验证p值
N <- 2572976 # 背景基因总数(TFxTarget)
K <- 154636   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 241 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.0001645635

#Subsidiary
#读取tsv文件
library(vroom)
DeepSEM_Subsidiary<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/Subsidiary.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
DeepSEM_Subsidiary$EdgeWeight.zscore <- transform_zscore(DeepSEM_Subsidiary$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_Subsidiary <- biadjacency_matrix(DeepSEM_Subsidiary$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_Subsidiary <- ifelse(DeepSEM_biadjacency_Subsidiary == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepSEM_Subsidiary <- DeepSEM_Subsidiary[apply(DeepSEM_biadjacency_Subsidiary, 1, function(x) any(x == 1)), ]
save(DeepSEM_Subsidiary,file='DeepSEM_Subsidiary.Rdata')
# 保留前两列
DeepSEM_Subsidiary <- DeepSEM_Subsidiary[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(266重叠调控对)
DeepSEM_overlap_networks_Subsidiary <- merge(DeepSEM_Subsidiary, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(112个TFs)
DeepSEM_unique_TF_Subsidiary <- length(unique(DeepSEM_Subsidiary$TF))
#统计不重复的 Target数量(22292个Targets)
DeepSEM_unique_Target_Subsidiary <- length(unique(DeepSEM_Subsidiary$Target))
#使用超几何分布来验证p值
N <- 2496707 # 背景基因总数(TFxTarget)
K <- 141870   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 266 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 4.193465e-10

#bulk数据
#读取tsv文件
library(vroom)
Deep_SEM_leaf<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/leaf.tsv", show_col_types = FALSE)
#利用z分值筛选基因调控网络
transform_zscore <- function(vec){
  vec.mean <- mean(vec, na.rm = TRUE)
  vec.sd <- sd(vec, na.rm = TRUE)
  vec.zscore <- (vec - vec.mean) / vec.sd
  return(vec.zscore)
}
biadjacency_matrix <- function(zscore_matrix, pvalue.cutoff = 0.05){
  zscore_cutoff <- -qnorm(pvalue.cutoff)
  biadjacency <- as.matrix(zscore_matrix > zscore_cutoff)
  return(biadjacency)
}
length(unique(Deep_SEM_leaf$Target))
Deep_SEM_leaf$EdgeWeight.zscore <- transform_zscore(Deep_SEM_leaf$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
biadjacency_leaf <- biadjacency_matrix(Deep_SEM_leaf$EdgeWeight.zscore, pvalue.cutoff = 0.05)
biadjacency_leaf <- ifelse(biadjacency_leaf == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
Deep_SEM_leaf <- Deep_SEM_leaf[apply(biadjacency_leaf, 1, function(x) any(x == 1)), ]
