setwd("D:/玉米多模态数据的基因调控网络研究/DeepRig")
#读取DeepRig推断出的GRN筛选Z分值大于0.5的关系
#Bundle sheath
DeepRig_Bun <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Bun.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_Bun) <- c("TF", "Target","EdgeWeight")
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
DeepRig_Bun$EdgeWeight.zscore <- transform_zscore(DeepRig_Bun$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_Bun <- biadjacency_matrix(DeepRig_Bun$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_Bun <- ifelse(DeepRig_biadjacency_Bun == TRUE, 1, 0)
DeepRig_Bun <- DeepRig_Bun[apply(DeepRig_biadjacency_Bun, 1, function(x) any(x == 1)), ]
save(DeepRig_Bun,file='DeepRig_Bun.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_Bun <- DeepRig_Bun[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(288重叠调控对)
DeepRig_overlap_networks_Bun <- merge(DeepRig_Bun, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(99个TFs)
DeepRig_unique_TF_Bun <- length(unique(DeepRig_Bun$TF))   
#统计不重复的 Target数量(6042个Targets)
DeepRig_unique_Target_Bun <- length(unique(DeepRig_Bun$Target))
#使用超几何分布来验证p值
N <- 598158 # 背景基因总数(TFxTarget)
K <- 185958   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 288 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 1

#Mesoohyll
DeepRig_Mesophyll <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Mesophyll.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_Mesophyll) <- c("TF", "Target","EdgeWeight")
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
DeepRig_Mesophyll$EdgeWeight.zscore <- transform_zscore(DeepRig_Mesophyll$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_Mesophyll <- biadjacency_matrix(DeepRig_Mesophyll$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_Mesophyll <- ifelse(DeepRig_biadjacency_Mesophyll == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepRig_Mesophyll <- DeepRig_Mesophyll[apply(DeepRig_biadjacency_Mesophyll, 1, function(x) any(x == 1)), ]
save(DeepRig_Mesophyll,file='DeepRig_Mesophyll.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_Mesophyll <- DeepRig_Mesophyll[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(1068重叠调控对)
DeepRig_overlap_networks_Mesophyll <- merge(DeepRig_Mesophyll, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(92个TFs)
DeepRig_unique_TF_Mesophyll <- length(unique(DeepRig_Mesophyll$TF))
#统计不重复的 Target数量(17846个Targets)
DeepRig_unique_Target_Mesophyll <- length(unique(DeepRig_Mesophyll$Target))
#使用超几何分布来验证p值
N <- 1641832 # 背景基因总数(TFxTarget)
K <- 685609   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 1068 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 1

#Pavement
DeepRig_Pavement <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Pavement.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_Pavement) <- c("TF", "Target","EdgeWeight")
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
DeepRig_Pavement$EdgeWeight.zscore <- transform_zscore(DeepRig_Pavement$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_Pavement <- biadjacency_matrix(DeepRig_Pavement$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_Pavement <- ifelse(DeepRig_biadjacency_Pavement == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepRig_Pavement <- DeepRig_Pavement[apply(DeepRig_biadjacency_Pavement, 1, function(x) any(x == 1)), ]
save(DeepRig_Pavement,file='DeepRig_Pavement.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_Pavement <- DeepRig_Pavement[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(没有重叠调控对)
DeepRig_overlap_networks_Pavement <- merge(DeepRig_Pavement, chip_seq, by = c("TF", "Target"))

#Guard
DeepRig_Guard <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Guard.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_Guard) <- c("TF", "Target","EdgeWeight")
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
DeepRig_Guard$EdgeWeight.zscore <- transform_zscore(DeepRig_Guard$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_Guard <- biadjacency_matrix(DeepRig_Guard$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_Guard <- ifelse(DeepRig_biadjacency_Guard == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepRig_Guard <- DeepRig_Guard[apply(DeepRig_biadjacency_Guard, 1, function(x) any(x == 1)), ]
save(DeepRig_Guard,file='DeepRig_Guard.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_Guard <- DeepRig_Guard[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(587重叠调控对)
DeepRig_overlap_networks_Guard <- merge(DeepRig_Guard, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(107个TFs)
DeepRig_unique_TF_Guard <- length(unique(DeepRig_Guard$TF))
#统计不重复的 Target数量(19525个Targets)
DeepRig_unique_Target_Guard <- length(unique(DeepRig_Guard$Target))
#使用超几何分布来验证p值
N <- 2089175 # 背景基因总数(TFxTarget)
K <- 426514   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 587 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.996914

#Subsidiary
DeepRig_Subsidiary <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Subsidiary.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_Subsidiary) <- c("TF", "Target","EdgeWeight")
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
DeepRig_Subsidiary$EdgeWeight.zscore <- transform_zscore(DeepRig_Subsidiary$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_Subsidiary <- biadjacency_matrix(DeepRig_Subsidiary$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_Subsidiary <- ifelse(DeepRig_biadjacency_Subsidiary == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepRig_Subsidiary <- DeepRig_Subsidiary[apply(DeepRig_biadjacency_Subsidiary, 1, function(x) any(x == 1)), ]
save(DeepRig_Subsidiary,file='DeepRig_Subsidiary.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_Subsidiary <- DeepRig_Subsidiary[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(没有重叠调控对)
DeepRig_overlap_networks_Subsidiary <- merge(DeepRig_Subsidiary, chip_seq, by = c("TF", "Target"))

#bulk数据
DeepRig_leaf <- read.csv("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_leaf.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(DeepRig_leaf) <- c("TF", "Target","EdgeWeight")
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
DeepRig_leaf$EdgeWeight.zscore <- transform_zscore(DeepRig_leaf$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepRig_biadjacency_leaf <- biadjacency_matrix(DeepRig_leaf$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepRig_biadjacency_leaf <- ifelse(DeepRig_biadjacency_leaf == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepRig_leaf <- DeepRig_leaf[apply(DeepRig_biadjacency_leaf, 1, function(x) any(x == 1)), ]
length(unique(DeepRig_leaf$Target))
length(unique(DeepRig_leaf$TF))
DeepRig_leaf<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_leaf.csv", show_col_types = FALSE)
