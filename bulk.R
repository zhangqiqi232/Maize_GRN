setwd("D:/玉米多模态数据的基因调控网络研究/DeepSEM")
#bulk数据网络评估
#DeepSEM
#读取tsv文件
library(vroom)
DeepSEM_leaf<- vroom("D:/玉米多模态数据的基因调控网络研究/DeepSEM/leaf.tsv", show_col_types = FALSE)
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
DeepSEM_leaf$EdgeWeight.zscore <- transform_zscore(DeepSEM_leaf$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
DeepSEM_biadjacency_leaf <- biadjacency_matrix(DeepSEM_leaf$EdgeWeight.zscore, pvalue.cutoff = 0.05)
DeepSEM_biadjacency_leaf <- ifelse(DeepSEM_biadjacency_leaf == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
DeepSEM_leaf <- DeepSEM_leaf[apply(DeepSEM_biadjacency_leaf, 1, function(x) any(x == 1)), ]
save(DeepSEM_leaf,file='DeepSEM_leaf.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepSEM_leaf <- DeepSEM_leaf[, 1:2]
#读取chip-seq数据
library(readxl)
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(11527重叠调控对)
DeepSEM_overlap_networks_leaf <- merge(DeepSEM_leaf, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(1527个TFs)
DeepSEM_unique_TF_leaf <- length(unique(DeepSEM_leaf$TF))
#统计不重复的 Target数量(24675个Targets)
DeepSEM_unique_Target_leaf <- length(unique(DeepSEM_leaf$Target))
# 提取重叠调控对中的所有 TF
DeepSEM_overlap_TFs <- unique(DeepSEM_overlap_networks_leaf$TF)
# 从 chip-seq 数据中筛选出这些 TF 的所有调控对（22296）
DeepSEM_TF_chipseq <- subset(chip_seq, TF %in% DeepSEM_overlap_TFs)
#使用超几何分布来验证p值
N <- 37678725 # 背景基因总数(TFxTarget)
K <- 1988844   # GRN推断的调控对
n <- 222946  # 被验证的TFs的chip-seq调控对
x <- 11527 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.9892933  

#DeepRig
#读取tsv文件
library(vroom)
DeepRig_leaf <- DeepRig_leaf_processed
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
save(DeepRig_leaf,file='DeepRig_leaf.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
DeepRig_leaf <- DeepRig_leaf[, 1:2]
colnames(DeepRig_leaf) <- c("TF", "Target")
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(18519重叠调控对)
DeepRig_overlap_networks_leaf <- merge(DeepRig_leaf, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(1527个TFs)
DeepRig_unique_TF_leaf <- length(unique(DeepRig_leaf$TF))
#统计不重复的 Target数量(1527个Targets)
DeepRig_unique_Target_leaf <- length(unique(DeepRig_leaf$Target))
# 提取重叠调控对中的所有 TF
DeepRig_overlap_TFs <- unique(DeepRig_overlap_networks_leaf$TF)
# 从 chip-seq 数据中筛选出这些 TF 的所有调控对（222967）
DeepRig_TF_chipseq <- subset(chip_seq, TF %in% DeepRig_overlap_TFs)
#使用超几何分布来验证p值
N <- 2331729 # 背景基因总数(TFxTarget)
K <- 2167411   # GRN推断的调控对
n <- 222967  # 被验证的TFs的chip-seq调控对
x <- 18519 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 1


#KBoost
#将权重矩阵进行转置
expr_leaf <- t(leaf_GRN)
wm_leaf.zscore <- transform_zscore(expr_leaf)
#创建一个功能函数，根据z-score矩阵生成生成一个仅包含0和1的矩阵
biadjacency_matrix <- function(mat.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  mat.biadjacency <- as.matrix(mat.zscore > zscore.cutoff) * 1
  return(mat.biadjacency)
}
leaf.biadjacency <- biadjacency_matrix(wm_leaf.zscore)
library(igraph)
library(openxlsx)
#使用igraph可视化基因调控网络
graph_leaf <- graph_from_biadjacency_matrix(leaf.biadjacency)
plot(graph_Mesophyll, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
KBoost_leaf <- as_edgelist(graph_leaf)
KBoost_leaf <- as.data.frame(KBoost_leaf)
colnames(KBoost_leaf) <- c("TF", "Target")
save(KBoost_leaf,file='KBoost_leaf.Rdata')
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(415重叠调控对)
KBoost_overlap_networks_leaf <- merge(KBoost_leaf, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(962个TFs)
KBoost_unique_TF_leaf <- length(unique(KBoost_leaf$TF))
#统计不重复的 Target数量(24053个Targets)
KBoost_unique_Target_leaf <- length(unique(KBoost_leaf$Target))
# 提取重叠调控对中的所有 TF
KBoost_overlap_TFs <- unique(KBoost_overlap_networks_leaf$TF)
# 从 chip-seq 数据中筛选出这些 TF 的所有调控对（171218）
KBoost_TF_chipseq <- subset(chip_seq, TF %in% KBoost_overlap_TFs)
#使用超几何分布来验证p值
N <- 23138986 # 背景基因总数(TFxTarget)
K <- 28494   # GRN推断的调控对
n <- 171218  # 被验证的TFs的chip-seq调控对
x <- 415 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 6.931347e-36

# 提取调控边及其权重
extract_weights_from_graph <- function(graph, weight_matrix) {
  # 获取调控边
  edges <- as_edgelist(graph)
  
  # 提取权重
  weights <- sapply(1:nrow(edges), function(i) {
    tf <- edges[i, 1]
    target <- edges[i, 2]
    
    # 使用match来查找行名和列名对应的索引
    weight <- weight_matrix[match(tf, rownames(weight_matrix)), match(target, colnames(weight_matrix))]
    return(weight)
  })
  
  # 创建一个新的数据框
  edges_with_weights <- cbind(edges, Weight = weights)
  return(as.data.frame(edges_with_weights))
}
# 运行示例（此处为示例，需使用实际的 weight_matrix）
edges_with_weights_leaf <- extract_weights_from_graph(graph_leaf, expr_leaf)
write.xlsx(edges_with_weights_leaf, "edges_with_weights_leaf.xlsx")

#GENIE3
transform_zscore <- function(mat){
  
  mat.mean <- mean(mat[!is.na(mat)])
  mat.sd <- sd(mat[!is.na(mat)])
  mat.zscore <- (mat - mat.mean)/mat.sd
  return(mat.zscore)
}
wm_leaf.zscore <- transform_zscore(wm_leaf)
#创建一个功能函数，根据z-score矩阵生成生成一个仅包含0和1的矩阵
biadjacency_matrix <- function(mat.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  mat.biadjacency <- as.matrix(mat.zscore > zscore.cutoff) * 1
  return(mat.biadjacency)
}
leaf.biadjacency <- biadjacency_matrix(wm_leaf.zscore)
library(igraph)
#使用igraph可视化基因调控网络
graph_leaf <- graph_from_biadjacency_matrix(leaf.biadjacency)
plot(graph_Mesophyll, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
GENIE3_leaf <- as_edgelist(graph_leaf)
GENIE3_leaf <- as.data.frame(GENIE3_leaf)
colnames(GENIE3_leaf) <- c("TF", "Target")
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(18468重叠调控对)
GENIE3_overlap_networks_leaf <- merge(GENIE3_leaf, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(1526个TFs)
GENIE3_unique_TF_leaf <- length(unique(GENIE3_leaf$TF))
#统计不重复的 Target数量(24675个Targets)
GENIE3_unique_Target_leaf <- length(unique(GENIE3_leaf$Target))
# 提取重叠调控对中的所有 TF
GENIE3_overlap_TFs <- unique(GENIE3_overlap_networks_leaf$TF)
# 从 chip-seq 数据中筛选出这些 TF 的所有调控对（222967）
GENIE3_TF_chipseq <- subset(chip_seq, TF %in% GENIE3_overlap_TFs)
#使用超几何分布来验证p值
N <- 37654050 # 背景基因总数(TFxTarget)
K <- 1480855   # GRN推断的调控对
n <- 222967  # 被验证的TFs的chip-seq调控对
x <- 18468 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0
library(openxlsx)
write.xlsx(linklist_Mesophyll, "linklist_Mesophyll.xlsx")
# 提取调控边及其权重
extract_weights_from_graph <- function(graph, weight_matrix) {
  # 获取调控边
  edges <- as_edgelist(graph)
  
  # 提取权重
  weights <- sapply(1:nrow(edges), function(i) {
    tf <- edges[i, 1]
    target <- edges[i, 2]
    
    # 使用match来查找行名和列名对应的索引
    weight <- weight_matrix[match(tf, rownames(weight_matrix)), match(target, colnames(weight_matrix))]
    return(weight)
  })
  
  # 创建一个新的数据框
  edges_with_weights <- cbind(edges, Weight = weights)
  return(as.data.frame(edges_with_weights))
}
# 运行示例（此处为示例，需使用实际的 weight_matrix）
edges_with_weights_leaf <- extract_weights_from_graph(graph_leaf, wm_leaf)
write.xlsx(edges_with_weights_leaf, "edges_with_weights_leaf.xlsx")

#GRNBoost2
GRNBoost2_leaf <- read.csv("D:/玉米多模态数据的基因调控网络研究/玉米叶片bulk/推断GRN/GRNBoost2/network_leaf_output.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(GRNBoost2_leaf) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_leaf$EdgeWeight.zscore <- transform_zscore(GRNBoost2_leaf$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_Guard <- biadjacency_matrix(GRNBoost2_leaf$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_Guard <- ifelse(GRNBoost2_biadjacency_Guard == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_leaf <- GRNBoost2_leaf[apply(GRNBoost2_biadjacency_Guard, 1, function(x) any(x == 1)), ]
save(GRNBoost2_leaf,file='GRNBoost2_leaf.Rdata')
# 保留前两列
GRNBoost2_leaf <- GRNBoost2_leaf[, 1:2]
#读取chip-seq数据
library(readxl)
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(3212重叠调控对)
GRNBoost2_overlap_networks_leaf <- merge(GRNBoost2_leaf, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(1476个TFs)
GRNBoost2_unique_TF_leaf <- length(unique(GRNBoost2_leaf$TF))
#统计不重复的 Target数量(23978个Targets)
GRNBoost2_unique_Target_leaf <- length(unique(GRNBoost2_leaf$Target))
# 提取重叠调控对中的所有 TF
GRNBoost2_overlap_TFs <- unique(GRNBoost2_overlap_networks_leaf$TF)
# 从 chip-seq 数据中筛选出这些 TF 的所有调控对（202964）
GRNBoost2_TF_chipseq <- subset(chip_seq, TF %in% GRNBoost2_overlap_TFs)
#使用超几何分布来验证p值
N <- 35391528 # 背景基因总数(TFxTarget)
K <- 175576   # GRN推断的调控对
n <- 202964  # 被验证的TFs的chip-seq调控对
x <- 3212 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0
