setwd("D:/玉米多模态数据的基因调控网络研究/KBoost")
#创建一个功能函数，将目标中的每个元素标准化z-score，输出一个新的矩阵一个新的矩阵mat.zscore，其中每个元素是对应原矩阵元素的z-score
transform_zscore <- function(mat){
  
  mat.mean <- mean(mat[!is.na(mat)])
  mat.sd <- sd(mat[!is.na(mat)])
  mat.zscore <- (mat - mat.mean)/mat.sd
  return(mat.zscore)
}
#Bundle sheath
#将权重矩阵进行转置
expr_Bun_GRN <- t(expr_Bun_GRN)
wm_Bun.zscore <- transform_zscore(expr_Bun_GRN)
#创建一个功能函数，根据z-score矩阵生成生成一个仅包含0和1的矩阵
biadjacency_matrix <- function(mat.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  mat.biadjacency <- as.matrix(mat.zscore > zscore.cutoff) * 1
  return(mat.biadjacency)
}
Bun.biadjacency <- biadjacency_matrix(wm_Bun.zscore)
library(igraph)
library(openxlsx)
#使用igraph可视化基因调控网络
graph_Bun <- graph_from_biadjacency_matrix(Bun.biadjacency)
plot(graph_Mesophyll, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
KBoost_Bun <- as_edgelist(graph_Bun)
colnames(KBoost_Bun) <- c("TF", "Target")
save(KBoost_Bun,file='KBoost_Bun.Rdata')
write.xlsx(KBoost_Bun, "KBoost_Bun.xlsx")
#使用chip-seq数据验证GRN
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(53重叠调控对)
KBoost_overlap_networks_Bun <- merge(KBoost_Bun, chip_seq, by = c("TF", "Target"))
# 将KBoost_Bun转换为数据框
KBoost_Bun <- as.data.frame(KBoost_Bun)
#统计不重复的 TF 数量(86个TFs)
KBoost_unique_TF_Bun <- length(unique(KBoost_Bun$TF))
#统计不重复的 Target数量(18527个Targets)
KBoost_unique_Target_Bun <- length(unique(KBoost_Bun$Target))
#使用超几何分布来验证p值
N <- 1593322 # 背景基因总数(TFxTarget)
K <- 24266   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 53 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value:  0.2681303

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
edges_with_weights_Bun <- extract_weights_from_graph(graph_Bun, expr_Bun_GRN)
write.xlsx(edges_with_weights_Bun, "edges_with_weights_Bun.xlsx")

#Mesophyll
#将权重矩阵进行转置
expr_Mesophyll_GRN <- t(expr_Mesophyll_GRN)
wm_Mesophyll.zscore <- transform_zscore(expr_Mesophyll_GRN)
Mesophyll.biadjacency <- biadjacency_matrix(wm_Mesophyll.zscore)
library(igraph)
#使用igraph可视化基因调控网络
graph_Mesophyll <- graph_from_biadjacency_matrix(Mesophyll.biadjacency)
plot(graph_Mesophyll, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
KBoost_Mesophyll <- as_edgelist(graph_Mesophyll)
colnames(KBoost_Mesophyll) <- c("TF", "Target")
save(KBoost_Mesophyll,file='KBoost_Mesophyll.Rdata')
write.xlsx(KBoost_Mesophyll, "KBoost_Mesophyll.xlsx")
#使用chip-seq数据验证GRN
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(26重叠调控对)
KBoost_overlap_networks_Mesophyll <- merge(KBoost_Mesophyll, chip_seq, by = c("TF", "Target"))
# 将KBoost_Bun转换为数据框
KBoost_Mesophyll <- as.data.frame(KBoost_Mesophyll)
#统计不重复的 TF 数量(109个TFs)
KBoost_unique_TF_Mesophyll <- length(unique(KBoost_Mesophyll$TF))
#统计不重复的 Target数量(24708个Targets)
KBoost_unique_Target_Mesophyll <- length(unique(KBoost_Mesophyll$Target))
#使用超几何分布来验证p值
N <- 2693172 # 背景基因总数(TFxTarget)
K <- 24975   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 26 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value:  0.7626473

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
edges_with_weights_Mesophyll <- extract_weights_from_graph(graph_Mesophyll, expr_Mesophyll_GRN)
write.xlsx(edges_with_weights_Mesophyll, "edges_with_weights_Mesophyll.xlsx")

#Guard
#将权重矩阵进行转置
expr_Guard_GRN <- t(expr_Guard_GRN)
wm_Guard.zscore <- transform_zscore(expr_Guard_GRN)
Guard.biadjacency <- biadjacency_matrix(wm_Guard.zscore)
graph_Guard <- graph_from_biadjacency_matrix(Guard.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
KBoost_Guard <- as_edgelist(graph_Guard)
KBoost_Guard <- as.data.frame(KBoost_Guard)
colnames(KBoost_Guard) <- c("TF", "Target")
save(KBoost_Guard,file='KBoost_Guard.Rdata')
write.xlsx(KBoost_Guard, "KBoost_Guard.xlsx")
#使用chip-seq数据验证GRN
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(没有重叠调控对)
KBoost_overlap_networks_Guard <- merge(KBoost_Guard, chip_seq, by = c("TF", "Target"))

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
edges_with_weights_Guard <- extract_weights_from_graph(graph_Guard, expr_Guard_GRN)
write.xlsx(edges_with_weights_Guard, "edges_with_weights_Guard.xlsx")


#Pavement
expr_Pavement_GRN <- t(expr_Pavement_GRN)
wm_Pavement.zscore <- transform_zscore(expr_Pavement_GRN)
Pavement.biadjacency <- biadjacency_matrix(wm_Pavement.zscore)
graph_Pavement <- graph_from_biadjacency_matrix(Pavement.biadjacency)
#获取连接关系
KBoost_Pavement <- as_edgelist(graph_Pavement)
KBoost_Pavement <- as.data.frame(KBoost_Pavement)
colnames(KBoost_Pavement) <- c("TF", "Target")
save(KBoost_Pavement,file='KBoost_Pavement.Rdata')
write.xlsx(KBoost_Pavement, "KBoost_Pavement.xlsx")
#使用chip-seq数据验证GRN
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(20重叠调控对)
KBoost_overlap_networks_Pavement <- merge(KBoost_Pavement, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(107个TFs)
KBoost_unique_TF_Pavement <- length(unique(KBoost_Pavement$TF))
#统计不重复的 Target数量(22910个Targets)
KBoost_unique_Target_Pavement <- length(unique(KBoost_Pavement$Target))
#使用超几何分布来验证p值
N <- 2451370 # 背景基因总数(TFxTarget)
K <- 23469   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 20 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value:  0.9818303

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
edges_with_weights_Pavement <- extract_weights_from_graph(graph_Pavement, expr_Pavement_GRN)
write.xlsx(edges_with_weights_Pavement, "edges_with_weights_Pavement.xlsx")

#Subsidiary
expr_Subsidiary_GRN <- t(expr_Subsidiary_GRN)
wm_Subsidiary.zscore <- transform_zscore(expr_Subsidiary_GRN)
Subsidiary.biadjacency <- biadjacency_matrix(wm_Subsidiary.zscore)
graph_Subsidiary <- graph_from_biadjacency_matrix(Subsidiary.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
KBoost_Subsidiary <- as_edgelist(graph_Subsidiary)
KBoost_Subsidiary <- as.data.frame(KBoost_Subsidiary)
colnames(KBoost_Subsidiary) <- c("TF", "Target")
save(KBoost_Subsidiary,file='KBoost_Subsidiary.Rdata')
write.xlsx(KBoost_Subsidiary, "KBoost_Subsidiary.xlsx")
#使用chip-seq数据验证GRN
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(43重叠调控对)
KBoost_overlap_networks_Subsidiary <- merge(KBoost_Subsidiary, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(105个TFs)
KBoost_unique_TF_Subsidiary <- length(unique(KBoost_Subsidiary$TF))
#统计不重复的 Target数量(22381个Targets)
KBoost_unique_Target_Subsidiary <- length(unique(KBoost_Subsidiary$Target))
#使用超几何分布来验证p值
N <- 2350005 # 背景基因总数(TFxTarget)
K <- 22639   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 43 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value:  0.01896779 

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
edges_with_weights_Subsidiary <- extract_weights_from_graph(graph_Subsidiary, expr_Subsidiary_GRN)
write.xlsx(edges_with_weights_Subsidiary, "edges_with_weights_Subsidiary.xlsx")
