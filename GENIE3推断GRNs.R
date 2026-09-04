setwd("D:/玉米多模态数据的基因调控网络研究/GENIE3")
load("D:/玉米多模态数据的基因调控网络研究/单细胞测序/zeamaysleaf.Rdata")

#提取每个细胞类型的表达矩阵
#Mesophyll
Mesophyll_expression_matrix <- subset(zeamaysleaf, idents = "Mesophyll")@assays$RNA$data
#将Mesophyll_expression_matrix转换为数据框
Mesophyll_expr=data.frame(Mesophyll_expression_matrix)
#转换为矩阵
Mesophyll_expr <- as.matrix(Mesophyll_expr)
#将表达值为0的基因过滤
expr_Mesophyll <- Mesophyll_expr[rowSums(Mesophyll_expr != 0) > 0, ]
dim(expr_Mesophyll)
write.csv(expr_Mesophyll, file = "expr_Mesophyll.csv", row.names = TRUE)
save(expr_Mesophyll,file='expr_Mesophyll.Rdata')
#Guard
Guard_expression_matrix <- subset(zeamaysleaf, idents = "Guard")@assays$RNA$data
#将Guard_expression_matrix转换为数据框
Guard_expr=data.frame(Guard_expression_matrix)
#转换为矩阵
Guard_expr <- as.matrix(Guard_expr)
#将表达值为0的基因过滤
expr_Guard <- Guard_expr[rowSums(Guard_expr != 0) > 0, ]
dim(expr_Guard)
save(expr_Guard,file='expr_Guard.Rdata')
write.csv(expr_Guard, file = "expr_Guard.csv", row.names = TRUE)
#Pavement
Pavement_expression_matrix <- subset(zeamaysleaf, idents = "Pavement")@assays$RNA$data
#将Pavement_expression_matrix转换为数据框
Pavement_expr=data.frame(Pavement_expression_matrix)
#转换为矩阵
Pavement_expr <- as.matrix(Pavement_expr)
#将表达值为0的基因过滤
expr_Pavement <- Pavement_expr[rowSums(Pavement_expr != 0) > 0, ]
dim(Pavement_expr)
dim(expr_Pavement)
save(expr_Pavement,file='expr_Pavement.Rdata')
write.csv(expr_Pavement, file = "expr_Pavement.csv", row.names = TRUE)
#Subsidiary
Subsidiary_expression_matrix <- subset(zeamaysleaf, idents = "Subsidiary")@assays$RNA$data
#将Subsidiary_expression_matrix转换为数据框
Subsidiary_expr=data.frame(Subsidiary_expression_matrix)
#转换为矩阵
Subsidiary_expr <- as.matrix(Subsidiary_expr)
#将表达值为0的基因过滤
expr_Subsidiary <- Subsidiary_expr[rowSums(Subsidiary_expr != 0) > 0, ]
save(expr_Subsidiary,file='expr_Subsidiary.Rdata')
dim(expr_Subsidiary)
save(expr_Subsidiary,file='expr_Subsidiary.Rdata')
write.csv(expr_Subsidiary, file = "expr_Subsidiary.csv", row.names = TRUE)
#Bundle sheath
Bun_expression_matrix <- subset(zeamaysleaf, idents = "Bundle sheath")@assays$RNA$data
#将Bun_expression_matrix转换为数据框
Bun_expr=data.frame(Bun_expression_matrix)
#转换为矩阵
Bun_expr <- as.matrix(Bun_expr)
#将表达值为0的基因过滤
expr_Bun <- Bun_expr[rowSums(Bun_expr != 0) > 0, ]
save(expr_Bun,file='expr_Bun.Rdata')
write.csv(expr_Bun, file = "expr_Bun.csv", row.names = TRUE)
dim(expr_Bun)

#筛选每个细胞类型的转录因子
library(readxl)
#Mesophyll
tf_all<- read_excel("D:/玉米多模态数据的基因调控网络研究/Zeamays_TFs.xlsx")
#将玉米的转录因子转换为一个向量
tf_all <- tf_all$`gene id`
#筛选Mesophyll细胞类型中的转录因子
tf_Mesophyll <- intersect(rownames(expr_Mesophyll), tf_all)
write.csv(tf_Mesophyll, file = "tf_Mesophyll.csv", row.names = TRUE)
save(tf_Mesophyll,file='tf_Mesophyll.Rdata')
#将向量转换为一个数据框格式
tf_Mesophyll_df <- data.frame(TF = tf_Mesophyll)
write.xlsx(tf_Mesophyll_df, file = "tf_Mesophyll.xlsx")
#Guard
tf_Guard <- intersect(rownames(expr_Guard), tf_all)
write.csv(tf_Guard, file = "tf_Guard.csv", row.names = TRUE)
#将向量转换为一个数据框格式
tf_Guard_df <- data.frame(TF = tf_Guard)
write.xlsx(tf_Guard_df, file = "tf_Guard_df.xlsx")
#Pavement
tf_Pavement <- intersect(rownames(expr_Pavement), tf_all)
write.csv(tf_Pavement, file = "tf_Pavement.csv", row.names = TRUE)
#将向量转换为一个数据框格式
tf_Pavement_df <- data.frame(TF = tf_Pavement)
write.xlsx(tf_Pavement_df, file = "tf_Pavement.xlsx")
#Subsidiary
tf_Subsidiary <- intersect(rownames(expr_Subsidiary), tf_all)
write.csv(tf_Subsidiary, file = "tf_Subsidiary.csv", row.names = TRUE)
#将向量转换为一个数据框格式
tf_Subsidiary_df <- data.frame(TF = tf_Subsidiary)
write.xlsx(tf_Subsidiary_df, file = "tf_Subsidiary.xlsx")
save(tf_Subsidiary,file='tf_Subsidiary.Rdata')
#Bundle sheath
tf_Bun <- intersect(rownames(expr_Bun), tf_all)
write.csv(tf_Bun, file = "tf_Bun.csv", row.names = TRUE)
#将向量转换为一个数据框格式
tf_Bun_df <- data.frame(TF = tf_Bun)
write.xlsx(tf_Bun_df, file = "tf_Bun.xlsx")

#使用GENIE3推断每个细胞类型的基因调控网络
library(GENIE3)
#Mesophyll
wm_Mesophyll <- GENIE3(expr_Mesophyll, nCores = 2, regulators = tf_Mesophyll, K = "sqrt")
#将Mesophyll的权重矩阵保存
save(wm_Mesophyll,file='wm_Mesophyll.Rdata')
#创建一个功能函数，将目标中的每个元素标准化z-score，输出一个新的矩阵一个新的矩阵mat.zscore，其中每个元素是对应原矩阵元素的z-score
transform_zscore <- function(mat){
  
  mat.mean <- mean(mat[!is.na(mat)])
  mat.sd <- sd(mat[!is.na(mat)])
  mat.zscore <- (mat - mat.mean)/mat.sd
  return(mat.zscore)
}
wm_Mesophyll.zscore <- transform_zscore(wm_Mesophyll)
#创建一个功能函数，根据z-score矩阵生成生成一个仅包含0和1的矩阵
biadjacency_matrix <- function(mat.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  mat.biadjacency <- as.matrix(mat.zscore > zscore.cutoff) * 1
  return(mat.biadjacency)
}
Mesophyll.biadjacency <- biadjacency_matrix(wm_Mesophyll.zscore)
library(igraph)
#使用igraph可视化基因调控网络
graph_Mesophyll <- graph_from_biadjacency_matrix(Mesophyll.biadjacency)
plot(graph_Mesophyll, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
GENIE3_Mesophyll <- as_edgelist(graph_Mesophyll)
colnames(GENIE3_Mesophyll) <- c("TF", "Target")
GENIE3_Mesophyll <- as.data.frame(GENIE3_Mesophyll)
save(GENIE3_Mesophyll,file='GENIE3_Mesophyll.Rdata')
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
edges_with_weights_Mesophyll <- extract_weights_from_graph(graph_Mesophyll, wm_Mesophyll)
write.xlsx(edges_with_weights_Mesophyll, "edges_with_weights_Mesophyll.xlsx")
#使用CHIP-seq数据验证基因调控网络
linklist_mesophyll <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Mesophyll.xlsx")
# 指定目标转录因子
target_TF <- "Zm00001d015468"
mesophyll_targets <- linklist_mesophyll[linklist_mesophyll$TF == target_TF, "Target"]
mesophyll_targets <- mesophyll_targets$Target
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_mesophyll<- intersect(mesophyll_targets, CHIP_targets)
# 超几何分布的参数
N <- 3036011 # 背景基因总数
K <- 1226438   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 2195 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 4.814566e-235  

#Guard
wm_Guard <- GENIE3(expr_Guard, nCores = 2, regulators = tf_Guard, K = "sqrt")
save(wm_Guard,file='wm_Guard.Rdata')
wm_Guard.zscore <- transform_zscore(wm_Guard)
Guard.biadjacency <- biadjacency_matrix(wm_Guard.zscore)
graph_Guard <- graph_from_biadjacency_matrix(Guard.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Guard <- as_edgelist(graph_Guard)
GENIE3_Guard <- as_edgelist(graph_Guard)
colnames(GENIE3_Guard) <- c("TF", "Target")
GENIE3_Guard <- as.data.frame(GENIE3_Guard)
save(GENIE3_Guard,file='GENIE3_Guard.Rdata')
write.xlsx(linklist_Guard, "linklist_Guard.xlsx")
linklist_Guard <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Guard.xlsx")
#获取调控边的权重值
edges_with_weights_Guard <- extract_weights_from_graph(graph_Guard, wm_Guard)
write.xlsx(edges_with_weights_Guard, "edges_with_weights_Guard.xlsx")
# 指定目标转录因子
target_TF <- "Zm00001d015468"
# 筛选出目标转录因子的靶基因
target_genes <- linklist_Guard[linklist_Guard$TF == target_TF, "Target"]
guard_targets <- target_genes$Target
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_guard<- intersect(guard_targets, CHIP_targets)
# 超几何分布的参数
N <- 2302640 # 背景基因总数
K <- 81908   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 114 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.04219126 

#Pavement
wm_Pavement <- GENIE3(expr_Pavement, nCores = 2, regulators = tf_Pavement, K = "sqrt")
save(wm_Pavement,file='wm_Pavement.Rdata')
wm_Pavement.zscore <- transform_zscore(wm_Pavement)
Pavement.biadjacency <- biadjacency_matrix(wm_Pavement.zscore)
graph_Pavement <- graph_from_biadjacency_matrix(Pavement.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Pavement <- as_edgelist(graph_Pavement)
GENIE3_Pavement <- as_edgelist(graph_Pavement)
colnames(GENIE3_Pavement) <- c("TF", "Target")
GENIE3_Pavement <- as.data.frame(GENIE3_Pavement)
save(GENIE3_Pavement,file='GENIE3_Pavement.Rdata')
write.xlsx(linklist_Pavement, "linklist_Pavement.xlsx")
#获取调控边的权重值
edges_with_weights_Pavement <- extract_weights_from_graph(graph_Pavement, wm_Pavement)
write.xlsx(edges_with_weights_Pavement, "edges_with_weights_Pavement.xlsx")

#Subsidiary
wm_Subsidiary <- GENIE3(expr_Subsidiary, nCores = 2, regulators = tf_Subsidiary, K = "sqrt")
wm_Subsidiary.zscore <- transform_zscore(wm_Subsidiary)
Subsidiary.biadjacency <- biadjacency_matrix(wm_Subsidiary.zscore)
graph_Subsidiary <- graph_from_biadjacency_matrix(Subsidiary.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Subsidiary <- as_edgelist(graph_Subsidiary)
GENIE3_Subsidiary <- as_edgelist(graph_Subsidiary)
colnames(GENIE3_Subsidiary) <- c("TF", "Target")
GENIE3_Subsidiary <- as.data.frame(GENIE3_Subsidiary)
save(GENIE3_Subsidiary,file='GENIE3_Subsidiary.Rdata')
write.xlsx(linklist_Subsidiary, "linklist_Subsidiary.xlsx")
#获取调控边的权重值
edges_with_weights_Subsidiary <- extract_weights_from_graph(graph_Subsidiary, wm_Subsidiary)
write.xlsx(edges_with_weights_Subsidiary, "edges_with_weights_Subsidiary.xlsx")

#Bundle Sheath
wm_Bun <- GENIE3(expr_Bun, nCores = 2, regulators = tf_Bun, K = "sqrt")
save(wm_Bun,file='wm_Bun.Rdata')
wm_Bun.zscore <- transform_zscore(wm_Bun)
Bun.biadjacency <- biadjacency_matrix(wm_Bun.zscore)
graph_Bun <- graph_from_biadjacency_matrix(Bun.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Bun <- as_edgelist(graph_Bun)
GENIE3_Bun <- as_edgelist(graph_Bun)
colnames(GENIE3_Bun) <- c("TF", "Target")
GENIE3_Bun <- as.data.frame(GENIE3_Bun)
save(GENIE3_Bun,file='GENIE3_Bun.Rdata')
write.xlsx(linklist_Bun, "linklist_Bun.xlsx")
#获取调控边的权重值
edges_with_weights_Bun <- extract_weights_from_graph(graph_Bun, wm_Bun)
write.xlsx(edges_with_weights_Bun, "edges_with_weights_Bun.xlsx")

#使用CHIP-seq数据验证基因调控网络
linklist_Bun <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Bun.xlsx")
# 指定目标转录因子
target_TF <- "Zm00001d015468"
# 筛选出目标转录因子的靶基因
bun_targets <- linklist_Bun[linklist_Bun$TF == target_TF, "Target"]
bun_targets <- bun_targets$Target
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_bun<- intersect(bun_targets, CHIP_targets)
# 超几何分布的参数
N <- 1829124 # 背景基因总数
K <- 77957   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 107 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.9954872  

#使用GENIE3的另一种机器学习算法来推断GRN
#Bundle Sheath
set.seed(123)
weightMat_Bun <- GENIE3(expr_Bun, nCores = 2, treeMethod="ET", K="sqrt", regulators = tf_Bun, nTrees=1000)
save(weightMat_Bun,file='weightMat_Bun.Rdata')
weightMat_Bun.zscore <- transform_zscore(weightMat_Bun)
weightMat_Bun.biadjacency <- biadjacency_matrix(weightMat_Bun.zscore)
graph_Bun_ET <- graph_from_biadjacency_matrix(weightMat_Bun.biadjacency)
plot(graph_Guard, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Bun <- as_edgelist(graph_Bun_ET)
write.xlsx(linklist_Bun, "linklist_Bun.xlsx")
#获取调控边的权重值
edges_with_weights_Bun_ET <- extract_weights_from_graph(graph_Bun_ET, weightMat_Bun)
write.xlsx(edges_with_weights_Bun_ET, "edges_with_weights_Bun_ET.xlsx")

#Guard
set.seed(123)
weightMat_Guard <- GENIE3(expr_Guard, nCores = 2, treeMethod="ET", K="sqrt", regulators = tf_Guard, nTrees=1000)
save(weightMat_Guard,file='weightMat_Guard.Rdata')
weightMat_Guard.zscore <- transform_zscore(weightMat_Guard)
weightMat_Guard.biadjacency <- biadjacency_matrix(weightMat_Guard.zscore)
graph_Guard_ET <- graph_from_biadjacency_matrix(weightMat_Guard.biadjacency)
plot(graph_Guard_ET, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Guard <- as_edgelist(graph_Guard_ET)
write.xlsx(linklist_Guard, "linklist_Guard.xlsx")
#获取调控边的权重值
edges_with_weights_Guard_ET <- extract_weights_from_graph(graph_Guard_ET, weightMat_Guard)
write.xlsx(edges_with_weights_Guard_ET, "edges_with_weights_Guard_ET.xlsx")

#Pavement
set.seed(123)
weightMat_Pavement <- GENIE3(expr_Pavement, nCores = 2, treeMethod="ET", K="sqrt", regulators = tf_Pavement, nTrees=1000)
save(weightMat_Pavement,file='weightMat_Pavement.Rdata')
weightMat_Pavement.zscore <- transform_zscore(weightMat_Pavement)
weightMat_Pavement.biadjacency <- biadjacency_matrix(weightMat_Pavement.zscore)
graph_Pavement_ET <- graph_from_biadjacency_matrix(weightMat_Pavement.biadjacency)
plot(graph_Pavement_ET, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Pavement <- as_edgelist(graph_Pavement_ET)
write.xlsx(linklist_Pavement, "linklist_Pavement.xlsx")
#获取调控边的权重值
edges_with_weights_Pavement_ET <- extract_weights_from_graph(graph_Pavement_ET, weightMat_Pavement)
write.xlsx(edges_with_weights_Pavement_ET, "edges_with_weights_Pavement_ET.xlsx")

#Mesophyll
set.seed(123)
weightMat_Mesophyll <- GENIE3(expr_Mesophyll, nCores = 2, treeMethod="ET", K="sqrt", regulators = tf_Mesophyll, nTrees=1000)
save(weightMat_Mesophyll,file='weightMat_Mesophyll.Rdata')
weightMat_Mesophyll.zscore <- transform_zscore(weightMat_Mesophyll)
weightMat_Mesophyll.biadjacency <- biadjacency_matrix(weightMat_Mesophyll.zscore)
graph_Mesophyll_ET <- graph_from_biadjacency_matrix(weightMat_Mesophyll.biadjacency)
plot(graph_Mesophyll_ET, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Mesophyll <- as_edgelist(graph_Mesophyll_ET)
write.xlsx(linklist_Mesophyll, "linklist_Mesophyll.xlsx")
#获取调控边的权重值
edges_with_weights_Mesophyll_ET <- extract_weights_from_graph(graph_Mesophyll_ET, weightMat_Mesophyll)
write.xlsx(edges_with_weights_Mesophyll_ET, "edges_with_weights_Mesophyll_ET.xlsx")

#Subsidiary
set.seed(123)
weightMat_Subsidiary <- GENIE3(tf_Subsidiary, nCores = 2, treeMethod="ET", K="sqrt", regulators = tf_Subsidiary, nTrees=1000)
weightMat_Subsidiary.zscore <- transform_zscore(weightMat_Subsidiary)
weightMat_Subsidiary.biadjacency <- biadjacency_matrix(weightMat_Subsidiary.zscore)
graph_Subsidiary_ET <- graph_from_biadjacency_matrix(weightMat_Subsidiary.biadjacency)
plot(graph_Subsidiary_ET, vertex.size = 5, vertex.label.cex = 0.7)
#获取连接关系
linklist_Subsidiary <- as_edgelist(graph_Subsidiary_ET)
write.xlsx(linklist_Subsidiary, "linklist_Subsidiary.xlsx")
#获取调控边的权重值
edges_with_weights_Subsidiary_ET <- extract_weights_from_graph(graph_Subsidiary_ET, weightMat_Subsidiary)
write.xlsx(edges_with_weights_Subsidiary_ET, "edges_with_weights_Subsidiary_ET.xlsx")

#Bundle sheath
Bun_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3预测的可验证靶基因-ET/Bun.xlsx")
Bun_targets <- Bun_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_Bun <- intersect(Bun_targets, CHIP_targets)

#Guard
guard_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3预测的可验证靶基因-ET/Guard.xlsx")
guard_targets <- guard_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/验证的CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_guard<- intersect(guard_targets, CHIP_targets)

#pavement
pavement_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3预测的可验证靶基因-ET/Pavement.xlsx")
pavement_targets <- pavement_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468的靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_pavement<- intersect(pavement_targets, CHIP_targets)

#Subsidiary
Subsidiary_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3预测的可验证靶基因-ET/Subsidiary.xlsx")
Subsidiary_targets <- Subsidiary_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_Subsidiary<- intersect(Subsidiary_targets, CHIP_targets)

#Mesophyll
Mesophyll_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3预测的可验证靶基因/Mesophyll预测的可验证靶基因.xlsx")
Mesophyll_targets <- Mesophyll_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_Mesophyll<- intersect(Mesophyll_targets, CHIP_targets)

#使用CHIP-seq数据验证基因调控网络
#使用igraph包来获取重叠的调控对
#Mesophyll
genie3_graph_mesophyll <- graph_from_data_frame(linklist_Mesophyll, directed = TRUE)
CHIP_seq_104 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/104个转录因子的靶基因.xlsx")
chipseq_graph <- graph_from_data_frame(CHIP_seq_104, directed = TRUE)
genie3_edges_mesophyll <- as.data.frame(as_edgelist(genie3_graph_mesophyll))
chipseq_edges <- as.data.frame(as_edgelist(chipseq_graph))
colnames(chipseq_edges) <- c("TF", "Target")
colnames(genie3_edges_mesophyll) <- c("TF", "Target")
overlapping_pairs <- merge(genie3_edges_mesophyll, chipseq_edges, by = "Target")
Mesophyll_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Mesophyll推断靶基因.xlsx")

mesophyll_targets <- Mesophyll_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468的靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_mesophyll <- intersect(mesophyll_targets, CHIP_targets)
# 超几何分布的参数
N <- 10906164 # 背景基因总数
K <- 826721   # GRN推断的调控对
n <- 181375  # 被验证的TFs的chip-seq调控对
x <- 7331 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")

#Bundle sheath
Bun_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Bun预测的可验证靶基因.xlsx")
Bun_targets <- Bun_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/验证的CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene 
common_targets_Bun <- intersect(Bun_targets, CHIP_targets)

#guard
guard_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Guard预测的可验证靶基因.xlsx")
guard_targets <- guard_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/验证的CHIP-seq靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_guard<- intersect(guard_targets, CHIP_targets)

#pavement
pavement_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Pavement预测的可验证靶基因.xlsx")
pavement_targets <- pavement_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468的靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_pavement<- intersect(pavement_targets, CHIP_targets)

#subsidiary
subsidiary_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Subsidiary预测的可验证靶基因.xlsx")
subsidiary_targets <- subsidiary_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468的靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_subsidiary<- intersect(subsidiary_targets, CHIP_targets)

#Mesophyll
Mesophyll_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Mesophyll预测的可验证靶基因.xlsx")
Mesophyll_targets <- Mesophyll_targets$gene
CHIP_targets <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/Zm00001d015468的靶基因.xlsx")
CHIP_targets <- CHIP_targets$gene
common_targets_Mesophyll<- intersect(Mesophyll_targets, CHIP_targets)









