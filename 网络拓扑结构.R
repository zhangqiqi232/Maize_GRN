library(dplyr)
library(readxl)
#聚类系数
RF1 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Bun.xlsx")
RF2 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Guard.xlsx")
RF3 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Mesophyll.xlsx")
RF4 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Pavement.xlsx")
RF5 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3推断的调控关系/linklist_Subsidiary.xlsx")
# 合并这五个数据框
RF <- bind_rows(RF1, RF2, RF3, RF4, RF5)
library(igraph)
# 1. 构建预测网络
# 从列表创建有向图
predicted_network <- graph_from_data_frame(RF, directed = TRUE)
# 计算预测网络的平均聚类系数
predicted_clustering <- transitivity(predicted_network, type = "average")
cat("预测网络的聚类系数:", predicted_clustering, "\n")
#预测网络的聚类系数:  0.6039737
# 2. 随机化并计算聚类系数
set.seed(123)  # 设置随机种子保证结果可复现
random_clustering_list <- numeric(100)  # 用于存储随机网络的聚类系数
for (i in 1:100 ) {
  # 随机生成网络，保持节点数和边数不变
  random_network <- erdos.renyi.game(
    n = vcount(predicted_network),       # 节点数量与预测网络一致
    p.or.m = ecount(predicted_network), # 边数量与预测网络一致
    type = "gnm",                       # 固定边数模型
    directed = TRUE                     # 保持网络有向性
  )
  
  # 计算随机网络的聚类系数
  random_clustering_list[i] <- transitivity(random_network, type = "average")
}
# 3. 统计随机网络聚类系数大于预测网络聚类系数的次数
greater_count <- sum(random_clustering_list > predicted_clustering)
p_value <- greater_count / 100  # 计算 p 值
# 输出结果
cat("随机网络中聚类系数大于预测网络的数量:", greater_count, "\n")
cat("p 值:", p_value, "\n")

ET1 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3-ET调控关系/linklist_Bun.xlsx")
ET2 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3-ET调控关系/linklist_Guard.xlsx")
ET3 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3-ET调控关系/linklist_Mesophyll.xlsx")
ET4 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3-ET调控关系/linklist_Pavement.xlsx")
ET5 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/GENIE3-ET调控关系/linklist_Subsidiary.xlsx")
# 合并这五个数据框
ET <- bind_rows(ET1, ET2, ET3, ET4, ET5)
# 1. 构建预测网络
# 从列表创建有向图
predicted_network_ET <- graph_from_data_frame(ET, directed = TRUE)
# 计算预测网络的平均聚类系数
predicted_clustering_ET <- transitivity(predicted_network_ET, type = "average")
cat("预测网络的平均聚类系数:", predicted_clustering_ET, "\n")
#预测网络的平均聚类系数:  0.5775238
# 2. 随机化并计算聚类系数
set.seed(123)  # 设置随机种子保证结果可复现
random_clustering_list_ET <- numeric(100)  # 用于存储随机网络的聚类系数
for (i in 1:100 ) {
  # 随机生成网络，保持节点数和边数不变
  random_network_ET <- erdos.renyi.game(
    n = vcount(predicted_network_ET),       # 节点数量与预测网络一致
    p.or.m = ecount(predicted_network_ET), # 边数量与预测网络一致
    type = "gnm",                       # 固定边数模型
    directed = TRUE                     # 保持网络有向性
  )
  
  # 计算随机网络的聚类系数
  random_clustering_list_ET[i] <- transitivity(random_network_ET, type = "average")
}
# 3. 统计随机网络聚类系数大于预测网络聚类系数的次数
greater_count_ET <- sum(random_clustering_list_ET > predicted_clustering_ET)
p_value <- greater_count_ET / 100  # 计算 p 值
# 输出结果
cat("随机网络中聚类系数大于预测网络的数量:", greater_count_ET, "\n")
cat("p 值:", p_value, "\n")

GRNBoost1 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Bundle.xlsx")
GRNBoost2 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Guard.xlsx")
GRNBoost3 <- read.csv("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Mesophyll.csv", header = TRUE, stringsAsFactors = FALSE)
GRNBoost4 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/pavement.xlsx")
GRNBoost5 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/subsidiary.xlsx")
# 合并这五个数据框
GRNBoost <- bind_rows(GRNBoost1, GRNBoost2, GRNBoost3, GRNBoost4, GRNBoost5)
# 删除第三列
GRNBoost <- GRNBoost[, -3]
# 1. 构建预测网络
# 从列表创建有向图
predicted_network_GRNBoost <- graph_from_data_frame(GRNBoost, directed = TRUE)
# 计算预测网络的平均聚类系数
predicted_clustering_GRNBoost <- transitivity(predicted_network_GRNBoost, type = "average")
cat("预测网络的平均聚类系数:", predicted_clustering_GRNBoost, "\n")
#预测网络的平均聚类系数:   0.0113282 
# 2. 随机化并计算聚类系数
set.seed(123)  # 设置随机种子保证结果可复现
random_clustering_list_GRNBoost <- numeric(100)  # 用于存储随机网络的聚类系数
for (i in 1:100 ) {
  # 随机生成网络，保持节点数和边数不变
  random_network_GRNBoost <- erdos.renyi.game(
    n = vcount(predicted_network_GRNBoost),       # 节点数量与预测网络一致
    p.or.m = ecount(predicted_network_GRNBoost), # 边数量与预测网络一致
    type = "gnm",                       # 固定边数模型
    directed = TRUE                     # 保持网络有向性
  )
  
  # 计算随机网络的聚类系数
  random_clustering_list_GRNBoost[i] <- transitivity(random_network_GRNBoost, type = "average")
}
# 3. 统计随机网络聚类系数大于预测网络聚类系数的次数
greater_count_GRNBoost <- sum(random_clustering_list_GRNBoost > predicted_clustering_GRNBoost)
p_value <- greater_count_GRNBoost / 100  # 计算 p 值
# 输出结果
cat("随机网络中聚类系数大于预测网络的数量:", greater_count_GRNBoost, "\n")
cat("p 值:", p_value, "\n")

#平均路径长度
# 计算预测网络的平均路径长度
predicted_avg_path_length_RF <- average.path.length(predicted_network, directed = TRUE)
cat("预测网络的平均路径长度:", predicted_avg_path_length_RF, "\n")
#预测网络的平均路径长度: 2.600353
# 设置随机种子保证结果可复现
set.seed(123)
# 用于存储随机网络的平均路径长度
random_avg_path_length_list <- numeric(100)

set.seed(123)  # 设置随机种子保证结果可复现
random_avg_path_length_list <- numeric(100)  # 用于存储随机网络的聚类系数
for (i in 1:100 ) {
  # 随机生成网络，保持节点数和边数不变
  random_network <- erdos.renyi.game(
    n = vcount(predicted_network),       # 节点数量与预测网络一致
    p.or.m = ecount(predicted_network), # 边数量与预测网络一致
    type = "gnm",                       # 固定边数模型
    directed = TRUE                     # 保持网络有向性
  )
  
  # 计算随机网络的聚类系数
  random_avg_path_length_list[i] <- average.path.length(random_network, directed = TRUE)
}
# 统计随机网络平均路径长度大于预测网络的次数
greater_count <- sum(random_avg_path_length_list > predicted_avg_path_length_RF)
p_value <- greater_count / 100  # 计算 p 值

# 输出结果
cat("随机网络中平均路径长度大于预测网络的数量:", greater_count, "\n")
cat("p 值:", p_value, "\n")


#RF
# 1. 计算预测网络的特性
predicted_clustering <- transitivity(predicted_network, type = "average")  # 聚类系数
predicted_avg_path_length <- average.path.length(predicted_network, directed = TRUE)  # 平均路径长度
predicted_density <- edge_density(predicted_network)  # 网络密度
cat("预测网络的特性:\n")
cat("聚类系数:", predicted_clustering, "\n")
cat("平均路径长度:", predicted_avg_path_length, "\n")
cat("网络密度:", predicted_density, "\n\n")
#聚类系数: 0.6039737
#平均路径长度: 2.600353
#网络密度: 0.0007607521 
# 设置随机种子保证结果可复现
set.seed(123)

# 初始化存储随机网络的特性
random_clustering_list <- numeric(100)
random_avg_path_length_list <- numeric(100)
random_density_list <- numeric(100)

# 2. 生成随机网络并计算特性
for (i in 1:100) {
  # 随机生成网络，保持节点数和边数不变
  random_network <- sample_gnm(
    n = vcount(predicted_network),       # 节点数量
    m = ecount(predicted_network),       # 边数量
    directed = TRUE                      # 有向网络
  )
  
  # 计算随机网络的特性
  random_clustering_list[i] <- transitivity(random_network, type = "average")
  random_avg_path_length_list[i] <- average.path.length(random_network, directed = TRUE)
  random_density_list[i] <- edge_density(random_network)
}

# 3. 比较预测网络与随机网络
# 聚类系数
greater_clustering <- sum(random_clustering_list > predicted_clustering)
p_value_clustering <- greater_clustering / 100

# 平均路径长度
greater_avg_path_length <- sum(random_avg_path_length_list > predicted_avg_path_length)
p_value_avg_path_length <- greater_avg_path_length / 100

# 网络密度
greater_density <- sum(random_density_list > predicted_density)
p_value_density <- greater_density / 100
# 4. 输出结果
cat("比较结果:\n")
cat("随机网络中聚类系数大于预测网络的数量:", greater_clustering, "\n")
cat("聚类系数的 p 值:", p_value_clustering, "\n\n")

cat("随机网络中平均路径长度大于预测网络的数量:", greater_avg_path_length, "\n")
cat("平均路径长度的 p 值:", p_value_avg_path_length, "\n\n")

cat("随机网络中网络密度大于预测网络的数量:", greater_density, "\n")
cat("网络密度的 p 值:", p_value_density, "\n\n")

#判断预测网络是否为小世界网络
# 4. 比较预测网络与随机网络的特性
cat("\n随机网络聚类系数均值:", mean(random_clustering_list), "\n")
cat("随机网络平均路径长度均值:", mean(random_avg_path_length_list), "\n")
# 判断是否为小世界网络
is_small_world <- predicted_clustering > mean(random_clustering_list) &&
  abs(predicted_avg_path_length - mean(random_avg_path_length_list)) < 0.1 * mean(random_avg_path_length_list)

if (is_small_world) {
  cat("\n预测网络具有小世界网络特性。\n")
} else {
  cat("\n预测网络不具有小世界网络特性。\n")
}

#ET
# 1. 计算预测网络的特性
# 从数据框创建预测网络（有向图）
predicted_network_ET <- graph_from_data_frame(ET, directed = TRUE)
predicted_clustering_ET <- transitivity(predicted_network_ET, type = "average")  # 聚类系数
predicted_avg_path_length_ET <- average.path.length(predicted_network_ET, directed = TRUE)  # 平均路径长度
predicted_density_ET <- edge_density(predicted_network_ET)  # 网络密度
cat("预测网络的特性:\n")
cat("聚类系数:", predicted_clustering_ET, "\n")
cat("平均路径长度:", predicted_avg_path_length_ET, "\n")
cat("网络密度:", predicted_density_ET, "\n\n")
#聚类系数: 0.5775238
#平均路径长度: 2.686938 
#网络密度:  0.0006906535 
# 设置随机种子保证结果可复现
set.seed(123)

# 初始化存储随机网络的特性
random_clustering_list_ET <- numeric(100)
random_avg_path_length_list_ET <- numeric(100)
random_density_list_ET <- numeric(100)

# 2. 生成随机网络并计算特性
for (i in 1:100) {
  # 随机生成网络，保持节点数和边数不变
  random_network_ET <- sample_gnm(
    n = vcount(predicted_network_ET),       # 节点数量
    m = ecount(predicted_network_ET),       # 边数量
    directed = TRUE                      # 有向网络
  )
  
  # 计算随机网络的特性
  random_clustering_list_ET[i] <- transitivity(random_network_ET, type = "average")
  random_avg_path_length_list_ET[i] <- average.path.length(random_network_ET, directed = TRUE)
  random_density_list_ET[i] <- edge_density(random_network_ET)
}

# 3. 比较预测网络与随机网络
# 聚类系数
greater_clustering_ET <- sum(random_clustering_list_ET > predicted_clustering_ET)
p_value_clustering_ET <- greater_clustering_ET / 100

# 平均路径长度
greater_avg_path_length_ET <- sum(random_avg_path_length_list_ET > predicted_avg_path_length_ET)
p_value_avg_path_length_ET <- greater_avg_path_length_ET / 100

# 网络密度
greater_density_ET <- sum(random_density_list_ET > predicted_density_ET)
p_value_density_ET <- greater_density_ET / 100
# 4. 输出结果
cat("比较结果:\n")
cat("随机网络中聚类系数大于预测网络的数量:", greater_clustering_ET, "\n")
cat("聚类系数的 p 值:", p_value_clustering_ET, "\n\n")

cat("随机网络中平均路径长度大于预测网络的数量:", greater_avg_path_length_ET, "\n")
cat("平均路径长度的 p 值:", p_value_avg_path_length_ET, "\n\n")

cat("随机网络中网络密度大于预测网络的数量:", greater_density_ET, "\n")
cat("网络密度的 p 值:", p_value_density_ET, "\n\n")

#GRNBoost2
# 1. 计算预测网络的特性
# 从数据框创建预测网络（有向图）
predicted_network_GRNBoost <- graph_from_data_frame(GRNBoost, directed = TRUE)
predicted_clustering_GRNBoost <- transitivity(predicted_network_GRNBoost, type = "average")  # 聚类系数
predicted_avg_path_length_GRNBoost <- average.path.length(predicted_network_GRNBoost, directed = TRUE)  # 平均路径长度
predicted_density_GRNBoost <- edge_density(predicted_network_GRNBoost)  # 网络密度
cat("预测网络的特性:\n")
cat("聚类系数:", predicted_clustering_GRNBoost, "\n")
cat("平均路径长度:", predicted_avg_path_length_GRNBoost, "\n")
cat("网络密度:", predicted_density_GRNBoost, "\n\n")
#聚类系数: 0.9380575 
#平均路径长度: 1.355252 
#网络密度:  0.007695095  
# 设置随机种子保证结果可复现
set.seed(123)

# 初始化存储随机网络的特性
random_clustering_list_GRNBoost <- numeric(100)
random_avg_path_length_list_GRNBoost <- numeric(100)
random_density_list_GRNBoost <- numeric(100)

# 2. 生成随机网络并计算特性
for (i in 1:100) {
  # 随机生成网络，保持节点数和边数不变
  random_network_GRNBoost <- sample_gnm(
    n = vcount(predicted_network_GRNBoost),       # 节点数量
    m = ecount(predicted_network_GRNBoost),       # 边数量
    directed = TRUE                      # 有向网络
  )
  
  # 计算随机网络的特性
  random_clustering_list_GRNBoost[i] <- transitivity(random_network_GRNBoost, type = "average")
  random_avg_path_length_list_GRNBoost[i] <- average.path.length(random_network_GRNBoost, directed = TRUE)
  random_density_list_GRNBoost[i] <- edge_density(random_network_GRNBoost)
}
# 3. 比较预测网络与随机网络
# 聚类系数
greater_clustering_GRNBoost <- sum(random_clustering_list_GRNBoost > predicted_clustering_GRNBoost)
p_value_clustering_GRNBoost <- greater_clustering_GRNBoost / 100

# 平均路径长度
greater_avg_path_length_GRNBoost <- sum(random_avg_path_length_list_GRNBoost > predicted_avg_path_length_GRNBoost)
p_value_avg_path_length_GRNBoost <- greater_avg_path_length_GRNBoost / 100

# 网络密度
greater_density_GRNBoost <- sum(random_density_list_GRNBoost > predicted_density_GRNBoost)
p_value_density_GRNBoost <- greater_density_GRNBoost/ 100
# 4. 输出结果
cat("比较结果:\n")
cat("随机网络中聚类系数大于预测网络的数量:", greater_clustering_GRNBoost, "\n")
cat("聚类系数的 p 值:", greater_clustering_GRNBoost, "\n\n")

cat("随机网络中平均路径长度大于预测网络的数量:", greater_avg_path_length_GRNBoost, "\n")
cat("平均路径长度的 p 值:", greater_avg_path_length_GRNBoost, "\n\n")

cat("随机网络中网络密度大于预测网络的数量:", greater_density_GRNBoost, "\n")
cat("网络密度的 p 值:", greater_density_GRNBoost, "\n\n")
