setwd("D:/玉米多模态数据的基因调控网络研究/GRNBoost2")
#Guard
library(readxl)
library(dplyr)
GRNBoost2_Guard <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Guard.xlsx")
colnames(GRNBoost2_Guard) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_Guard$EdgeWeight.zscore <- transform_zscore(GRNBoost2_Guard$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_Guard <- biadjacency_matrix(GRNBoost2_Guard$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_Guard <- ifelse(GRNBoost2_biadjacency_Guard == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_Guard <- GRNBoost2_Guard[apply(GRNBoost2_biadjacency_Guard, 1, function(x) any(x == 1)), ]
save(GRNBoost2_Guard,file='GRNBoost2_Guard.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
GRNBoost2_Guard <- GRNBoost2_Guard[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(82重叠调控对)
GRNBoost2_overlap_networks_Guard <- merge(GRNBoost2_Guard, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(107个TFs)
GRNBoost2_unique_TF_Guard <- length(unique(GRNBoost2_Guard$TF))
#统计不重复的 Target数量(21098个Targets)
GRNBoost2_unique_Target_Guard <- length(unique(GRNBoost2_Guard$Target))
#使用超几何分布来验证p值
N <- 2257486 # 背景基因总数(TFxTarget)
K <- 48992   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 82 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.06497081

#Bundle sheath
library(readxl)
library(dplyr)
GRNBoost2_Bun <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Bundle.xlsx")
colnames(GRNBoost2_Bun) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_Bun$EdgeWeight.zscore <- transform_zscore(GRNBoost2_Bun$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_Bun <- biadjacency_matrix(GRNBoost2_Bun$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_Bun <- ifelse(GRNBoost2_biadjacency_Bun == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_Bun <- GRNBoost2_Bun[apply(GRNBoost2_biadjacency_Bun, 1, function(x) any(x == 1)), ]
save(GRNBoost2_Bun,file='GRNBoost2_Bun.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
GRNBoost2_Bun <- GRNBoost2_Bun[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(64重叠调控对)
GRNBoost2_overlap_networks_Bun <- merge(GRNBoost2_Bun, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(99个TFs)
GRNBoost2_unique_TF_Bun <- length(unique(GRNBoost2_Bun$TF))
#统计不重复的 Target数量(17950个Targets)
GRNBoost2_unique_Target_Bun <- length(unique(GRNBoost2_Bun$Target))
#使用超几何分布来验证p值
N <- 1777050 # 背景基因总数(TFxTarget)
K <- 46204   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 64 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 0.9858698
 
#pavement
library(readxl)
library(dplyr)
GRNBoost2_pavement <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/pavement.xlsx")
colnames(GRNBoost2_pavement) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_pavement$EdgeWeight.zscore <- transform_zscore(GRNBoost2_pavement$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_pavement <- biadjacency_matrix(GRNBoost2_pavement$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_pavement <- ifelse(GRNBoost2_biadjacency_pavement == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_pavement <- GRNBoost2_pavement[apply(GRNBoost2_biadjacency_pavement, 1, function(x) any(x == 1)), ]
save(GRNBoost2_pavement,file='GRNBoost2_pavement.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
GRNBoost2_pavement <- GRNBoost2_pavement[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(157重叠调控对)
GRNBoost2_overlap_networks_pavement <- merge(GRNBoost2_pavement, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(112个TFs)
GRNBoost2_unique_TF_pavement <- length(unique(GRNBoost2_pavement$TF))
#统计不重复的 Target数量(22459个Targets)
GRNBoost2_unique_Target_pavement <- length(unique(GRNBoost2_pavement$Target))
#使用超几何分布来验证p值
N <- 2515408 # 背景基因总数(TFxTarget)
K <- 51533   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 157 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 8.76501e-23 

#Mesophyll
library(readxl)
library(dplyr)
GRNBoost2_Mesophyll <- read.csv("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Mesophyll.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(GRNBoost2_Mesophyll) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_Mesophyll$EdgeWeight.zscore <- transform_zscore(GRNBoost2_Mesophyll$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_Mesophyll <- biadjacency_matrix(GRNBoost2_Mesophyll$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_Mesophyll <- ifelse(GRNBoost2_biadjacency_Mesophyll == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_Mesophyll <- GRNBoost2_Mesophyll[apply(GRNBoost2_biadjacency_Mesophyll, 1, function(x) any(x == 1)), ]
save(GRNBoost2_Mesophyll,file='GRNBoost2_Mesophyll.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
GRNBoost2_Mesophyll <- GRNBoost2_Mesophyll[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(146重叠调控对)
GRNBoost2_overlap_networks_Mesophyll <- merge(GRNBoost2_Mesophyll, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(121个TFs)
GRNBoost2_unique_TF_Mesophyll <- length(unique(GRNBoost2_Mesophyll$TF))
#统计不重复的 Target数量(24608个Targets)
GRNBoost2_unique_Target_Mesophyll <- length(unique(GRNBoost2_Mesophyll$Target))
#使用超几何分布来验证p值
N <- 2977568 # 背景基因总数(TFxTarget)
K <- 66291   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 146 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value: 1.223197e-15  

#subsidiary
library(readxl)
library(dplyr)
GRNBoost2_Subsidiary <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/subsidiary.xlsx")
colnames(GRNBoost2_Subsidiary) <- c("TF", "Target","EdgeWeight")
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
GRNBoost2_Subsidiary$EdgeWeight.zscore <- transform_zscore(GRNBoost2_Subsidiary$EdgeWeight)
# 将 biadjacency 矩阵中 z-score 大于阈值的值转换为 1，小于等于阈值的值转换为 0
# 根据 z-score 值进行筛选
GRNBoost2_biadjacency_Subsidiary <- biadjacency_matrix(GRNBoost2_Subsidiary$EdgeWeight.zscore, pvalue.cutoff = 0.05)
GRNBoost2_biadjacency_Subsidiary <- ifelse(GRNBoost2_biadjacency_Subsidiary == TRUE, 1, 0)
# 从原始数据框中保留 biadjacency 矩阵中值为 1 的行
GRNBoost2_Subsidiary <- GRNBoost2_Subsidiary[apply(GRNBoost2_biadjacency_Subsidiary, 1, function(x) any(x == 1)), ]
save(GRNBoost2_Subsidiary,file='GRNBoost2_Subsidiary.Rdata')
#使用chip-seq数据验证GRN
# 保留前两列
GRNBoost2_Subsidiary <- GRNBoost2_Subsidiary[, 1:2]
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 找到重叠的调控对(192重叠调控对)
GRNBoost2_overlap_networks_Subsidiary <- merge(GRNBoost2_Subsidiary, chip_seq, by = c("TF", "Target"))
#统计不重复的 TF 数量(112个TFs)
GRNBoost2_unique_TF_Subsidiary <- length(unique(GRNBoost2_Subsidiary$TF))
#统计不重复的 Target数量(21722个Targets)
GRNBoost2_unique_Target_Subsidiary <- length(unique(GRNBoost2_Subsidiary$Target))
#使用超几何分布来验证p值
N <- 2432864 # 背景基因总数(TFxTarget)
K <- 48367   # GRN推断的调控对
n <- 3174  # 被验证的TFs的chip-seq调控对
x <- 192 #重叠的调控对
# 计算超几何分布的 p 值
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
# 输出结果
cat("P-value:", p_value, "\n")
#P-value:3.774524e-40 