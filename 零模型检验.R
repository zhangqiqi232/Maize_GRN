#单细胞分辨率下零模型评估
#GENIE3
# 真实网络（chip-seq）
library(readxl)
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
colnames(chip_seq)[1:2] <- c("TF", "Target")

# 五个细胞类型GRN列表（改名为 GRNBoost2_list）
GENIE3_list <- list(
  Bun = GENIE3_Bun,
  Guard = GENIE3_Guard,
  Mesophyll = GENIE3_Mesophyll,
  Pavement = GENIE3_Pavement,
  Subsidiary = GENIE3_Subsidiary
)

# 统一列名
for (name in names(GENIE3_list)) {
  colnames(GENIE3_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(GENIE3_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = GENIE3_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 显示图
p_line

#GRNBoost2
# 加载依赖包（仅保留必要包）
library(tidyverse)
# ------------------------------------------------------------------------------
# 1. 加载GRNBoost2算法预测的GRN
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_Bun.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_Bun.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_Guard.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_Mesophyll.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_pavement.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/GRNBoost2_Subsidiary.Rdata")
# ------------------------------------------------------------------------------
# 真实网络（chip-seq）
colnames(chip_seq)[1:2] <- c("TF", "Target")

# 五个细胞类型GRN列表（改名为 GRNBoost2_list）
GRNBoost2_list <- list(
  Bun = GRNBoost2_Bun,
  Guard = GRNBoost2_Guard,
  Mesophyll = GRNBoost2_Mesophyll,
  Pavement = GRNBoost2_pavement,
  Subsidiary = GRNBoost2_Subsidiary
)

# 统一列名
for (name in names(GRNBoost2_list)) {
  colnames(GRNBoost2_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(GRNBoost2_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = GRNBoost2_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 显示图
p_line



#KBoost
# 1. 加载GRNBoost2算法预测的GRN
load("D:/玉米多模态数据的基因调控网络研究/KBoost/KBoost_Bun.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/KBoost/KBoost_Guard.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/KBoost/KBoost_Mesophyll.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/KBoost/KBoost_Subsidiary.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/KBoost/KBoost_Pavement.Rdata")
# ------------------------------------------------------------------------------
# 真实网络（chip-seq）
colnames(chip_seq)[1:2] <- c("TF", "Target")

# 五个细胞类型GRN列表（改名为 GRNBoost2_list）
KBoost_list <- list(
  Bun = KBoost_Bun,
  Guard = KBoost_Guard,
  Mesophyll = KBoost_Mesophyll,
  Pavement = KBoost_Pavement,
  Subsidiary = KBoost_Subsidiary
)

# 统一列名
for (name in names(KBoost_list)) {
  colnames(KBoost_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(KBoost_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = KBoost_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 显示图
p_line


#DeepSEM
# 1. 加载DeepSEM算法预测的GRN
load("D:/玉米多模态数据的基因调控网络研究/DeepSEM/DeepSEM_Bun.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepSEM/DeepSEM_Guard.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepSEM/DeepSEM_Mesophyll.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepSEM/DeepSEM_Pavement.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepSEM/DeepSEM_Subsidiary.Rdata")
# ------------------------------------------------------------------------------
# 真实网络（chip-seq）
colnames(chip_seq)[1:2] <- c("TF", "Target")

# 五个细胞类型GRN列表（改名为 GRNBoost2_list）
DeepSEM_list <- list(
  Bun = DeepSEM_Bun,
  Guard = DeepSEM_Guard,
  Mesophyll = DeepSEM_Mesophyll,
  Pavement = DeepSEM_Pavement,
  Subsidiary = DeepSEM_Subsidiary
)

# 统一列名
for (name in names(DeepSEM_list)) {
  colnames(DeepSEM_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(DeepSEM_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = DeepSEM_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 显示图
p_line

#DeepRIG
# 1. 加载DeepSEM算法预测的GRN
load("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Bun.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Guard.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Mesophyll.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Pavement.Rdata")
load("D:/玉米多模态数据的基因调控网络研究/DeepRig/DeepRig_Subsidiary.Rdata")
# ------------------------------------------------------------------------------
library(readxl)
#读取chip-seq数据
chip_seq <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 真实网络（chip-seq）
colnames(chip_seq)[1:2] <- c("TF", "Target")

# 五个细胞类型GRN列表
DeepRig_list <- list(
  Bun = DeepRig_Bun,
  Guard = DeepRig_Guard,
  Mesophyll = DeepRig_Mesophyll,
  Pavement = DeepRig_Pavement,
  Subsidiary = DeepRig_Subsidiary
)

# 统一列名
for (name in names(DeepRig_list)) {
  colnames(DeepRig_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(DeepRig_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = DeepRig_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 显示图
p_line



#LEAF
# 五个细胞类型GRN列表（改名为 GRNBoost2_list）
leaf_list <- list(
  GENIE3 = GENIE3_leaf,
  GRNBoost2 = GRNBoost2_leaf,
  KBoost = KBoost_leaf,
  DeepSEM = DeepSEM_leaf,
  DeepRig = DeepRig_leaf
)

# 统一列名
for (name in names(leaf_list)) {
  colnames(leaf_list[[name]])[1:2] <- c("TF", "Target")
}

# ------------------------------------------------------------------------------
# 2. overlap计算函数
# ------------------------------------------------------------------------------

calc_overlap <- function(pred_df, chip_df, top_n = Inf) {
  pred_edges <- paste(pred_df$TF, pred_df$Target)
  chip_edges <- paste(chip_df$TF, chip_df$Target)
  sum(pred_edges %in% chip_edges)
}

# ------------------------------------------------------------------------------
# 3. 零模型函数
# ------------------------------------------------------------------------------

null_model_overlap <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  true_overlap <- calc_overlap(pred_df, chip_df, top_n)
  null_dist <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    shuffled <- pred_df %>%
      mutate(TF = sample(TF), Target = sample(Target)) %>%
      distinct(TF, Target, .keep_all = TRUE)
    
    null_dist[i] <- calc_overlap(shuffled, chip_df, top_n)
  }
  
  p_val = (sum(null_dist >= true_overlap) + 1) / (n_perm + 1)
  
  return(list(
    true_overlap = true_overlap,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    p_value = p_val,
    null = null_dist
  ))
}

# ------------------------------------------------------------------------------
# 4. 包装函数
# ------------------------------------------------------------------------------

evaluate_with_null <- function(pred_df, chip_df, top_n = Inf, n_perm = 1000) {
  
  null_res <- null_model_overlap(pred_df, chip_df, top_n, n_perm)
  
  return(list(
    true_overlap = null_res$true_overlap,
    null_mean = null_res$null_mean,
    null_sd = null_res$null_sd,
    p_value = null_res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 5. 批量运行（五个细胞类型）
# ------------------------------------------------------------------------------

results <- data.frame()

for (name in names(leaf_list)) {
  cat("正在评估：", name, "\n")
  
  res <- evaluate_with_null(
    pred_df = leaf_list[[name]],
    chip_df = chip_seq,
    top_n = Inf,
    n_perm = 1000
  )
  
  results <- rbind(results, data.frame(
    CellType = name,
    Overlap_real = res$true_overlap,
    Random_mean = res$null_mean,
    Random_sd = res$null_sd,
    P_value = res$p_value
  ))
}

# ------------------------------------------------------------------------------
# 6. 显著性标注
# ------------------------------------------------------------------------------

results$significance <- ifelse(results$P_value < 0.001, "***",
                               ifelse(results$P_value < 0.01, "**",
                                      ifelse(results$P_value < 0.05, "*", "ns")))

# 查看结果
print(results)
library(ggplot2)

# 设置细胞类型顺序（避免乱序）
results$CellType <- factor(results$CellType,
                           levels = c("Bun","Guard","Mesophyll","Pavement","Subsidiary"))

# 可视化
library(ggplot2)

# 构建和你原来一致的数据格式
final_result <- results
final_result$x <- 1:nrow(final_result)
final_result$Label <- final_result$CellType

p_line <- ggplot(final_result, aes(x = x)) +
  
  # 红线：真实网络
  geom_line(aes(y = Overlap_real, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Overlap_real, color = "Predicted Network"), size = 2.5) +
  
  # 蓝线：随机网络
  geom_line(aes(y = Random_mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_mean, color = "Random Mean"), size = 2.5) +
  
  # 显著性标注
  geom_text(aes(y = Overlap_real + 10, label = significance),
            size = 5, color = "black") +
  
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  
  scale_color_manual(
    values = c(
      "Predicted Network" = "#ee2a25",   # ⚠️ 必须加 #
      "Random Mean" = "#0076aa"
    )
  ) +
  
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Cell Type",
    y = "Overlap count",
    color = ""
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

