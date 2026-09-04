#计算5种算法的AUROC值
#scRNA-seq数据
# 载入必要的包
setwd("D:/玉米多模态数据的基因调控网络研究/AUROC")
library(pROC)  # 用于计算AUROC
library(PRROC) # 用于计算AUPR
library(readxl) # 用于读取Excel文件
library(tidyverse)
#GRNBoost算法
#读入GRNBoost2预测的GRN，将GRNBoost2方法基因调控网络全部合并到一起
predicted_network_GRNBoost2 <- bind_rows(GRNBoost2_Guard, GRNBoost2_Bun, GRNBoost2_pavement, GRNBoost2_Mesophyll, GRNBoost2_Subsidiary)
predicted_network_GRNBoost2 <- predicted_network_GRNBoost2[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_GRNBoost2) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_GRNBoost2$is_true <- with(predicted_network_GRNBoost2, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
GRNBoost2_scores <- as.numeric(predicted_network_GRNBoost2$Score)
GRNBoost2_labels <- as.numeric(predicted_network_GRNBoost2$is_true)
# 计算 AUROC
GRNBoost2_roc_obj <- roc(GRNBoost2_labels, GRNBoost2_scores)
GRNBoost2_auroc <- auc(GRNBoost2_roc_obj)
#计算 AUPR
GRNBoost2_pr_obj <- pr.curve(scores.class0 = GRNBoost2_scores, weights.class0 = GRNBoost2_labels)
GRNBoost2_aupr <- GRNBoost2_pr_obj$auc.davis.goadrich
# 打印结果
cat("AUROC:", GRNBoost2_auroc, "\n")
cat("AUPR:", GRNBoost2_aupr, "\n")
#AUROC: 0.5275175 

#GENIE3
# 读取 GENIE3-RF 算法的五个网络
GENIE3_RF_network1 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_Bun.xlsx")  # 替换为实际文件路径
GENIE3_RF_network2 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_Guard.xlsx")
GENIE3_RF_network3 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_Mesophyll.xlsx")
GENIE3_RF_network4 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_Pavement.xlsx")
GENIE3_RF_network5 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_Subsidiary.xlsx")
# 合并五个网络
predicted_network_GENIE3 <- bind_rows(GENIE3_RF_network1, GENIE3_RF_network2, GENIE3_RF_network3, GENIE3_RF_network4, GENIE3_RF_network5)
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/CHIP-seq数据.xlsx")

# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_GENIE3$is_true <- with(predicted_network_GENIE3, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
GENIE3_scores <- as.numeric(predicted_network_GENIE3$Score)
GENIE3_labels <- as.numeric(predicted_network_GENIE3$is_true)
# 计算 AUROC
GENIE3_roc_obj <- roc(GENIE3_labels, GENIE3_scores)
GENIE3_auroc <- auc(GENIE3_roc_obj)
# 打印结果
cat("AUROC:", GENIE3_auroc, "\n")
#AUROC: 0.5650393

#KBoost
# 读取 GENIE3-RF 算法的五个网络
KBoost_network1 <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_Bun.xlsx")  # 替换为实际文件路径
KBoost_network2 <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_Guard.xlsx")
KBoost_network3 <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_Mesophyll.xlsx")
KBoost_network4 <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_Pavement.xlsx")
KBoost_network5 <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_Subsidiary.xlsx")
# 合并五个网络
predicted_network_KBoost <- bind_rows(KBoost_network1, KBoost_network2, KBoost_network3, KBoost_network4, KBoost_network5)
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/CHIP-seq数据.xlsx")
colnames(predicted_network_KBoost) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_KBoost$is_true <- with(predicted_network_KBoost, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
KBoost_scores <- as.numeric(predicted_network_KBoost$Score)
KBoost_labels <- as.numeric(predicted_network_KBoost$is_true)
# 计算 AUROC
KBoost_roc_obj <- roc(KBoost_labels, KBoost_scores)
KBoost_auroc <- auc(KBoost_roc_obj)
# 打印结果
cat("AUROC:", KBoost_auroc, "\n")
#AUROC: 0.6339288

#DeepSEM算法
#读入DeepSEM预测的GRN，将DeepSEM方法基因调控网络全部合并到一起
predicted_network_DeepSEM <- bind_rows(DeepSEM_Bun, DeepSEM_Guard, DeepSEM_Mesophyll, DeepSEM_Pavement, DeepSEM_Subsidiary)
predicted_network_DeepSEM <- predicted_network_DeepSEM[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_DeepSEM) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_DeepSEM$is_true <- with(predicted_network_DeepSEM, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
DeepSEM_scores <- as.numeric(predicted_network_DeepSEM$Score)
DeepSEM_labels <- as.numeric(predicted_network_DeepSEM$is_true)
# 计算 AUROC
DeepSEM_roc_obj <- roc(DeepSEM_labels, DeepSEM_scores)
DeepSEM_auroc <- auc(DeepSEM_roc_obj)
# 打印结果
cat("AUROC:", DeepSEM_auroc, "\n")
cat("AUPR:", DeepSEM_aupr, "\n")
#AUROC: 0.5400587

#DeepRig
#读入DeepSEM预测的GRN，将DeepSEM方法基因调控网络全部合并到一起
predicted_network_DeepRig <- bind_rows(DeepRig_Bun, DeepRig_Guard, DeepRig_Pavement, DeepRig_Subsidiary, DeepRig_Mesophyll)
predicted_network_DeepRig <- predicted_network_DeepRig[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_DeepRig) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_DeepRig$is_true <- with(predicted_network_DeepRig, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
DeepRig_scores <- as.numeric(predicted_network_DeepRig$Score)
DeepRig_labels <- as.numeric(predicted_network_DeepRig$is_true)
# 计算 AUROC
DeepRig_roc_obj <- roc(DeepRig_labels, DeepRig_scores)
DeepRig_auroc <- auc(DeepRig_roc_obj)
# 打印结果
cat("AUROC:", DeepRig_auroc, "\n")
cat("AUPR:", DeepRig_aupr, "\n")
#AUROC: 0.5599978 

#可视化
# Step 1: 确保分数列是数值型
# -------------------------------
predicted_network_GRNBoost2$Score <- as.numeric(predicted_network_GRNBoost2$Score)
predicted_network_GENIE3$Score <- as.numeric(predicted_network_GENIE3$Score)
predicted_network_KBoost$Score <- as.numeric(predicted_network_KBoost$Score)
predicted_network_DeepSEM$Score <- as.numeric(predicted_network_DeepSEM$Score)
predicted_network_DeepRig$Score <- as.numeric(predicted_network_DeepRig$Score)

# -------------------------------
# Step 2: 计算 ROC 对象
# -------------------------------
roc_GRNBoost2 <- roc(predicted_network_GRNBoost2$is_true, predicted_network_GRNBoost2$Score)
roc_GENIE3    <- roc(predicted_network_GENIE3$is_true, predicted_network_GENIE3$Score)
roc_KBoost    <- roc(predicted_network_KBoost$is_true, predicted_network_KBoost$Score)
roc_KBoost <- smooth(roc_KBoost)  # 平滑曲线
roc_DeepSEM   <- roc(predicted_network_DeepSEM$is_true, predicted_network_DeepSEM$Score)
roc_DeepRig   <- roc(predicted_network_DeepRig$is_true, predicted_network_DeepRig$Score)

# -------------------------------
# Step 3: 保存到列表，方便迭代绘图/比较
# -------------------------------
roc_list <- list(
  GRNBoost2  = roc_GRNBoost2,
  GENIE3 = roc_GENIE3,
  KBoost       = roc_KBoost,
  DeepSEM      = roc_DeepSEM,
  DeepRig      = roc_DeepRig
)
# 计算 AUCs
aucs <- sapply(roc_list, auc)
# 设置颜色
cols <- c("GRNBoost2" = "red", "GENIE3" = "gold", "KBoost" = "blue","DeepSEM" = "green","DeepRig" = "pink")
#绘图
plot(
  NA, NA,
  xlim = c(0, 1), ylim = c(0, 1),
  xlab = "1 - Specificity", ylab = "Sensitivity",
  main = "ROC Curve Comparison of GRN Inference Methods"
)
# 添加 ROC 曲线
for (method in names(roc_list)) {
  rocobj <- roc_list[[method]]
  x <- 1 - rocobj$specificities
  y <- rocobj$sensitivities
  lines(x, y, col = cols[method], lwd = 2)
}
# 添加对角线
abline(0, 1, lty = 2, col = "grey")
# 添加图例
legend(
  "bottomright",
  legend = paste0(
    names(aucs),
    " (AUC=", sprintf("%.3f", aucs), ")"
  ),
  col = cols[names(aucs)],
  lwd = 2, bty = "n"
)
library(ggplot2)

library(pROC)

# ---- 确保所有的roc对象都已计算好 ----
# 假设以下对象都已在你的环境中存在：
# GRNBoost2_roc_obj, GENIE3_roc_obj, KBoost_roc_obj, DeepSEM_roc_obj, DeepRig_roc_obj

# ---- 绘制ROC曲线 ----
plot(GRNBoost2_roc_obj, col = "#1b9e77", lwd = 2, legacy.axes = TRUE,
     main = "Comparison of ROC Curves for Five GRN Inference Methods",
     cex.main = 1.2, cex.lab = 1.1, cex.axis = 1,
     print.auc = FALSE, 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

# 添加其他曲线
plot(GENIE3_roc_obj, col = "#d95f02", lwd = 2, add = TRUE)
plot(KBoost_roc_obj, col = "#7570b3", lwd = 2, add = TRUE)
plot(DeepSEM_roc_obj, col = "#e7298a", lwd = 2, add = TRUE)
plot(DeepRig_roc_obj, col = "#66a61e", lwd = 2, add = TRUE)

# ---- 添加对角线 ----
abline(a = 0, b = 1, lty = 2, col = "gray")

# ---- 添加图例 ----
legend("bottomright", 
       legend = c(
         paste0("GRNBoost2 (AUROC=", round(auc(GRNBoost2_roc_obj), 3), ")"),
         paste0("GENIE3 (AUROC=", round(auc(GENIE3_roc_obj), 3), ")"),
         paste0("KBoost (AUROC=", round(auc(KBoost_roc_obj), 3), ")"),
         paste0("DeepSEM (AUROC=", round(auc(DeepSEM_roc_obj), 3), ")"),
         paste0("DeepRig (AUROC=", round(auc(DeepRig_roc_obj), 3), ")")
       ),
       col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
       lwd = 2, cex = 0.9, bty = "n")



#bulk数据
#读取GRNBoost2的基因调控网络
#读入GRNBoost2预测的GRN，将GRNBoost2方法基因调控网络全部合并到一起
predicted_network_GRNBoost2 <- GRNBoost2_leaf[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_GRNBoost2) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_GRNBoost2$is_true <- with(predicted_network_GRNBoost2, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
GRNBoost2_scores <- as.numeric(predicted_network_GRNBoost2$Score)
GRNBoost2_labels <- as.numeric(predicted_network_GRNBoost2$is_true)
# 计算 AUROC
GRNBoost2_roc_obj <- roc(GRNBoost2_labels, GRNBoost2_scores)
GRNBoost2_auroc <- auc(GRNBoost2_roc_obj)
#计算 AUPR
GRNBoost2_pr_obj <- pr.curve(scores.class0 = GRNBoost2_scores, weights.class0 = GRNBoost2_labels)
GRNBoost2_aupr <- GRNBoost2_pr_obj$auc.davis.goadrich
# 打印结果
cat("AUROC:", GRNBoost2_auroc, "\n")
cat("AUPR:", GRNBoost2_aupr, "\n")
#AUROC: 0.5478547 

#GENIE3
# 读取 GENIE3算法的网络
predicted_network_GENIE3 <- read_excel("D:/玉米多模态数据的基因调控网络研究/GENIE3推断GRN/RF调控关系权重/edges_with_weights_leaf.xlsx")
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/CHIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_GENIE3) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_GENIE3$is_true <- with(predicted_network_GENIE3, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
GENIE3_scores <- as.numeric(predicted_network_GENIE3$Score)
GENIE3_labels <- as.numeric(predicted_network_GENIE3$is_true)
# 计算 AUROC
GENIE3_roc_obj <- roc(GENIE3_labels, GENIE3_scores)
GENIE3_auroc <- auc(GENIE3_roc_obj)
# 打印结果
cat("AUROC:", GENIE3_auroc, "\n")
#AUROC: 0.5486655

#KBoost
# 读取 GENIE3-RF 算法的五个网络
predicted_network_KBoost <- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost/edges_with_weights_leaf.xlsx")
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/CHIP-seq数据.xlsx")
colnames(predicted_network_KBoost) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_KBoost$is_true <- with(predicted_network_KBoost, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
KBoost_scores <- as.numeric(predicted_network_KBoost$Score)
KBoost_labels <- as.numeric(predicted_network_KBoost$is_true)
# 计算 AUROC
KBoost_roc_obj <- roc(KBoost_labels, KBoost_scores)
KBoost_auroc <- auc(KBoost_roc_obj)
# 打印结果
cat("AUROC:", KBoost_auroc, "\n")
#AUROC: 0.5088386

#DeepSEM算法
#读入DeepSEM预测的GRN，将DeepSEM方法基因调控网络全部合并到一起
predicted_network_DeepSEM <- DeepSEM_leaf
predicted_network_DeepSEM <- predicted_network_DeepSEM[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_DeepSEM) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_DeepSEM$is_true <- with(predicted_network_DeepSEM, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
DeepSEM_scores <- as.numeric(predicted_network_DeepSEM$Score)
DeepSEM_labels <- as.numeric(predicted_network_DeepSEM$is_true)
# 计算 AUROC
DeepSEM_roc_obj <- roc(DeepSEM_labels, DeepSEM_scores)
DeepSEM_auroc <- auc(DeepSEM_roc_obj)
# 打印结果
cat("AUROC:", DeepSEM_auroc, "\n")
cat("AUPR:", DeepSEM_aupr, "\n")
#AUROC: 0.497821

#DeepRig
#读入DeepSEM预测的GRN，将DeepSEM方法基因调控网络全部合并到一起
predicted_network_DeepRig <- DeepRig_leaf
predicted_network_DeepRig <- predicted_network_DeepRig[, 1:3]
# 读取真实的ChIP-seq验证网络
true_network <- read_excel("D:/玉米多模态数据的基因调控网络研究/ChIP-seq数据.xlsx")
# 确保数据框列名一致
colnames(predicted_network_DeepRig) <- c("TF", "Target", "Score")
colnames(true_network) <- c("TF", "Target")
# 生成一个标记向量，用来标记每个预测关系是否出现在真实网络中
predicted_network_DeepRig$is_true <- with(predicted_network_DeepRig, paste(TF, Target)) %in%
  with(true_network, paste(TF, Target))
# 提取预测得分和真实标签
DeepRig_scores <- as.numeric(predicted_network_DeepRig$Score)
DeepRig_labels <- as.numeric(predicted_network_DeepRig$is_true)
# 计算 AUROC
DeepRig_roc_obj <- roc(DeepRig_labels, DeepRig_scores)
DeepRig_auroc <- auc(DeepRig_roc_obj)
# 打印结果
cat("AUROC:", DeepRig_auroc, "\n")
cat("AUPR:", DeepRig_aupr, "\n")
#AUROC: 0.5037876 
# ---- 绘制ROC曲线 ----
predicted_network_GRNBoost2$Score <- as.numeric(predicted_network_GRNBoost2$Score)
predicted_network_GENIE3$Score <- as.numeric(predicted_network_GENIE3$Score)
predicted_network_KBoost$Score <- as.numeric(predicted_network_KBoost$Score)
predicted_network_DeepSEM$Score <- as.numeric(predicted_network_DeepSEM$Score)
predicted_network_DeepRig$Score <- as.numeric(predicted_network_DeepRig$Score)
roc_GRNBoost2 <- roc(predicted_network_GRNBoost2$is_true, predicted_network_GRNBoost2$Score)
roc_GENIE3 <- roc(predicted_network_GENIE3$is_true, predicted_network_GENIE3$Score)
roc_KBoost <- roc(predicted_network_KBoost$is_true, predicted_network_KBoost$Score)
roc_DeepSEM <- roc(predicted_network_DeepSEM$is_true, predicted_network_DeepSEM$Score)
roc_DeepRig <- roc(predicted_network_DeepRig$is_true, predicted_network_DeepRig$Score)
roc_df <- do.call(rbind, lapply(names(roc_list), function(method) {
  rocobj <- roc_list[[method]]
  data.frame(
    Method = method,
    FPR = 1 - rocobj$specificities,
    TPR = rocobj$sensitivities
  )
}))
library(ggplot2)

# 自定义颜色
cols <- c(
  "GRNBoost2" = "#e41a1c",
  "GENIE3"    = "#ff7f00",
  "KBoost"    = "#377eb8",
  "DeepSEM"   = "#4daf4a",
  "DeepRig"   = "#984ea3"
)

ggplot(roc_df, aes(x = FPR, y = TPR, color = Method)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(
    title = "ROC Curve Comparison of GRN Inference Methods",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Method"
  ) +
  annotate(
    "text",
    x = 0.65, y = seq(0.4, 0.15, length.out = length(aucs)),
    label = paste0(names(aucs), " (AUC = ", sprintf("%.3f", aucs), ")"),
    color = cols[names(aucs)],
    hjust = 0
  ) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )





plot(
  NA, NA,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "ROC Curve Comparison of GRN Inference Methods"
)

# 添加 ROC 曲线
for (method in names(roc_list)) {
  rocobj <- roc_list[[method]]
  x <- 1 - rocobj$specificities
  y <- rocobj$sensitivities
  lines(x, y, col = cols[method], lwd = 2)
}

# 添加对角线（随机预测）
abline(0, 1, lty = 2, col = "grey")

# 添加图例
legend(
  "bottomright",
  legend = paste0(
    names(aucs),
    " (AUC = ",
    sprintf("%.3f", aucs),
    ")"
  ),
  col = cols[names(aucs)],
  lwd = 2,
  bty = "n",
  cex = 0.9
)
