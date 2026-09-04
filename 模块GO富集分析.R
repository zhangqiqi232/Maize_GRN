#用柱状图展示不同方法识别出来的不同细胞类型的模块
library(tidyverse)
# 原始数据
module_df <- data.frame(
  CellType = c("Bundle sheath", "Guard", "Mesophyll", "Pavement", "Subsidiary"),
  infomap = c(79, 98, 109, 107, 105),
  louvain = c(40, 82, 97, 88, 93),
  walktrap = c(54, 79, 84, 79, 78),
  eign = c(54, 72, 103, 93, 93),
  prop = c(84, 98, 108, 105, 105),
  FN = c(38, 83, 100, 88, 97),
  MCL = c(81, 98, 109, 107, 105),
  MCODE = c(4, 9, 6, 6, 5)
)
# 转为长表（ggplot 推荐格式）
module_long <- module_df %>%
  pivot_longer(
    cols = -CellType,
    names_to = "Method",
    values_to = "ModuleNumber"
  )
ggplot(module_long, aes(x = Method, y = ModuleNumber, fill = CellType)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  # ⭐ 在柱子上方加数值
  geom_text(
    aes(label = ModuleNumber),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(values = c(
    "#b5aed5",  # Bundle sheath
    "#b2e6fd",  # Guard
    "#b8d2cc",  # Mesophyll
    "#e8b2a7",  # Pavement
    "#feebb9"   # Subsidiary
  )) +
  labs(
    x = "Community detection method",
    y = "Number of modules",
    title = "Comparison of module numbers across methods and cell types"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  # 给 y 轴留一点空间，避免数字被截掉
  expand_limits(y = max(module_long$ModuleNumber) * 1.1)


#识别MCL方法中含有关键TFs最多的前五个模块
#Bundle sheath
# 1. 统计每个模块中关键 TF 的数量
# 统计每个模块的 TF 数量和模块大小
module_TF_count_Bun <- data.frame(
  Module = seq_along(Bun_Cluster_MCL),
  TF_number = sapply(
    Bun_Cluster_MCL,
    function(genes) sum(genes %in% Bun_key_TFs)
  ),
  Module_size = sapply(Bun_Cluster_MCL, length)
)

# ⭐ 核心：排序 + 取前 5
top5_modules_Bun <- module_TF_count_Bun %>%
  arrange(desc(TF_number), desc(Module_size)) %>%  # 先 TF，再模块大小
  slice(1:5)
Bun_MCL_modules <- Bun_Cluster_MCL[top5_modules_Bun$Module]
library(openxlsx)
# 提取第 1 个模块的基因
Bun_module1<- Bun_MCL_modules[[1]]
# 转为数据框
Bun_module1 <- data.frame(Gene = Bun_module1)
Bun_module1 <- Bun_module1[grepl("^Zm", Bun_module1[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun")
# 保存为 Excel
write.xlsx(
  Bun_module1,
  file = "Bun_module1.xlsx",
  rowNames = FALSE
)

# 提取第2 个模块的基因
Bun_module2<- Bun_MCL_modules[[2]]
# 转为数据框
Bun_module2 <- data.frame(Gene = Bun_module2)
Bun_module2 <- Bun_module1[grepl("^Zm", Bun_module1[[2]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun")
# 保存为 Excel
write.xlsx(
  Bun_module2,
  file = "Bun_module2.xlsx",
  rowNames = FALSE
)

# 提取第 3 个模块的基因
Bun_module3<- Bun_MCL_modules[[3]]
# 转为数据框
Bun_module3 <- data.frame(Gene = Bun_module3)
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun")
# 保存为 Excel
write.xlsx(
  Bun_module3,
  file = "Bun_module3.xlsx",
  rowNames = FALSE
)
# 提取第 4个模块的基因
Bun_module4 <- Bun_MCL_modules[[4]]
# 转为数据框
Bun_module4 <- data.frame(Gene = Bun_module4)
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun")
# 保存为 Excel
write.xlsx(
  Bun_module4,
  file = "Bun_module4.xlsx",
  rowNames = FALSE
)
# 提取第 5个模块的基因
Bun_module5 <- Bun_MCL_modules[[5]]
# 转为数据框
Bun_module5 <- data.frame(Gene = Bun_module5)
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun")
# 保存为 Excel
write.xlsx(
  Bun_module5,
  file = "Bun_module5.xlsx",
  rowNames = FALSE
)

#Guard
# 1. 统计每个模块中关键 TF 的数量
# 统计每个模块的 TF 数量和模块大小
module_TF_count_Guard <- data.frame(
  Module = seq_along(Guard_Cluster_MCL),
  TF_number = sapply(
    Guard_Cluster_MCL,
    function(genes) sum(genes %in% Guard_key_TFs)
  ),
  Module_size = sapply(Guard_Cluster_MCL, length)
)

# ⭐ 核心：排序 + 取前 5
top5_modules_Guard <- module_TF_count_Guard %>%
  arrange(desc(TF_number), desc(Module_size)) %>%  # 先 TF，再模块大小
  slice(1:5)
Guard_MCL_modules <- Guard_Cluster_MCL[top5_modules_Guard$Module]
library(openxlsx)
# 提取第 1 个模块的基因
Guard_module1<- Guard_MCL_modules[[1]]
# 转为数据框
Guard_module1 <- data.frame(Gene = Guard_module1)
Guard_module1 <- Guard_module1[grepl("^Zm", Guard_module1[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard")
# 保存为 Excel
write.xlsx(
  Guard_module1,
  file = "Guard_module1.xlsx",
  rowNames = FALSE
)
# 提取第 2个模块的基因
Guard_module2<- Guard_MCL_modules[[2]]
# 转为数据框
Guard_module2 <- data.frame(Gene = Guard_module2)
Guard_module2 <- Guard_module2[grepl("^Zm", Guard_module2[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard")
# 保存为 Excel
write.xlsx(
  Guard_module2,
  file = "Guard_module2.xlsx",
  rowNames = FALSE
)
# 提取第3个模块的基因
Guard_module3<- Guard_MCL_modules[[3]]
# 转为数据框
Guard_module3 <- data.frame(Gene = Guard_module3)
Guard_module3 <- Guard_module3[grepl("^Zm", Guard_module3[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard")
# 保存为 Excel
write.xlsx(
  Guard_module3,
  file = "Guard_module3.xlsx",
  rowNames = FALSE
)
# 提取第4个模块的基因
Guard_module4<- Guard_MCL_modules[[4]]
# 转为数据框
Guard_module4 <- data.frame(Gene = Guard_module4)
Guard_module4 <- Guard_module4[grepl("^Zm", Guard_module4[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard")
# 保存为 Excel
write.xlsx(
  Guard_module4,
  file = "Guard_module4.xlsx",
  rowNames = FALSE
)
# 提取第5个模块的基因
Guard_module5<- Guard_MCL_modules[[5]]
# 转为数据框
Guard_module5 <- data.frame(Gene = Guard_module5)
Guard_module5 <- Guard_module5[grepl("^Zm", Guard_module5[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard")
# 保存为 Excel
write.xlsx(
  Guard_module5,
  file = "Guard_module5.xlsx",
  rowNames = FALSE
)

#Mesophyll
# 1. 统计每个模块中关键 TF 的数量
# 统计每个模块的 TF 数量和模块大小
module_TF_count_Mesophyll <- data.frame(
  Module = seq_along(Mesophyll_Cluster_MCL),
  TF_number = sapply(
    Mesophyll_Cluster_MCL,
    function(genes) sum(genes %in% Mesophyll_key_TFs)
  ),
  Module_size = sapply(Mesophyll_Cluster_MCL, length)
)

# ⭐ 核心：排序 + 取前 5
top5_modules_Mesophyll <- module_TF_count_Mesophyll %>%
  arrange(desc(TF_number), desc(Module_size)) %>%  # 先 TF，再模块大小
  slice(1:5)
Mesophyll_MCL_modules <- Mesophyll_Cluster_MCL[top5_modules_Mesophyll$Module]
library(openxlsx)
# 提取第 1 个模块的基因
Mesophyll_module1<- Mesophyll_MCL_modules[[1]]
# 转为数据框
Mesophyll_module1 <- data.frame(Gene = Mesophyll_module1)
Mesophyll_module1 <- Mesophyll_module1[grepl("^Zm", Mesophyll_module1[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll")
# 保存为 Excel
write.xlsx(
  Mesophyll_module1,
  file = "Mesophyll_module1.xlsx",
  rowNames = FALSE
)
# 提取第 2 个模块的基因
Mesophyll_module2<- Mesophyll_MCL_modules[[2]]
# 转为数据框
Mesophyll_module2 <- data.frame(Gene = Mesophyll_module2)
Mesophyll_module2 <- Mesophyll_module2[grepl("^Zm", Mesophyll_module2[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll")
# 保存为 Excel
write.xlsx(
  Mesophyll_module2,
  file = "Mesophyll_module2.xlsx",
  rowNames = FALSE
)
# 提取第 3 个模块的基因
Mesophyll_module3<- Mesophyll_MCL_modules[[3]]
# 转为数据框
Mesophyll_module3 <- data.frame(Gene = Mesophyll_module3)
Mesophyll_module3 <- Mesophyll_module3[grepl("^Zm", Mesophyll_module3[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll")
# 保存为 Excel
write.xlsx(
  Mesophyll_module3,
  file = "Mesophyll_module3.xlsx",
  rowNames = FALSE
)
# 提取第 4 个模块的基因
Mesophyll_module4<- Mesophyll_MCL_modules[[4]]
# 转为数据框
Mesophyll_module4 <- data.frame(Gene = Mesophyll_module4)
Mesophyll_module4 <- Mesophyll_module4[grepl("^Zm", Mesophyll_module4[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll")
# 保存为 Excel
write.xlsx(
  Mesophyll_module4,
  file = "Mesophyll_module4.xlsx",
  rowNames = FALSE
)
# 提取第 5 个模块的基因
Mesophyll_module5<- Mesophyll_MCL_modules[[5]]
# 转为数据框
Mesophyll_module5 <- data.frame(Gene = Mesophyll_module5)
Mesophyll_module5 <- Mesophyll_module5[grepl("^Zm", Mesophyll_module5[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll")
# 保存为 Excel
write.xlsx(
  Mesophyll_module5,
  file = "Mesophyll_module5.xlsx",
  rowNames = FALSE
)

#Pavement
# 1. 统计每个模块中关键 TF 的数量
# 统计每个模块的 TF 数量和模块大小
module_TF_count_Pavement <- data.frame(
  Module = seq_along(Pavement_Cluster_MCL),
  TF_number = sapply(
    Pavement_Cluster_MCL,
    function(genes) sum(genes %in% Pavement_key_TFs)
  ),
  Module_size = sapply(Pavement_Cluster_MCL, length)
)

# ⭐ 核心：排序 + 取前 5
top5_modules_Pavement<- module_TF_count_Pavement %>%
  arrange(desc(TF_number), desc(Module_size)) %>%  # 先 TF，再模块大小
  slice(1:5)
Pavement_MCL_modules <- Pavement_Cluster_MCL[top5_modules_Pavement$Module]
library(openxlsx)
# 提取第 1 个模块的基因
Pavement_module1<- Pavement_MCL_modules[[1]]
# 转为数据框
Pavement_module1 <- data.frame(Gene = Pavement_module1)
Pavement_module1 <- Pavement_module1[grepl("^Zm", Pavement_module1[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement")
# 保存为 Excel
write.xlsx(
  Pavement_module1,
  file = "Pavement_module1.xlsx",
  rowNames = FALSE
)
# 提取第 2 个模块的基因
Pavement_module2<- Pavement_MCL_modules[[2]]
# 转为数据框
Pavement_module2 <- data.frame(Gene = Pavement_module2)
Pavement_module2 <- Pavement_module2[grepl("^Zm", Pavement_module2[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement")
# 保存为 Excel
write.xlsx(
  Pavement_module2,
  file = "Pavement_module2.xlsx",
  rowNames = FALSE
)
# 提取第 3 个模块的基因
Pavement_module3<- Pavement_MCL_modules[[3]]
# 转为数据框
Pavement_module3 <- data.frame(Gene = Pavement_module3)
Pavement_module3 <- Pavement_module3[grepl("^Zm", Pavement_module3[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement")
# 保存为 Excel
write.xlsx(
  Pavement_module3,
  file = "Pavement_module3.xlsx",
  rowNames = FALSE
)
# 提取第 4 个模块的基因
Pavement_module4<- Pavement_MCL_modules[[4]]
# 转为数据框
Pavement_module4 <- data.frame(Gene = Pavement_module4)
Pavement_module4 <- Pavement_module4[grepl("^Zm", Pavement_module4[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement")
# 保存为 Excel
write.xlsx(
  Pavement_module4,
  file = "Pavement_module4.xlsx",
  rowNames = FALSE
)
# 提取第 5 个模块的基因
Pavement_module5<- Pavement_MCL_modules[[5]]
# 转为数据框
Pavement_module5 <- data.frame(Gene = Pavement_module5)
Pavement_module5 <- Pavement_module5[grepl("^Zm", Pavement_module5[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement")
# 保存为 Excel
write.xlsx(
  Pavement_module5,
  file = "Pavement_module5.xlsx",
  rowNames = FALSE
)

#Subsidiary
# 1. 统计每个模块中关键 TF 的数量
# 统计每个模块的 TF 数量和模块大小
module_TF_count_Subsidiary <- data.frame(
  Module = seq_along(Subsidiary_Cluster_MCL),
  TF_number = sapply(
    Subsidiary_Cluster_MCL,
    function(genes) sum(genes %in% Subsidiary_key_TFs)
  ),
  Module_size = sapply(Subsidiary_Cluster_MCL, length)
)

# ⭐ 核心：排序 + 取前 5
top5_modules_Subsidiary<- module_TF_count_Subsidiary %>%
  arrange(desc(TF_number), desc(Module_size)) %>%  # 先 TF，再模块大小
  slice(1:5)
Subsidiary_MCL_modules <- Subsidiary_Cluster_MCL[top5_modules_Subsidiary$Module]
library(openxlsx)
# 提取第 1 个模块的基因
Subsidiary_module1<- Subsidiary_MCL_modules[[1]]
# 转为数据框
Subsidiary_module1 <- data.frame(Gene = Subsidiary_module1)
Subsidiary_module1 <- Subsidiary_module1[grepl("^Zm", Subsidiary_module1[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary")
# 保存为 Excel
write.xlsx(
  Subsidiary_module1,
  file = "Subsidiary_module1.xlsx",
  rowNames = FALSE
)
# 提取第 2 个模块的基因
Subsidiary_module2<- Subsidiary_MCL_modules[[2]]
# 转为数据框
Subsidiary_module2 <- data.frame(Gene = Subsidiary_module2)
Subsidiary_module2 <- Subsidiary_module2[grepl("^Zm", Subsidiary_module2[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary")
# 保存为 Excel
write.xlsx(
  Subsidiary_module2,
  file = "Subsidiary_module2.xlsx",
  rowNames = FALSE
)
# 提取第 3 个模块的基因
Subsidiary_module3<- Subsidiary_MCL_modules[[3]]
# 转为数据框
Subsidiary_module3 <- data.frame(Gene = Subsidiary_module3)
Subsidiary_module3 <- Subsidiary_module3[grepl("^Zm", Subsidiary_module3[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary")
# 保存为 Excel
write.xlsx(
  Subsidiary_module3,
  file = "Subsidiary_module3.xlsx",
  rowNames = FALSE
)
# 提取第4 个模块的基因
Subsidiary_module4<- Subsidiary_MCL_modules[[4]]
# 转为数据框
Subsidiary_module4 <- data.frame(Gene = Subsidiary_module4)
Subsidiary_module4 <- Subsidiary_module4[grepl("^Zm", Subsidiary_module4[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary")
# 保存为 Excel
write.xlsx(
  Subsidiary_module4,
  file = "Subsidiary_module4.xlsx",
  rowNames = FALSE
)
# 提取第 5 个模块的基因
Subsidiary_module5<- Subsidiary_MCL_modules[[5]]
# 转为数据框
Subsidiary_module5 <- data.frame(Gene = Subsidiary_module5)
Subsidiary_module5 <- Subsidiary_module5[grepl("^Zm", Subsidiary_module5[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary")
# 保存为 Excel
write.xlsx(
  Subsidiary_module5,
  file = "Subsidiary_module5.xlsx",
  rowNames = FALSE
)
#对每个细胞类型的模块进行GO富集分析
library("AnnotationHub")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
hub <- AnnotationHub()
query(hub, "zea mays")
zeamaize <- hub[['AH114308']]
#使用maizeGDB数据库进行基因id转换
#GO富集分析
#Bundle sheath模块
library(readxl)
library(clusterProfiler)
#模块1
Bun_module1<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun/转换后的id/Bun_module1.xlsx")
Bun_module1 <- Bun_module1$`Gene`
Bun_module1_GO <- enrichGO(gene = Bun_module1,#我们上面定义了
                    OrgDb=zeamaize,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Bun_module1_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块2
Bun_module2<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun/转换后的id/Bun_module2.xlsx")
Bun_module2 <- Bun_module2$`Gene`
Bun_module2_GO <- enrichGO(gene = Bun_module2,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Bun_module2_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块3
Bun_module3<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun/转换后的id/Bun_module3.xlsx")
Bun_module3 <- Bun_module3$`Gene`
Bun_module3_GO <- enrichGO(gene = Bun_module3,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
#模块4
Bun_module4<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun/转换后的id/Bun_module4_id.xlsx")
Bun_module4 <- Bun_module4$`Gene`
Bun_module4_GO <- enrichGO(gene = Bun_module4,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Bun_module4_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块5
Bun_module5<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Bun/转换后的id/Bun_module5_id.xlsx")
Bun_module5 <- Bun_module5$`Gene`
Bun_module5_GO <- enrichGO(gene = Bun_module5,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
#Guard模块
#模块1
Guard_module1<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard/id/Guard_module1_id.xlsx")
Guard_module1 <- Guard_module1$`Gene`
Guard_module1_GO <- enrichGO(gene = Guard_module1,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Guard_module1_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块2
Guard_module2<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard/id/Guard_module2_id.xlsx")
Guard_module2 <- Guard_module2$`Gene`
Guard_module2_GO <- enrichGO(gene = Guard_module2,#我们上面定义了
                             OrgDb=zeamaize,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Guard_module2_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO

#模块3
Guard_module3<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard/id/Guard_module3_id.xlsx")
Guard_module3 <- Guard_module3$`Gene`
Guard_module3_GO <- enrichGO(gene = Guard_module3,#我们上面定义了
                             OrgDb=zeamaize,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Guard_module3_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块4
Guard_module4<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard/id/Guard_module4_id.xlsx")
Guard_module4 <- Guard_module4$`Gene`
Guard_module4_GO <- enrichGO(gene = Guard_module4,#我们上面定义了
                             OrgDb=zeamaize,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
#模块5
Guard_module5<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Guard/id/Guard_module5_id.xlsx")
Guard_module5 <- Guard_module5$`Gene`
Guard_module5_GO <- enrichGO(gene = Guard_module5,#我们上面定义了
                             OrgDb=zeamaize,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Guard_module5_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#Mesophyll模块
#模块1
Mesophyll_module1<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll/id/Mesophyll_module1_id.xlsx")
Mesophyll_module1 <- Mesophyll_module1$`Gene`
Mesophyll_module1_GO <- enrichGO(gene = Mesophyll_module1,#我们上面定义了
                             OrgDb=zeamaize,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Mesophyll_module1_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块2
Mesophyll_module2<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll/id/Mesophyll_module2_id.xlsx")
Mesophyll_module2 <- Mesophyll_module2$`Gene`
Mesophyll_module2_GO <- enrichGO(gene = Mesophyll_module2,#我们上面定义了
                                 OrgDb=zeamaize,
                                 keyType = "ENTREZID",
                                 ont = "ALL",#富集的GO类型
                                 pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,#P值可以取0.05
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Mesophyll_module2_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块3
Mesophyll_module3<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll/id/Mesophyll_module3_id.xlsx")
Mesophyll_module3 <- Mesophyll_module3$`Gene`
Mesophyll_module3_GO <- enrichGO(gene = Mesophyll_module3,#我们上面定义了
                                 OrgDb=zeamaize,
                                 keyType = "ENTREZID",
                                 ont = "ALL",#富集的GO类型
                                 pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,#P值可以取0.05
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Mesophyll_module3_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块4
Mesophyll_module4<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll/id/Mesophyll_module4_id.xlsx")
Mesophyll_module4 <- Mesophyll_module4$`Gene`
Mesophyll_module4_GO <- enrichGO(gene = Mesophyll_module4,#我们上面定义了
                                 OrgDb=zeamaize,
                                 keyType = "ENTREZID",
                                 ont = "ALL",#富集的GO类型
                                 pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,#P值可以取0.05
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Mesophyll_module4_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块5
Mesophyll_module5<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Mesophyll/id/Mesophyll_module5_id.xlsx")
Mesophyll_module5 <- Mesophyll_module5$`Gene`
Mesophyll_module5_GO <- enrichGO(gene = Mesophyll_module5,#我们上面定义了
                                 OrgDb=zeamaize,
                                 keyType = "ENTREZID",
                                 ont = "ALL",#富集的GO类型
                                 pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,#P值可以取0.05
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
#Pavement模块
#模块1
Pavement_module1<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement/id/Pavement_module1_id.xlsx")
Pavement_module1 <- Pavement_module1$`Gene`
Pavement_module1_GO <- enrichGO(gene = Pavement_module1,#我们上面定义了
                                 OrgDb=zeamaize,
                                 keyType = "ENTREZID",
                                 ont = "ALL",#富集的GO类型
                                 pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,#P值可以取0.05
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
#模块2
Pavement_module2<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement/id/Pavement_module2_id.xlsx")
Pavement_module2 <- Pavement_module2$`Gene`
Pavement_module2_GO <- enrichGO(gene = Pavement_module2,#我们上面定义了
                                OrgDb=zeamaize,
                                keyType = "ENTREZID",
                                ont = "ALL",#富集的GO类型
                                pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                minGSSize = 1,
                                pvalueCutoff = 0.05,#P值可以取0.05
                                qvalueCutoff = 0.05,
                                readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Pavement_module2_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块3
Pavement_module3<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement/id/Pavement_module3_id.xlsx")
Pavement_module3 <- Pavement_module3$`Gene`
Pavement_module3_GO <- enrichGO(gene = Pavement_module3,#我们上面定义了
                                OrgDb=zeamaize,
                                keyType = "ENTREZID",
                                ont = "ALL",#富集的GO类型
                                pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                minGSSize = 1,
                                pvalueCutoff = 0.05,#P值可以取0.05
                                qvalueCutoff = 0.05,
                                readable = TRUE)
#模块4
Pavement_module4<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement/id/Pavement_module4_id.xlsx")
Pavement_module4 <- Pavement_module4$`Gene`
Pavement_module4_GO <- enrichGO(gene = Pavement_module4,#我们上面定义了
                                OrgDb=zeamaize,
                                keyType = "ENTREZID",
                                ont = "ALL",#富集的GO类型
                                pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                minGSSize = 1,
                                pvalueCutoff = 0.05,#P值可以取0.05
                                qvalueCutoff = 0.05,
                                readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Pavement_module4_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块5
Pavement_module5<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Pavement/id/Pavement_module5_id.xlsx")
Pavement_module5 <- Pavement_module5$`Gene`
Pavement_module5_GO <- enrichGO(gene = Pavement_module5,#我们上面定义了
                                OrgDb=zeamaize,
                                keyType = "ENTREZID",
                                ont = "ALL",#富集的GO类型
                                pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                minGSSize = 1,
                                pvalueCutoff = 0.05,#P值可以取0.05
                                qvalueCutoff = 0.05,
                                readable = TRUE)
#Subsidiary模块
#模块1
Subsidiary_module1<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary/id/Subsidiary_module1_id.xlsx")
Subsidiary_module1 <- Subsidiary_module1$`Gene`
Subsidiary_module1_GO <- enrichGO(gene = Subsidiary_module1,#我们上面定义了
                                OrgDb=zeamaize,
                                keyType = "ENTREZID",
                                ont = "ALL",#富集的GO类型
                                pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                minGSSize = 1,
                                pvalueCutoff = 0.05,#P值可以取0.05
                                qvalueCutoff = 0.05,
                                readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Subsidiary_module1_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块2
Subsidiary_module2<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary/id/Subsidairy_module2_id.xlsx")
Subsidiary_module2 <- Subsidiary_module2$`Gene`
Subsidiary_module2_GO <- enrichGO(gene = Subsidiary_module2,#我们上面定义了
                                  OrgDb=zeamaize,
                                  keyType = "ENTREZID",
                                  ont = "ALL",#富集的GO类型
                                  pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                  minGSSize = 1,
                                  pvalueCutoff = 0.05,#P值可以取0.05
                                  qvalueCutoff = 0.05,
                                  readable = TRUE)
#模块3
Subsidiary_module3<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary/id/Subsidiary_module3_id.xlsx")
Subsidiary_module3 <- Subsidiary_module3$`Gene`
Subsidiary_module3_GO <- enrichGO(gene = Subsidiary_module3,#我们上面定义了
                                  OrgDb=zeamaize,
                                  keyType = "ENTREZID",
                                  ont = "ALL",#富集的GO类型
                                  pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                  minGSSize = 1,
                                  pvalueCutoff = 0.05,#P值可以取0.05
                                  qvalueCutoff = 0.05,
                                  readable = TRUE)
#模块4
Subsidiary_module4<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary/id/Subsidiary_module4_id.xlsx")
Subsidiary_module4 <- Subsidiary_module4$`Gene`
Subsidiary_module4_GO <- enrichGO(gene = Subsidiary_module4,#我们上面定义了
                                  OrgDb=zeamaize,
                                  keyType = "ENTREZID",
                                  ont = "ALL",#富集的GO类型
                                  pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                  minGSSize = 1,
                                  pvalueCutoff = 0.05,#P值可以取0.05
                                  qvalueCutoff = 0.05,
                                  readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(Subsidiary_module4_GO)

# 按pvalue排序并取前20条
dt2 <- dt2[order(dt2$pvalue), ]
dt2 <- dt2[1:20, ]

# 计算GeneRatio
dt2$GeneRatio <- sapply(dt2$GeneRatio, function(x){
  eval(parse(text = x))
})

# 绘制GO富集气泡图
p_GO <- ggplot(data = dt2, aes(x = Description, y = GeneRatio, fill = -log10(pvalue))) +
  geom_point(aes(size = Count), shape = 21, colour = "black") +
  labs(
    x = NULL,
    y = "GeneRatio",
    title = "GO Enrichment",
    fill = bquote("-"~Log[10]~"(Pvalue)"),
    size = "Gene Count"
  ) +
  scale_fill_gradientn(colours = c("#f7ca64", "#46bac2","#7e62a3")) +
  scale_size_continuous(range = c(4, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, color = "darkred"),
    axis.title.x = element_text(size = 10, color = "darkgreen"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "purple"),
    legend.text = element_text(size = 10, color = "brown")
  ) +
  scale_x_discrete(labels = function(dat) stringr::str_wrap(dat, width = 50))
# 显示图
p_GO
#模块5
Subsidiary_module5<- read_excel("D:/玉米多模态数据的基因调控网络研究/KBoost下游分析/模块GO富集分析/Subsidiary/id/Subsidiary_module5_id.xlsx")
Subsidiary_module5 <- Subsidiary_module5$`Gene`
Subsidiary_module5_GO <- enrichGO(gene = Subsidiary_module5,#我们上面定义了
                                  OrgDb=zeamaize,
                                  keyType = "ENTREZID",
                                  ont = "ALL",#富集的GO类型
                                  pAdjustMethod = "BH",#这个不用管，一般都用的BH
                                  minGSSize = 1,
                                  pvalueCutoff = 0.05,#P值可以取0.05
                                  qvalueCutoff = 0.05,
                                  readable = TRUE)


library(dplyr)
library(tidyr)
library(pheatmap)
GO_all <- list(
  # Bundle sheath
  Bun_Module1 = Bun_module1_GO,
  Bun_Module2 = Bun_module2_GO,
  Bun_Module3 = Bun_module3_GO,
  Bun_Module4 = Bun_module4_GO,
  Bun_Module5 = Bun_module5_GO,
  
  # Guard
  Guard_Module1 = Guard_module1_GO,
  Guard_Module2 = Guard_module2_GO,
  Guard_Module3 = Guard_module3_GO,
  Guard_Module4 = Guard_module4_GO,
  Guard_Module5 = Guard_module5_GO,
  
  # Mesophyll
  Mesophyll_Module1 = Mesophyll_module1_GO,
  Mesophyll_Module2 = Mesophyll_module2_GO,
  Mesophyll_Module3 = Mesophyll_module3_GO,
  Mesophyll_Module4 = Mesophyll_module4_GO,
  Mesophyll_Module5 = Mesophyll_module5_GO,
  
  # Pavement
  Pavement_Module1 = Pavement_module1_GO,
  Pavement_Module2 = Pavement_module2_GO,
  Pavement_Module3 = Pavement_module3_GO,
  Pavement_Module4 = Pavement_module4_GO,
  Pavement_Module5 = Pavement_module5_GO,
  
  # Subsidiary
  Subsidiary_Module1 = Subsidiary_module1_GO,
  Subsidiary_Module2 = Subsidiary_module2_GO,
  Subsidiary_Module3 = Subsidiary_module3_GO,
  Subsidiary_Module4 = Subsidiary_module4_GO,
  Subsidiary_Module5 = Subsidiary_module5_GO
)
extract_GO <- function(go_obj, module_name, top_n = 20) {
  
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    return(NULL)
  }
  
  as.data.frame(go_obj) %>%
    filter(ONTOLOGY == "BP") %>%      # 只用 Biological Process
    arrange(qvalue) %>%
    slice_head(n = top_n) %>%
    mutate(
      Module = module_name,
      score = -log10(qvalue)
    ) %>%
    select(Description, Module, score)
}
GO_long <- bind_rows(
  lapply(names(GO_all), function(m) {
    extract_GO(GO_all[[m]], m, top_n = 5)
  })
)
GO_mat <- GO_long %>%
  pivot_wider(
    names_from = Module,
    values_from = score,
    values_fill = 0
  )

GO_mat <- as.data.frame(GO_mat)
rownames(GO_mat) <- GO_mat$Description
GO_mat$Description <- NULL
pheatmap(
  GO_mat,
  scale = "row",   # 按 GO term 标准化，突出模块差异
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "#FDB863", "#B2182B"))(100),
  fontsize_row = 7,
  fontsize_col = 9,
  main = "GO enrichment landscape across cell-type-specific regulatory modules"
)



# 从列名中提取细胞类型
cell_type <- sub("_.*", "", colnames(GO_mat))

# 构建列注释数据框
annotation_col <- data.frame(
  CellType = cell_type
)

rownames(annotation_col) <- colnames(GO_mat)
# 根据细胞类型排序列
order_cols <- order(annotation_col$CellType)

GO_mat_ordered <- GO_mat[, order_cols]
annotation_col <- annotation_col[order_cols, , drop = FALSE]
pheatmap(
  GO_mat_ordered,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # ❗非常关键：禁止打乱列顺序
  annotation_col = annotation_col,  # 加列注释
  annotation_colors = list(
    CellType = c(
      Bun = "#1b9e77",
      Guard = "#d95f02",
      Mesophyll = "#7570b3",
      Pavement = "#e7298a",
      Subsidiary = "#66a61e"
    )
  ),
  color = colorRampPalette(c("#e9c6c6", "#f9efef", "#98cadd"))(100),
  fontsize_row = 7,
  fontsize_col = 9,
  main = "GO enrichment landscape across cell-type-specific regulatory modules"
)

GO_by_celltype <- list(
  Bun = list(
    Bun_module1_GO, Bun_module2_GO, Bun_module3_GO,
    Bun_module4_GO, Bun_module5_GO
  ),
  Guard = list(
    Guard_module1_GO, Guard_module2_GO, Guard_module3_GO,
    Guard_module4_GO, Guard_module5_GO
  ),
  Mesophyll = list(
    Mesophyll_module1_GO, Mesophyll_module2_GO, Mesophyll_module3_GO,
    Mesophyll_module4_GO, Mesophyll_module5_GO
  ),
  Pavement = list(
    Pavement_module1_GO, Pavement_module2_GO, Pavement_module3_GO,
    Pavement_module4_GO, Pavement_module5_GO
  ),
  Subsidiary = list(
    Subsidiary_module1_GO, Subsidiary_module2_GO, Subsidiary_module3_GO,
    Subsidiary_module4_GO, Subsidiary_module5_GO
  )
)
library(dplyr)
library(purrr)

GO_count_celltype <- map_df(
  names(GO_by_celltype),
  function(ct) {
    
    # 合并该细胞类型下所有模块的 GO 结果
    go_df <- map_df(GO_by_celltype[[ct]], function(x) {
      if (is.null(x) || nrow(as.data.frame(x)) == 0) return(NULL)
      as.data.frame(x)
    })
    
    if (nrow(go_df) == 0) {
      return(
        data.frame(
          CellType = ct,
          Category = c("BP", "non-BP"),
          Count = 0
        )
      )
    }
    
    # 去重（按 GO ID）
    go_df <- distinct(go_df, ID, .keep_all = TRUE)
    
    data.frame(
      CellType = ct,
      Category = c("BP", "non-BP"),
      Count = c(
        sum(go_df$ONTOLOGY == "BP"),
        sum(go_df$ONTOLOGY != "BP")
      )
    )
  }
)
ggplot(GO_count_celltype,
       aes(x = CellType, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c("BP" = "#4C72B0", "non-BP" = "#DCE3F1")
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(
      color = "grey80",
      linewidth = 0.5
    ),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )
