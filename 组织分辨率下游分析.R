setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析")
#筛选组织特异性GRN关键TFs
save(GENIE3_leaf,file='GENIE3_leaf.Rdata')
library(dplyr)
regulation_count_leaf<- GENIE3_leaf %>%
  group_by(TF) %>%        # TF 是转录因子列名
  summarise(n_targets = n())  # 统计每个转录因子的靶基因数
# 计算所有转录因子的平均调控数
leaf_average_targets <- mean(regulation_count_leaf$n_targets)
# 使用泊松分布计算每个转录因子的p值
leaf_ppois <- regulation_count_leaf %>%
  mutate(p_value = 1 - ppois(n_targets - 1, lambda = leaf_average_targets))  # -1 因为ppois是≤的累积分布函数
# 设置显著性水平，比如0.05
significance_level <- 0.05
# 筛选出p值小于显著性水平的关键转录因子
leaf_key_TFs <- leaf_ppois %>%
  filter(p_value < significance_level)
library(dplyr)
library(ggplot2)

# 取度中心性最高的前10个TF
top10_TFs <- regulation_count_leaf %>%
  arrange(desc(n_targets)) %>%
  slice(1:10)

# 按度中心性从高到低排序
top10_TFs$TF <- factor(top10_TFs$TF, levels = rev(top10_TFs$TF))

# 绘制棒棒糖图
ggplot(top10_TFs, aes(x = TF, y = n_targets)) +
  geom_segment(aes(x = TF, xend = TF, y = 0, yend = n_targets),
               linewidth = 1.2, color = "grey60") +
  geom_point(size = 6, color = "#fcb8c9") +
  coord_flip() +
  labs(x = "Transcription Factors",
       y = "Degree Centrality (Number of Targets)",
       title = "Top 10 TFs with Highest Degree Centrality in GENIE3 GRN") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
#计算每个靶基因的入度（被多少TF调控）
target_indegree <- GENIE3_leaf %>%
  group_by(Target) %>%
  summarise(n_TFs = n())
# 去掉指定基因
top_targets_filtered <- target_indegree %>%
  filter(Target != "6:108383457-108389721")

# 重新排序并取前10
top10_targets <- top_targets_filtered %>%
  arrange(desc(n_TFs)) %>%
  slice(1:10)

# 设置排序顺序（用于绘图从高到低）
top10_targets$Target <- factor(top10_targets$Target, levels = rev(top10_targets$Target))
# 按入度从高到低排序
top10_targets$Target <- factor(top10_targets$Target, levels = rev(top10_targets$Target))
# 绘制棒棒糖图
ggplot(top10_targets, aes(x = Target, y = n_TFs)) +
  geom_segment(aes(x = Target, xend = Target, y = 0, yend = n_TFs),
               linewidth = 1.2, color = "grey60") +
  geom_point(size = 6, color = "#a6d0f1") +
  coord_flip() +
  labs(x = "Target Genes",
       y = "In-degree (Number of TFs)",
       title = "Top 10 Target Genes with Highest In-degree in GENIE3 GRN") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
#绘制关键TFs柱状图
library(ggplot2)

tf_data <- data.frame(
  Category = c("Total TFs", "Key TFs"),
  Count = c(1526, 462)
)

ggplot(tf_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.3) +
  labs(x = "", y = "Number of TFs",
       title = "Identification of Key Transcription Factors") +
  scale_fill_manual(values = c("Total TFs" = "#89c9c8",
                               "Key TFs" = "#f9bebb")) +
  theme_bw() +
  theme(legend.position = "none")
#模块识别
library(miRspongeR)
set.seed(123)
#infomap
leaf_Cluster_infomap <- netModule(GENIE3_leaf,modulesize = 10, method = "infomap")
#FN
leaf_Cluster_FN <- netModule(GENIE3_leaf,modulesize = 10, method = "FN")
#eign
leaf_Cluster_eigen <- netModule(GENIE3_leaf,modulesize = 10, method = "eigen")
#louvain
leaf_Cluster_louvain <- netModule(GENIE3_leaf,modulesize = 10, method = "louvain")
#walktrap
leaf_Cluster_walktrap <- netModule(GENIE3_leaf,modulesize = 10, method = "walktrap")
#prop
library(miRspongeR)
set.seed(123)
leaf_Cluster_prop <- netModule(GENIE3_leaf,modulesize = 10, method = "prop")
leaf_Cluster_MCL <- netModule(GENIE3_leaf,modulesize = 10, method = "MCL")
leaf_Cluster_MCODE <- netModule(GENIE3_leaf,modulesize = 10, method = "MCODE")
leaf_Cluster_walktrap <- netModule(GENIE3_leaf,modulesize = 10, method = "walktrap")

library(ggplot2)

library(ggplot2)

module_method <- data.frame(
  Method = c("Infomap", "Walktrap", "MCL", "Louvain", "FN", "Eigen", "Prop", "MCODE"),
  ModuleNumber = c(477, 86, 7, 6, 3, 2, 1, 1)
)
module_method$Method <- factor(
  module_method$Method,
  levels = module_method$Method[order(module_method$ModuleNumber, decreasing = TRUE)]
)
ggplot(module_method, aes(x = Method, y = ModuleNumber, fill = Method)) +
  geom_col(width = 0.7) +
  
  geom_text(
    aes(label = ModuleNumber),
    vjust = -0.3,
    size = 4
  ) +
  
  scale_y_log10() +   # ⭐ 使用对数坐标
  
  scale_fill_manual(values = c(
    "#b5aed5",
    "#b2e6fd",
    "#b8d2cc",
    "#e8b2a7",
    "#feebb9",
    "#f4c2c2",
    "#c7e9c0",
    "#fdd49e"
  )) +
  
  labs(
    x = "Community detection method",
    y = "Number of modules (log10 scale)",
    title = "Comparison of module numbers identified by different community detection methods"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  
  expand_limits(y = max(module_method$ModuleNumber) * 1.3)
#统计各模块中关键TFs数量
key_TFs <- leaf_ppois %>%
  filter(p_value < significance_level) %>%
  pull(TF)
module_keyTF_count <- sapply(leaf_Cluster_infomap, function(genes){
  sum(genes %in% key_TFs)
})
module_keyTF_df <- data.frame(
  Module = paste0("Module_", seq_along(module_keyTF_count)),
  KeyTF_number = module_keyTF_count
)
library(dplyr)

module_gene_number <- sapply(leaf_Cluster_infomap, length)
module_keyTF_df$GeneNumber <- module_gene_number
top_modules <- module_keyTF_df %>%
  arrange(desc(KeyTF_number), desc(GeneNumber)) %>%
  slice(1:10)

library(dplyr)
library(tidyr)
library(ggplot2)

plot_data <- top_modules %>%
  select(Module, GeneNumber, KeyTF_number) %>%
  pivot_longer(
    cols = c(GeneNumber, KeyTF_number),
    names_to = "Type",
    values_to = "Count"
  )
ggplot(plot_data, aes(x = Module, y = Count, fill = Type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3.5
  ) +
  
  scale_fill_manual(values = c(
    "GeneNumber" = "#91abd2",
    "KeyTF_number" = "#f59694"
  )) +
  
  labs(
    x = "Module",
    y = "Number of genes",
    fill = "",
    title = "Comparison of total genes and key TFs in top modules"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#对关键TFs最多的前五个模块进行GO富集分析
library(openxlsx)
# 提取第 60个模块的基因
Leaf_module60<- leaf_Cluster_infomap[[60]]
# 转为数据框
Leaf_module60 <- data.frame(Gene = Leaf_module60)
Leaf_module60 <- Leaf_module60[grepl("^Zm", Leaf_module60[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析")
# 保存为 Excel
write.xlsx(
  Leaf_module60,
  file = "Leaf_module60.xlsx",
  rowNames = FALSE
)
# 提取第26个模块的基因
Leaf_module26<- leaf_Cluster_infomap[[26]]
# 转为数据框
Leaf_module26 <- data.frame(Gene = Leaf_module26)
Leaf_module26 <- Leaf_module26[grepl("^Zm", Leaf_module26[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析")
# 保存为 Excel
write.xlsx(
  Leaf_module26,
  file = "Leaf_module26.xlsx",
  rowNames = FALSE
)
# 提取第315个模块的基因
Leaf_module315<- leaf_Cluster_infomap[[315]]
# 转为数据框
Leaf_module315 <- data.frame(Gene = Leaf_module315)
Leaf_module315 <- Leaf_module315[grepl("^Zm", Leaf_module315[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析")
# 保存为 Excel
write.xlsx(
  Leaf_module315,
  file = "Leaf_module315.xlsx",
  rowNames = FALSE
)
# 提取第276个模块的基因
Leaf_module276<- leaf_Cluster_infomap[[276]]
# 转为数据框
Leaf_module276 <- data.frame(Gene = Leaf_module276)
Leaf_module276 <- Leaf_module276[grepl("^Zm", Leaf_module276[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析")
# 保存为 Excel
write.xlsx(
  Leaf_module276,
  file = "Leaf_module276.xlsx",
  rowNames = FALSE
)
# 提取第289个模块的基因
Leaf_module289<- leaf_Cluster_infomap[[289]]
# 转为数据框
Leaf_module289 <- data.frame(Gene = Leaf_module289)
Leaf_module289 <- Leaf_module289[grepl("^Zm", Leaf_module289[[1]]), , drop = FALSE]
setwd("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析")
# 保存为 Excel
write.xlsx(
  Leaf_module289,
  file = "Leaf_module289.xlsx",
  rowNames = FALSE
)
library(readxl)
library(clusterProfiler)
library("AnnotationHub")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
hub <- AnnotationHub()
query(hub, "zea mays")
zeamaize <- hub[['AH114308']]
#使用maizeGDB数据库进行基因id转换
#GO富集分析
#模块26
Module_26<- read_excel("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析/Module26.xlsx")
Module_26 <- Module_26$`Gene`
module26_GO <- enrichGO(gene = Module_26,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
#模块60
Module_60<- read_excel("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析/Module60.xlsx")
Module_60 <- Module_60$`Gene`
module60_GO <- enrichGO(gene = Module_60,#我们上面定义了
                        OrgDb=zeamaize,
                        keyType = "ENTREZID",
                        ont = "ALL",#富集的GO类型
                        pAdjustMethod = "BH",#这个不用管，一般都用的BH
                        minGSSize = 1,
                        pvalueCutoff = 0.05,#P值可以取0.05
                        qvalueCutoff = 0.05,
                        readable = TRUE)
#模块276
Module_276<- read_excel("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析/Module276.xlsx")
Module_276 <- Module_276$`Gene`
module276_GO <- enrichGO(gene = Module_276,#我们上面定义了
                        OrgDb=zeamaize,
                        keyType = "ENTREZID",
                        ont = "ALL",#富集的GO类型
                        pAdjustMethod = "BH",#这个不用管，一般都用的BH
                        minGSSize = 1,
                        pvalueCutoff = 0.05,#P值可以取0.05
                        qvalueCutoff = 0.05,
                        readable = TRUE)
#模块289
Module_289<- read_excel("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析/Module289.xlsx")
Module_289 <- Module_289$`Gene`
module289_GO <- enrichGO(gene = Module_289,#我们上面定义了
                           OrgDb=zeamaize,
                           keyType = "ENTREZID",
                           ont = "ALL",#富集的GO类型
                           pAdjustMethod = "BH",#这个不用管，一般都用的BH
                           minGSSize = 1,
                           pvalueCutoff = 0.05,#P值可以取0.05
                           qvalueCutoff = 0.05,
                           readable = TRUE)
# 将GO富集结果转换为数据框
dt2 <- as.data.frame(module289_GO)

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

#模块315
Module_315<- read_excel("D:/玉米多模态数据的基因调控网络研究/bulk数据下游分析/GO富集分析/Module315.xlsx")
Module_315 <- Module_315$`Gene`
module315_GO <- enrichGO(gene = Module_315,#我们上面定义了
                         OrgDb=zeamaize,
                         keyType = "ENTREZID",
                         ont = "ALL",#富集的GO类型
                         pAdjustMethod = "BH",#这个不用管，一般都用的BH
                         minGSSize = 1,
                         pvalueCutoff = 0.05,#P值可以取0.05
                         qvalueCutoff = 0.05,
                         readable = TRUE)

# 将GO富集结果转换为数据框
dt2 <- as.data.frame(module315_GO)

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
