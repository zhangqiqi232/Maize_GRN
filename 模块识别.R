#聚类
# Module identification from BRCA-related miRNA sponge interaction network identified by integrative method
setwd("D:/玉米多模态数据的基因调控网络研究/模块识别")
library(readxl)
set.seed(123)
# 读取六个基因调控网络
linklist_Bun <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Bundle.xlsx")
linklist_guard <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Guard.xlsx")
linklist_Mesophyll <- read.csv("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Mesophyll.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
linklist_Pavement <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/pavement.xlsx")
linklist_Subsidiary <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/subsidiary.xlsx")
linklist_tissue <- read.csv("D:/玉米多模态数据的基因调控网络研究/玉米叶片bulk/推断GRN/GRNBoost2/network_leaf_output.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

library(miRspongeR)
linklist_Bun <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Bundle.xlsx")
linklist_Bun <- linklist_Bun[ , -3]
Bundle_Cluster_FN <- netModule(linklist_Bun, method = "FN",modulesize = 10)
Bundle_Cluster_MCL <- netModule(linklist_Bun, method = "MCL", modulesize = 10)
Bundle_Cluster_LINKCOMM <- netModule(linklist_Bun, method = "LINKCOMM", modulesize = 10)
Bundle_Cluster_MCODE <- netModule(linklist_Bun, method = "MCODE", modulesize = 10)
Bundle_Cluster_louvain <- netModule(linklist_Bun, method = "louvain", modulesize = 10)
Bundle_Cluster_betweenness <- netModule(linklist_Bun, method = "betweenness", modulesize = 10)
Bundle_Cluster_infomap <- netModule(linklist_Bun, method = "infomap", modulesize = 10)
save(Bundle_Cluster_infomap,file='Bundle_Cluster_infomap.Rdata')
Bundle_Cluster_prop <- netModule(linklist_Bun, method = "prop", modulesize = 10)
Bundle_Cluster_eigen <- netModule(linklist_Bun, method = "eigen", modulesize = 10)
save(Bundle_Cluster_louvain,file='Bundle_Cluster_louvain.Rdata')

#读入Guard数据的基因调控网络
linklist_Guard <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Guard.xlsx")
linklist_Guard <- linklist_Guard[ , -3]
Guard_Cluster_MCL <- netModule(linklist_Guard, method = "MCL",modulesize = 10)
Guard_Cluster_LINKCOMM<- netModule(linklist_Guard, method = "LINKCOMM",modulesize = 10)
Guard_Cluster_FN <- netModule(linklist_Guard, method = "FN",modulesize = 10)
Guard_Cluster_MCODE <- netModule(linklist_Guard, method = "MCODE",modulesize = 10)
Guard_Cluster_louvain<- netModule(linklist_Guard, method = "louvain",modulesize = 10)
Guard_Cluster_infomap<- netModule(linklist_Guard, method = "infomap",modulesize = 10)
save(Guard_Cluster_infomap,file='Guard_Cluster_infomap.Rdata')
Guard_Cluster_prop<- netModule(linklist_Guard, method = "prop",modulesize = 10)
Guard_Cluster_eigen<- netModule(linklist_Guard, method = "eigen",modulesize = 10)
save(Guard_Cluster_louvain,file='Guard_Cluster_louvain.Rdata')

#读入Pavement数据的基因调控网络
linklist_Pavement <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/pavement.xlsx")
linklist_Pavement <- linklist_Pavement[ , -3]
Pavement_Cluster_MCL <- netModule(linklist_Pavement, method = "MCL",modulesize = 10)
Pavement_Cluster_LINKCOMM<- netModule(linklist_Pavement, method = "LINKCOMM",modulesize = 10)
Pavement_Cluster_FN <- netModule(linklist_Pavement, method = "FN",modulesize = 10)
Pavement_Cluster_MCODE <- netModule(linklist_Pavement, method = "MCODE",modulesize = 10)
Pavement_Cluster_louvain<- netModule(linklist_Pavement, method = "louvain",modulesize = 10)
Pavement_Cluster_infomap<- netModule(linklist_Pavement, method = "infomap",modulesize = 10)
save(Pavement_Cluster_infomap,file='Pavement_Cluster_infomap.Rdata')
Pavement_Cluster_prop<- netModule(linklist_Pavement, method = "prop",modulesize = 10)
Pavement_Cluster_eigen<- netModule(linklist_Pavement, method = "eigen",modulesize = 10)
save(Pavement_Cluster_louvain,file='Pavement_Cluster_louvain.Rdata')

#读入Subsidiary数据的基因调控网络
linklist_Subsidiary <- read_excel("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/subsidiary.xlsx")
linklist_Subsidiary <- linklist_Subsidiary[ , -3]
Subsidiary_Cluster_MCL <- netModule(linklist_Subsidiary, method = "MCL",modulesize = 10)
Subsidiary_Cluster_LINKCOMM<- netModule(linklist_Subsidiary, method = "LINKCOMM",modulesize = 10)
Subsidiary_Cluster_FN<- netModule(linklist_Subsidiary, method = "FN",modulesize = 10)
Subsidiary_Cluster_louvain<- netModule(linklist_Subsidiary, method = "louvain",modulesize = 10)
Subsidiary_Cluster_walktrap<- netModule(linklist_Subsidiary, method = "walktrap",modulesize = 10)
Subsidiary_Cluster_infomap<- netModule(linklist_Subsidiary, method = "infomap",modulesize = 10)
save(Subsidiary_Cluster_infomap,file='Subsidiary_Cluster_infomap.Rdata')
Subsidiary_Cluster_prop<- netModule(linklist_Subsidiary, method = "prop",modulesize = 10)
Subsidiary_Cluster_eigen<- netModule(linklist_Subsidiary, method = "eigen",modulesize = 10)
save(Subsidiary_Cluster_louvain,file='Subsidiary_Cluster_louvain.Rdata')

#读入Mesophyll数据的基因调控网络
linklist_Mesophyll <- read.csv("D:/玉米多模态数据的基因调控网络研究/GRNBoost2/Mesophyll.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
linklist_Mesophyll <- linklist_Mesophyll[ , -3]
Mesophyll_Cluster_MCL <- netModule(linklist_Mesophyll, method = "MCL",modulesize = 10)
Mesophyll_Cluster_LINKCOMM<- netModule(linklist_Mesophyll, method = "LINKCOMM",modulesize = 10)
Mesophyll_Cluster_FN<- netModule(linklist_Mesophyll, method = "FN",modulesize = 10)
Mesophyll_Cluster_louvain<- netModule(linklist_Mesophyll, method = "louvain",modulesize = 10)
Mesophyll_Cluster_infomap<- netModule(linklist_Mesophyll, method = "infomap",modulesize = 10)
save(Mesophyll_Cluster_infomap,file='Mesophyll_Cluster_infomap.Rdata')
Mesophyll_Cluster_prop<- netModule(linklist_Mesophyll, method = "prop",modulesize = 10)
Mesophyll_Cluster_eigen<- netModule(linklist_Mesophyll, method = "eigen",modulesize = 10)
save(Mesophyll_Cluster_louvain,file='Mesophyll_Cluster_louvain.Rdata')

#分别获取每个细胞类型的模块
#Bundle sheath
Bundle_module1 <- data.frame(Bundle_Cluster_infomap[[1]]) 
write.xlsx(Bundle_module1, "Bundle_module1.xlsx")
Bundle_module2 <- data.frame(Bundle_Cluster_louvain[[2]]) 
write.xlsx(Bundle_module2, "Bundle_module2.xlsx")
Bundle_module3 <- data.frame(Bundle_Cluster_louvain[[3]]) 
write.xlsx(Bundle_module3, "Bundle_module3.xlsx")
Bundle_module4 <- data.frame(Bundle_Cluster_louvain[[4]]) 
write.xlsx(Bundle_module4, "Bundle_module4.xlsx")
Bundle_module5 <- data.frame(Bundle_Cluster_louvain[[5]]) 
write.xlsx(Bundle_module5, "Bundle_module5.xlsx")
Bundle_module6 <- data.frame(Bundle_Cluster_louvain[[6]]) 
write.xlsx(Bundle_module6, "Bundle_module6.xlsx")
Bundle_module7 <- data.frame(Bundle_Cluster_louvain[[7]]) 
write.xlsx(Bundle_module7, "Bundle_module7.xlsx")
Bundle_module8 <- data.frame(Bundle_Cluster_louvain[[8]]) 
write.xlsx(Bundle_module8, "Bundle_module8.xlsx")
Bundle_module9 <- data.frame(Bundle_Cluster_louvain[[9]]) 
write.xlsx(Bundle_module9, "Bundle_module9.xlsx")
Bundle_module10 <- data.frame(Bundle_Cluster_louvain[[10]]) 
write.xlsx(Bundle_module10, "Bundle_module10.xlsx")
#对每一个模块进行GO富集分析
#module1
ensembl <- useMart("plants_mart", dataset = "zmays_eg_gene", host = "https://plants.ensembl.org")
# 你的基因 ID 列表
gene_ids <- c("Zm00001d014812", "Zm00001d035186", "Zm00001d031258")  # 添加更多 ID
# 获取对应的 Entrez Gene ID
converted <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "ensembl_gene_id",
                   values = gene_ids,
                   mart = ensembl)
query(hub, "zea")









#Guard
Guard_module1 <- data.frame(Guard_Cluster_louvain[[1]]) 
write.xlsx(Guard_module1, "Guard_module1.xlsx")
Guard_module2 <- data.frame(Guard_Cluster_louvain[[2]]) 
write.xlsx(Guard_module2, "Guard_module2.xlsx")
Guard_module3 <- data.frame(Guard_Cluster_louvain[[3]]) 
write.xlsx(Guard_module3, "Guard_module3.xlsx")
Guard_module4 <- data.frame(Guard_Cluster_louvain[[4]]) 
write.xlsx(Guard_module4, "Guard_module4.xlsx")
Guard_module5 <- data.frame(Guard_Cluster_louvain[[5]]) 
write.xlsx(Guard_module5, "Guard_module5.xlsx")
Guard_module6 <- data.frame(Guard_Cluster_louvain[[6]]) 
write.xlsx(Guard_module6, "Guard_module6.xlsx")
Guard_module7 <- data.frame(Guard_Cluster_louvain[[7]]) 
write.xlsx(Guard_module7, "Guard_module7.xlsx")
Guard_module8 <- data.frame(Guard_Cluster_louvain[[8]]) 
write.xlsx(Guard_module8, "Guard_module8.xlsx")
Guard_module9 <- data.frame(Guard_Cluster_louvain[[9]]) 
write.xlsx(Guard_module9, "Guard_module9.xlsx")
Guard_module10 <- data.frame(Guard_Cluster_louvain[[10]]) 
write.xlsx(Guard_module10, "Guard_module10.xlsx")
Guard_module11 <- data.frame(Guard_Cluster_louvain[[11]]) 
write.xlsx(Guard_module11, "Guard_module11.xlsx")
Guard_module12 <- data.frame(Guard_Cluster_louvain[[12]]) 
write.xlsx(Guard_module12, "Guard_module12.xlsx")
Guard_module13 <- data.frame(Guard_Cluster_louvain[[13]]) 
write.xlsx(Guard_module13, "Guard_module13.xlsx")

#Mesophyll
Guard_module13 <- data.frame(Guard_Cluster_louvain[[13]]) 
write.xlsx(Guard_module13, "Guard_module13.xlsx")


linklist_tissue <- linklist_tissue[ , -3]
Tissue_Cluster_infomap<- netModule(linklist_tissue, method = "infomap",modulesize = 10)
Tissue_Cluster_louvain<- netModule(linklist_tissue, method = "louvain",modulesize = 10)


#
收集六个细胞类型的模块列表
module_lists <- list(
  Bundle = Bundle_Cluster_infomap,
  Guard = Guard_Cluster_infomap,
  Pavement = Pavement_Cluster_infomap,
  Subsidiary = Subsidiary_Cluster_infomap,
  Mesophyll = Mesophyll_Cluster_infomap,
  Tissue = Tissue_Cluster_infomap
)
# 为每个模块生成唯一“签名” = 基因排序后 paste 成字符串
get_module_signatures <- function(modules, celltype) {
  signatures <- sapply(modules, function(genes) {
    paste(sort(genes), collapse = "_")
  })
  names(signatures) <- paste(celltype, seq_along(modules), sep = "_")
  return(signatures)
}
# 获取所有细胞类型的模块签名
all_signatures <- lapply(names(module_lists), function(celltype) {
  get_module_signatures(module_lists[[celltype]], celltype)
})
all_signatures <- unlist(all_signatures)
# 构建模块签名到细胞类型的映射
signature_to_celltype <- list()
for (i in seq_along(all_signatures)) {
  sig <- all_signatures[i]
  ct <- strsplit(names(all_signatures)[i], "_")[[1]][1]
  signature_to_celltype[[sig]] <- c(signature_to_celltype[[sig]], ct)
}
# 创建 presence/absence 二元矩阵
all_celltypes <- names(module_lists)
all_sigs <- names(signature_to_celltype)
binary_matrix <- matrix(0, nrow = length(all_sigs), ncol = length(all_celltypes))
rownames(binary_matrix) <- all_sigs
colnames(binary_matrix) <- all_celltypes
for (sig in all_sigs) {
  binary_matrix[sig, signature_to_celltype[[sig]]] <- 1
}
binary_df <- as.data.frame(binary_matrix)
upset(binary_df,
      sets = colnames(binary_df),
      sets.bar.color = "#F5B5D9",
      order.by = "freq",
      keep.order = TRUE,
      mainbar.y.label = "Number of Shared Modules (Identical Gene Set)",
      sets.x.label = "Modules per Cell Type")
write.xlsx(binary_df, "binary_df.xlsx")

upset(
  binary_df,
  sets = colnames(binary_df),
  sets.bar.color = "#F5B5D9",
  order.by = "freq",
  keep.order = TRUE,
  mainbar.y.label = "Number of Shared Modules (Identical Gene Set)",
  sets.x.label = "Modules per Cell Type",
  
  # 以下是新增优化参数：
  text.scale = c(1.8, 1.8, 1.5, 1.5, 1.8, 1.8),  # 增大所有字体
  point.size = 4.5,                            # 加大点的大小
  line.size = 2                                # 加粗连线
)
