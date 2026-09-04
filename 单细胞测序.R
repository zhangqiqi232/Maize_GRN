library(Seurat)
library(dplyr)
library(patchwork)
setwd("D:/玉米多模态数据的基因调控网络研究/单细胞测序")
h5_file <- "D:/玉米多模态数据的基因调控网络研究/原始数据/GSM4110694_filtered_feature_bc_matrix.h5"
zea_data <- Read10X_h5(file = h5_file)
#创建zeamaysleaf对象，并且去除了低质量的细胞和基因
zeamaysleaf <- CreateSeuratObject(counts = zea_data, project = "GSE138256", min.features = 200, min.cells = 3)
zeamaysleaf
#查看前面几行的基因名称，发现该数据中基因名称存在"gene:"的前缀
head(rownames(zeamaysleaf))
#将基因名称的前缀进行去除
rownames(zeamaysleaf) <- sub("^gene:", "", rownames(zeamaysleaf))
#计算玉米线粒体基因的比例
zeamaysleaf[["percent.mt"]] <- PercentageFeatureSet(zeamaysleaf, pattern = "^atp|^nad|^cox|^cob|^rps|^rpl|^mtt")
#小提琴图可视化
VlnPlot(zeamaysleaf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#绘制散点图，说明两两特征之间的关系
plot1 <- FeatureScatter(zeamaysleaf, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zeamaysleaf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#质控（根据绘制出来的小提琴图，仅保留每个细胞中检测到的基因数大于200并且小于5000的基因以及保留线粒体基因占比小于10%的细胞）
zeamaysleaf <- subset(zeamaysleaf, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
#对该数据进行标准化处理
zeamaysleaf <- NormalizeData(zeamaysleaf, normalization.method = "LogNormalize", scale.factor = 10000)                                                                                             
#识别该对象中的高变基因
zeamaysleaf <- FindVariableFeatures(zeamaysleaf, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(zeamaysleaf), 10)
plot1 <- VariableFeaturePlot(zeamaysleaf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#可视化其中变异度最高的前10个基因
plot1+plot2
#提取数据矩阵的行名称
all.genes <- rownames(zeamaysleaf)
#缩放数据
zeamaysleaf<- ScaleData(zeamaysleaf, features = all.genes) 
#只对挑选出来的2000个高变基因进行降维
zeamaysleaf <- RunPCA(zeamaysleaf, features = VariableFeatures(object = zeamaysleaf))
#可视化细胞在两个PC上的分布情况
DimPlot(zeamaysleaf, reduction = "pca")
#评估PCA中每个主成分的显著性，帮助决定哪些主成分应保留以描述数据的主要变化模式
zeamaysleaf <- JackStraw(zeamaysleaf, num.replicate = 100)
#Seurat结合JackStraw程序和置换检验对PC进行显著性分析，鉴定出显著的PC以进行后续分析
zeamaysleaf<- ScoreJackStraw(zeamaysleaf, dims = 1:15)
JackStrawPlot(zeamaysleaf, dims = 1:15)
ElbowPlot(zeamaysleaf)
#根据肘图选择前15个pc进行分析
zeamaysleaf <- FindNeighbors(zeamaysleaf, dims = 1:15)
#进行聚类
zeamaysleaf<- FindClusters(zeamaysleaf, resolution = 0.7) 
zeamaysleaf <- RunUMAP(zeamaysleaf, dims = 1:15)
DimPlot(zeamaysleaf,reduction = "umap")
#找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞
zeamaysleaf.markers <- FindAllMarkers(zeamaysleaf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   zeamaysleaf.markers %>%
  group_by(cluster) %>%
  slice_max(n=10,order_by = avg_log2FC)
zeamaysleaf.markers%>%
  group_by(cluster)%>%
  top_n(n=100,wt=avg_log2FC)->top100
write.csv(top100,file = "top100.csv")
#细胞类型注释
new.cluster.ids <- c("0" = "Mesophyll", "1" = "Pavement", "2" = "Subsidiary", "3" = "Subsidiary",  "4"= "Guard",  "5" = "Mesophyll", "6" = "Pavement" , "7" = "Subsidiary", "8" = "Mesophyll", "9" = "Pavement", "10" = "Bundle sheath", "11" = "Bundle sheath")
names(new.cluster.ids) <- levels(zeamaysleaf)
zeamaysleaf <- RenameIdents(zeamaysleaf, new.cluster.ids)
DimPlot(zeamaysleaf, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
levels(Idents(zeamaysleaf))
save(zeamaysleaf,file='zeamaysleaf.Rdata')
