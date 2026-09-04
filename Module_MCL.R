#单细胞测序数据聚类
#聚类
# Module identification from BRCA-related miRNA sponge interaction network identified by integrative method
library(miRspongeR)
set.seed(123)
#infomap
Bun_Cluster_infomap <- netModule(KBoost_Bun, method = "infomap")
Guard_Cluster_infomap <- netModule(KBoost_Guard, method = "infomap")
Mesophyll_Cluster_infomap <- netModule(KBoost_Mesophyll, method = "infomap")
Pavement_Cluster_infomap <- netModule(KBoost_Pavement, method = "infomap")
Subsidiary_Cluster_infomap <- netModule(KBoost_Subsidiary, method = "infomap")

#louvain
Bun_Cluster_louvain <- netModule(KBoost_Bun, method = "louvain")
Guard_Cluster_louvain <- netModule(KBoost_Guard, method = "louvain")
Mesophyll_Cluster_louvain <- netModule(KBoost_Mesophyll, method = "louvain")
Pavement_Cluster_louvain <- netModule(KBoost_Pavement, method = "louvain")
Subsidiary_Cluster_louvain <- netModule(KBoost_Subsidiary, method = "louvain")

#walktrap
Bun_Cluster_walktrap <- netModule(KBoost_Bun, method = "walktrap")
Guard_Cluster_walktrap <- netModule(KBoost_Guard, method = "walktrap")
Mesophyll_Cluster_walktrap <- netModule(KBoost_Mesophyll, method = "walktrap")
Pavement_Cluster_walktrap <- netModule(KBoost_Pavement, method = "walktrap")
Subsidiary_Cluster_walktrap <- netModule(KBoost_Subsidiary, method = "walktrap")

#prop
Bun_Cluster_prop <- netModule(KBoost_Bun, method = "prop")
Guard_Cluster_prop <- netModule(KBoost_Guard, method = "prop")
Mesophyll_Cluster_prop <- netModule(KBoost_Mesophyll, method = "prop")
Pavement_Cluster_prop <- netModule(KBoost_Pavement, method = "prop")
Subsidiary_Cluster_prop <- netModule(KBoost_Subsidiary, method = "prop")

#FN
Bun_Cluster_FN <- netModule(KBoost_Bun, method = "FN")
Guard_Cluster_FN <- netModule(KBoost_Guard, method = "FN")
Mesophyll_Cluster_FN <- netModule(KBoost_Mesophyll, method = "FN")
Pavement_Cluster_FN <- netModule(KBoost_Pavement, method = "FN")
Subsidiary_Cluster_FN <- netModule(KBoost_Subsidiary, method = "FN")

#eigen
Bun_Cluster_eigen <- netModule(KBoost_Bun, method = "eigen")
Guard_Cluster_eigen <- netModule(KBoost_Guard, method = "eigen")
Mesophyll_Cluster_eigen <- netModule(KBoost_Mesophyll, method = "eigen")
Pavement_Cluster_eigen <- netModule(KBoost_Pavement, method = "eigen")
Subsidiary_Cluster_eigen <- netModule(KBoost_Subsidiary, method = "eigen")

#bulk数据
#infomp(573)
leaf_Cluster_infomap <- netModule(GENIE3_leaf, method = "infomap")
#louvain(7)
leaf_Cluster_louvain <- netModule(GENIE3_leaf, method = "louvain")
#walktrap
leaf_Cluster_walktrap <- netModule(GENIE3_leaf, method = "walktrap")
#prop(1)
leaf_Cluster_prop <- netModule(GENIE3_leaf, method = "prop")
#eign(2)
leaf_Cluster_eigen <- netModule(GENIE3_leaf, method = "eigen")
#FN
leaf_Cluster_FN <- netModule(GENIE3_leaf, method = "FN")
#walktrap
