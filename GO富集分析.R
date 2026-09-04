#使用AnnotationHub来获得玉米的注释数据库
library("AnnotationHub")
hub <- AnnotationHub()
query(hub, "zea mays")
maize <- hub[['AH114308']]
query(hub, "maize")
maize <- hub[['AH115863']]
#使用maizeGDB数据库进行基因id转换
#GO富集分析(Bundle sheath)
gene_bun<- read_excel("D:/玉米多模态数据的基因调控网络研究/玉米叶片bulk/转录因子(Entrez).xlsx")
gene_bun <- gene_bun$`gene`
library(clusterProfiler)
ego_bun <- enrichGO(gene = gene_bun,#我们上面定义了
                    OrgDb=maize,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
bundle_sheath<-barplot(ego_bun, 
                          x = "GeneRatio", 
                          color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
                          showCategory =20, #只显示前20
                          label_format=100,
                          title = "leaf",#设置图片的标题
                          #split="ONTOLOGY") +  #以ONTOLOGY类型分开
                          #theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
                          #facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
)
#GO富集分析(Guard)
gene_guard<- read_excel("D:/玉米多模态数据的基因调控网络研究/GO富集分析/guard_Entrez.xlsx")
gene_guard <- gene_guard$`gene`
ego_guard <- enrichGO(gene = gene_guard,#我们上面定义了
                    OrgDb=maize,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
guard<-barplot(ego_guard, 
                       x = "GeneRatio", 
                       color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
                       showCategory =20, #只显示前20
                       label_format=100,
                       title = "Guard"#设置图片的标题
                       #split="ONTOLOGY") +  #以ONTOLOGY类型分开
                       #theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
                       #facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
)
#GO富集分析(mesophy)
gene_mesophyll<- read_excel("D:/玉米多模态数据的基因调控网络研究/GO富集分析/mesophyll_Entrez.xlsx")
gene_mesophyll <- gene_mesophyll$`gene`
ego_mesophyll <- enrichGO(gene = gene_mesophyll,#我们上面定义了
                      OrgDb=maize,
                      keyType = "ENTREZID",
                      ont = "ALL",#富集的GO类型
                      pAdjustMethod = "BH",#这个不用管，一般都用的BH
                      minGSSize = 1,
                      pvalueCutoff = 0.01,#P值可以取0.05
                      qvalueCutoff = 0.05,
                      readable = TRUE)
mesophyll<-barplot(ego_mesophyll, 
               x = "GeneRatio", 
               color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
               showCategory =20, #只显示前20
               label_format=100,
               title = "Mesophyll",#设置图片的标题
               #split="ONTOLOGY") +  #以ONTOLOGY类型分开
               #theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
               #facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
)
#GO富集分析(pavement)
gene_pavement<- read_excel("D:/玉米多模态数据的基因调控网络研究/GO富集分析/pavement_Entrez.xlsx")
gene_pavement <- gene_pavement$`gene`
ego_pavement <- enrichGO(gene = gene_pavement,#我们上面定义了
                   OrgDb=maize,
                   keyType = "ENTREZID",
                   ont = "ALL",#富集的GO类型
                   pAdjustMethod = "BH",#这个不用管，一般都用的BH
                   minGSSize = 1,
                   pvalueCutoff = 0.01,#P值可以取0.05
                   qvalueCutoff = 0.05,
                   readable = TRUE)
pavement<-barplot(ego_pavement, 
                   x = "GeneRatio", 
                   color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
                   showCategory =20, #只显示前20
                   label_format=100,
                   title = "Pavement",#设置图片的标题
                   #split="ONTOLOGY") +  #以ONTOLOGY类型分开
                   #theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
                   #facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
)
#GO富集分析(subsidiary)
gene_subsidiary<- read_excel("D:/玉米多模态数据的基因调控网络研究/GO富集分析/subsidiary_Entrez.xlsx")
gene_subsidiary <- gene_subsidiary$`gene`
ego_subsidiaty <- enrichGO(gene = gene_subsidiary ,#我们上面定义了
                         OrgDb=maize,
                         keyType = "ENTREZID",
                         ont = "BP",#富集的GO类型
                         pAdjustMethod = "BH",#这个不用管，一般都用的BH
                         minGSSize = 1,
                         pvalueCutoff = 0.01,#P值可以取0.05
                         qvalueCutoff = 0.05,
                         readable = TRUE)
subsidiary<-barplot(ego_subsidiaty, 
            x = "GeneRatio", 
            color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
            showCategory =20, #只显示前20
            label_format=100,
            title = "Subsidiary",#设置图片的标题
            #split="ONTOLOGY") +  #以ONTOLOGY类型分开
            #theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
            #facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
)
gene<- read_excel("C:/Users/张琦/Desktop/guard模块.xlsx")
#可视化气泡图
dotplot(ego_ALL, showCategory = 20,label_fomat = 30, title = "GO Enrichment Analysis Bubble Plot")
#可视化气泡图
dotplot(ego_ALL,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 35,
        size = NULL,
        split = NULL,
        font.size = 12,
        title="",
        orderBy = "x",
        label_format = 60)#方
