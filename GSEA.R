library(dplyr)
library(eoffice)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(Pi)

#clusterProfiler

library(msigdf)
H <- msigdf.human %>% filter(category_code == "hallmark") %>% select(geneset, symbol) %>% as.data.frame
c2 <- msigdf.human %>% filter(category_code == "c2") %>% select(geneset, symbol) %>% as.data.frame
c2$set<-substr(c2$geneset,1,4)
kegg<-filter(c2,set=='KEGG')%>%select(geneset, symbol) %>% as.data.frame
genset<-rbind(kegg,H)
#range
gsea_data<-select(DEG_IRGPI,symbol,logFC)
gsea_data.sort<-arrange(gsea_data, desc(logFC))
gsea_data<-gsea_data.sort$logFC
names(gsea_data)<-gsea_data.sort$symbol
head(gsea_data)
#gsea hallmark
gsea_hkegg <- GSEA(gsea_data, TERM2GENE = H, 
                   nPerm = 20000, minGSSize = 10, maxGSSize = 500,
                   verbose = TRUE, seed = FALSE, pvalueCutoff = 0.05,by = "DOSE",)
gsea_hkegg_order<-gsea_hkegg[order(gsea_hkegg$enrichmentScore,decreasing=T)]
gseaplot2(gsea_hkegg, row.names(gsea_hkegg_order)[1:5],base_size = 10)

nrow(gsea_hkegg_order)
gseaplot2(gsea_hkegg, row.names(gsea_hkegg_order)[24:28],base_size = 10)

#gsea in xPierGSEA

gene_df<- bitr(rownames(DEG_IRGPI),
               fromType = "ENSEMBL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
DEG_IRGPI$ENSEMBL<-rownames(DEG_IRGPI)
gene_df<-merge(gene_df,DEG_IRGPI,by="ENSEMBL")
gsea_data_entrezi<-select(gene_df,ENTREZID,logFC)
gsea_data_entrezi.sort<-arrange(gsea_data_entrezi, desc(logFC))
gsea_data_entrezi<-gsea_data_entrezi.sort$logFC
names(gsea_data_entrezi)<-gsea_data_entrezi.sort$ENTREZID
head(gsea_data_entrezi)



xPiergsea_data<-gsea_data.sort
xPiergsea_data<-xPiergsea_data[!duplicated(xPiergsea_data$symbol), ]
xPiergsea_data<-mutate(xPiergsea_data,rank = rank(-logFC), priority = logFC)
rownames(xPiergsea_data)<-xPiergsea_data$symbol
xPiergsea_data<-select(xPiergsea_data,priority,rank)

GSEA_pi_h<-
  Pi::xPierGSEA(
    xPiergsea_data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = "C:/Users/chnjqh/Documents/ontology_Rdata")


Pi::xGSEAdotplot(
  GSEA_pi_h,
  top = "HALLMARK_INTERFERON_ALPHA_RESPONSE",title="setID",
  leading = 1:30,max.overlaps=30,
  leading.edge.only = T,leading.arrow=F,leading.force=0.001,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  clab = 'logFC\n(IRGPI-high - IRGPI-low)',subtitle="both",
  signature = F
)+theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 10),
  plot.subtitle = element_text(hjust = 0.5, size = 9)
)

gsea_x_order<-GSEA_pi_h$df_summary
gsea_x_order<-arrange(gsea_x_order,es)

pi_gseaplot_up<-function(x){Pi::xGSEAdotplot(
  GSEA_pi_h,
  top = x, 
  title = 'setID',
  subtitle = 'none',
  compact = T,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  #clab = S'logFC\n(IRGPI-high - IRGPI-low)',
  signature = FALSE
)+ theme(plot.title = element_text(hjust = 1, vjust = -1,size = 10))}
up<-lapply(gsea_x_order$setID[49:45],pi_gseaplot_up)
gridExtra::grid.arrange(grobs = up, ncol = 1)


Pi::xGSEAdotplot(
  GSEA_pi_h,
  top = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",title="setID",
  leading = T,max.overlaps=30,
  leading.edge.only =F,leading.arrow=F,leading.force=0.003,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  clab = 'logFC\n(IRGPI-high - IRGPI-low)',subtitle="both",
  signature = F
)+theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 10),
  plot.subtitle = element_text(hjust = 0.5, size = 9)
)

pi_gseaplot_de<-function(x){Pi::xGSEAdotplot(
  GSEA_pi_h,
  top = x, 
  title = 'setID',
  subtitle = 'none',
  compact = T,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  #clab = 'logFC\n(IRGPI-high - IRGPI-low)',
  signature = FALSE
)+ theme(plot.title = element_text(hjust = 0, vjust = -23,size = 10))}
de<-lapply(gsea_x_order$setID[1:5],pi_gseaplot_de)
gridExtra::grid.arrange(grobs = de, ncol = 1)

#gsea in xPierGSEA(IRGPI_low)

load('DEG_IRGPI_low.Rdata')
gsea_data_2<-select(DEG_IRGPI_low,symbol,logFC)
gsea_data.sort_2<-arrange(gsea_data_2, desc(logFC))
xPiergsea_data_2<-gsea_data.sort_2
xPiergsea_data_2<-xPiergsea_data_2[!duplicated(xPiergsea_data_2$symbol), ]
xPiergsea_data_2<-mutate(xPiergsea_data_2,rank = rank(-logFC), priority = logFC)
rownames(xPiergsea_data_2)<-xPiergsea_data_2$symbol
xPiergsea_data_2<-select(xPiergsea_data_2,priority,rank)

GSEA_pi_h_2<-
  Pi::xPierGSEA(
    xPiergsea_data_2,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = "C:/Users/chnjqh/Documents/ontology_Rdata")


Pi::xGSEAdotplot(
  GSEA_pi_h_2,
  top = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",title="setID",
  leading = 1:30,max.overlaps=30,
  leading.edge.only = T,leading.arrow=F,leading.force=0.001,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  clab = 'logFC\n(IRGPI-low - IRGPI-high)',subtitle="both",
  signature = F
)+theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 10),
  plot.subtitle = element_text(hjust = 0.5, size = 9)
)

gsea_x_order_2<-GSEA_pi_h_2$df_summary
gsea_x_order_2<-arrange(gsea_x_order_2,es)

pi_gseaplot_up_2<-function(x){Pi::xGSEAdotplot(
  GSEA_pi_h_2,
  top = x, 
  title = 'setID',
  subtitle = 'none',
  compact = T,
  colormap = "spectral",
  zlim = c(-1,3),
  peak.color = 'black',
  #clab = S'logFC\n(IRGPI-high - IRGPI-low)',
  signature = FALSE
)+ theme(plot.title = element_text(hjust = 1, vjust = -1,size = 10))}
up_2<-lapply(gsea_x_order_2$setID[49:45],pi_gseaplot_up_2)
gridExtra::grid.arrange(grobs = up_2, ncol = 1)
