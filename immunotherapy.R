library(pheatmap)
dat <- read.table("expr_escc.txt", header=T, sep="\t", check.names=F)
load('IRGPI_ESCC.Rdata')
IRGPI<-clusterProfiler::select(ESCC_IRGPI,IRGPI)
head(IRGPI)
table(IRGPI$IRGPI)
ann<-IRGPI

#submap
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

skcm.immunotherapy.logNC <- read.table("skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) 
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) 
skcm.immunotherapy.info <- read.table("skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
tmp <- read.table("rna_escc.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T) 
colnames(tmp)[1]<-'gene'
tmp<-tmp[!duplicated(tmp$gene),]
rownames(tmp)<-tmp$gene
tmp<-tmp[,-1]
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) 

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

samples.High <- rownames(ann[which(ann$IRGPI == "High"),])
samples.Low <- rownames(ann[which(ann$IRGPI == "Low"),])

sam_info <- data.frame("IRGPI"=c(samples.High,samples.Low),row.names = c(samples.High,samples.Low))
sam_info$rank <- rep(c(1,2),times=c(length(samples.High),length(samples.Low))) 

gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) 
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

#heatmap
library(eoffice)
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"

tmp <- matrix(c(0.933,0.319,0.996,0.016,0.044,0.209,0.122,0.997, 
                1,1,1,0.135,0.359,1,0.975,1), 
              nrow = 4,byrow = T,dimnames = list(c("IRGPI-High_p","IRGPI-Low_p","IRGPI-High_b","IRGPI-Low_b"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry))
         )