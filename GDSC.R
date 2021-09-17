library(pRRophetic)
library(ggplot2)
library(cowplot)
dat <- read.table("rnaExpr_ESCC_symbol.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(dat)[1]<-'gene'
dat<-dat[!duplicated(dat$gene),]
rownames(dat)<-dat$gene
dat<-dat[,-1]
load('IRGPI_ESCC.Rdata')
IRGPI<-clusterProfiler::select(ESCC_IRGPI,IRGPI)
head(IRGPI)
table(IRGPI$IRGPI)
ann<-IRGPI

GCP.drug <- read.table("drug.txt") 
GCP.drug <- GCP.drug$V1
jco <- c("#DC0000B2", "#4DBBD5B2")

GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() 
plotp <- list()
for (drug in GCP.drug) {
  set.seed(1248103) 
  cat(drug," starts!\n") 
  
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) 
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} 
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "IRGPI"=ifelse(ann$IRGPI == "High","IRGPI-High","IRGPI-Low"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$IRGPI <- factor(predictedBoxdat[[drug]]$IRGPI,levels = c("IRGPI-High","IRGPI-Low"),ordered = T) 
  
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=IRGPI, y=est.ic50))
  p <- p + theme_classic()+geom_boxplot(aes(fill = IRGPI)) +
    scale_fill_manual(values = jco[1:length(unique(ann$IRGPI))]) + 
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") +theme(axis.text.x = element_blank())+
    ggtitle(drug) 
  
  plotp[[drug]] <- p 
  cat(drug," has been finished!\n")
}

p_2 <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$IRGPI %in% "IRGPI-High"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$IRGPI %in% "IRGPI-Low"),"est.ic50"]),alternative = "less")$p.value
  p_2 <- append(p_2,tmp) 
}
names(p_2) <- GCP.drug
print(p_2)
write.table(p_2,"output_pvalue.txt", quote = F, sep = "\t")

drug_select<- read.table("drug_select.txt",sep = "\t",header = F,stringsAsFactors = F,check.names = F)
drug_select <- drug_select$V1
plot_grid(plotlist=plotp[drug_select], nrow =2)
