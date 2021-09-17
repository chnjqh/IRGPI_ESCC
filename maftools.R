library(maftools)
library(dplyr)
library(ggsci)
library(eoffice)
rt <- read.maf(maf = 'TCGA.ESCA.mutect.7f8e1e7c-621c-4dfd-8fad-af07c739dbfc.DR-10.0.somatic.maf')
getSampleSummary(rt)
getGeneSummary(rt)
getFields(rt)
getClinicalData(rt)
maf_df = rt@data
save(rt,maf_df,file = "maf_esca.Rdata")
length(unique(maf_df$Tumor_Sample_Barcode))
length(unique(maf_df$Hugo_Symbol))

plotmafSummary(maf = rt, 
               rmOutlier = TRUE,  
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE)

oncoplot(maf = rt, 
         top = 30,   
         fontSize = 0.6,  
         showTumorSampleBarcodes = F)   

somaticInteractions(maf = rt, 
                    top = 20,   
                    pvalue = c(0.05, 0.01, 0.001),   
                    showCounts = FALSE, 
                    fontSize = 0.5)

oncoplot(maf=rt, top=10, fontSize=0.5,
         annotationDat=clinical, clinicalFeature=c("IRGPI"),
         sortByAnnotation=TRUE, removeNonMutated=FALSE,groupAnnotationBySize=TRUE, drawColBar=T, keepGeneOrder=TRUE)
}
clinical<-rt@clinical.data
clinical$sample<-substr(clinical$Tumor_Sample_Barcode,1,15)
clinical<-merge(clinical,ESCC_IRGPI,by="sample")
clinical<-clusterProfiler::select(clinical,Tumor_Sample_Barcode,IRGPI,IRGPIScore)
clinical_low<-subset(clinical,IRGPI=='Low')
clinical_high<-subset(clinical,IRGPI=='High')
rt_low<-subsetMaf(rt,tsb =clinical_low$Tumor_Sample_Barcode)
rt_high<-subsetMaf(rt,tsb =clinical_high$Tumor_Sample_Barcode)
highlow<-mafCompare(m1 =rt_high, m2 = rt_low, m1Name = 'high_risk', m2Name = 'low_risk', minMut = 5)
pal <- pal_npg("nrc", alpha=0.8)(8)
pal
vc_cols <-pal
names(vc_cols) = c('In_Frame_Ins','Frame_Shift_Del','Missense_Mutation','Frame_Shift_Ins','Splice_Site',
                   'Multi_Hit','In_Frame_Del','Nonsense_Mutation')
print(vc_cols)
oncoplot(maf = rt_high, colors = vc_cols,titleText="ESCC IRGPI High Group(40 samples)",
         top = 10,   
         fontSize = 0.9,  
         showTumorSampleBarcodes = F)   

oncoplot(maf = rt_low, colors = vc_cols,titleText="ESCC IRGPI Low Group(41 samples)", 
         top = 10,   
         fontSize = 0.9,    
         showTumorSampleBarcodes = F)   
