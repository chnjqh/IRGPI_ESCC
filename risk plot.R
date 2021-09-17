library(rms)
library(ggrisk)
fit1 <- rms::cph(Surv(OS.time,OS)~CLDN1+HCAR3+FNBP1L+BRCA2,clinical_ESCC_expr)
ggrisk(fit1,
       cutoff.value='median',
       cutoff.x = 47,
       cutoff.y = -0.8,
       code.0 = 'Still Alive',
       code.1 = 'Already Dead',
       code.highrisk = 'High IRGPI',
       code.lowrisk = 'Low IRGPI',
       title.A.ylab='IRGPI Score',
       title.B.ylab='Survival Time(Days)',
       title.A.legend='IRGPI Group',
       title.B.legend='Status',
       title.C.legend='Expression',
       color.A=c(low='#00A087B2',high='#DC0000B2'),
       color.B=c(code.0='#00A087B2',code.1='#DC0000B2'), 
       color.C=c(low='#4DBBD5B2',median='white',high='#E64B35B2') 
)

fun_fit1<-Function(fit1)
sascode(fun_fit1)

candidate_genes_for_cox2 <- c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
candidate_genes_for_cox2
FinalGeneExp <- clinical_ESCC_expr[,candidate_genes_for_cox2]
FinalGeneExp$IRGPIScore<-fit1$linear.predictors
outCol = c("OS.time", "OS")
dat_cox = cbind(clinical_ESCC_expr[,outCol],FinalGeneExp)


dat_y<-dat_cox
dat_y$OS.time<-dat_y$OS.time/365
survivalROC_helper <- function(t) {
  
  survivalROC(Stime=dat_y$OS.time, status=dat_y$OS, marker = dat_y$IRGPIScore, 
              
              predict.time =t, method="KM")}
library(tidyverse)
library(survivalROC)
survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
survivalROC_data1 <- survivalROC_data %>% 
    mutate(auc =sprintf("%.3f",auc))%>% 
    unite(year, t,auc,sep = " Year AUC= ")
AUC =factor(survivalROC_data1$year)
mytheme<-theme_bw()+theme(plot.title = element_text(size = rel(2),hjust = 0.5),axis.title = element_text(size = rel(1.2)),axis.text = element_text(size=rel(1.2)),panel.grid.major = element_line(color = "white"),panel.grid.minor = element_line(colour = "white"),axis.text.x = element_text(size = rel(1.2), color = "black"),axis.text.y = element_text(size = rel(1.2),color = "black"),legend.text = element_text(size = rel(1.2)),legend.title = element_blank())
survivalROC_data1 %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_path(aes(color= AUC))+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  mytheme + theme(legend.position = c(0.8,0.2))+
  xlab("1-Specificity(FPR)")+ylab("Sensitivity(TPR)")+scale_color_lancet()
