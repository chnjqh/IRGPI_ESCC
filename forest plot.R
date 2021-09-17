
Unicox<- function(x){
  form1 <- as.formula(paste0("Surv(OS.time,OS)~",x))
  cox_form1 <- coxph(form1,data = clinical_ESCC_expr)
  cox_Sum <- summary(cox_form1)
  HR=cox_Sum$conf.int[,"exp(coef)"]
  HR.95L=cox_Sum$conf.int[,"lower .95"]
  HR.95H=cox_Sum$conf.int[,"upper .95"]
  pvalue=cox_Sum$coefficients[,"Pr(>|z|)"]
  Unicox <- data.frame("Characteristics" = x,
                       "HR" = HR,
                       "HR.95L" = HR.95L,
                       "HR.95H" = HR.95H,
                       "pvalue" = pvalue)
  return(Unicox)
}
library(plyr)
Univar_gene <- lapply(candidate_genes_for_cox2, Unicox)
Univar_gene <- ldply(Univar_gene,data.frame)
Univar_gene[,2:ncol(Univar_gene)] <- as.numeric(unlist(Univar_gene[,2:ncol(Univar_gene)]))
hz <- paste(round(Univar_gene$HR,3),
            "(",round(Univar_gene$HR.95L,3),
            "-",round(Univar_gene$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"Gene",Univar_gene$Characteristics),
                   c(NA,"P value",round(Univar_gene$pvalue,3)),
                   c(NA,"Hazard Ratio(95% CI)",hz))
library(forestplot)
forestplot(labeltext=tabletext,
           graph.pos=2,  
           col=fpColors(box="#4DBBD5B2", lines="#3C5488B2", zero = "gray50"),
           mean=c(NA,NA,Univar_gene$HR),
           lower=c(NA,NA,Univar_gene$HR.95L), 
           upper=c(NA,NA,Univar_gene$HR.95H), 
           boxsize=0.3,lwd.ci=2,  
           ci.vertices.height = 0.1,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(20,"mm"),    
           xticks = c(0.5, 1,1.5), 
           lwd.xaxis=1,            
           lineheight = unit(2,"cm"), 
           graphwidth = unit(.3,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawDiamondCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), 
                           "7" = gpar(lwd=2, col="black")),
           mar=unit(rep(0.5, times = 4), "cm"),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")