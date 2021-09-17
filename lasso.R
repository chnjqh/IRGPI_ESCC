library(glmnet)
require(doMC)
registerDoMC(cores=6)
set.seed(1234)  
x=as.matrix(clinical_ESCC_expr[,c(12:558)]) 
y=data.matrix(Surv(clinical_ESCC_expr$OS.time,clinical_ESCC_expr$OS)) 
fit<-glmnet(x, y, family = "cox", parallel=TRUE) 
plot(fit, xvar = "lambda", label = TRUE)
topptx(filename = "fit.pptx")
cvfit <- cv.glmnet(x, y, family="cox", parallel=TRUE)
plot(cvfit)
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
coef = coef(fit, s = cvfit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 

