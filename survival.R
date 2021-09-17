library(survminer)
library(survival)
surv_rnaseq.cut <- surv_cutpoint(
  clinical_ESCC_expr,
  time = "OS.time",
  event = "OS",
  variables = c("BRCA2"))
summary(surv_rnaseq.cut)
plot(surv_rnaseq.cut, "BRCA2", palette = "npg")

surv_rnaseq.cat <- surv_categorize(surv_rnaseq.cut)
fit_1 <- survfit(Surv(OS.time, OS) ~ BRCA2,data = surv_rnaseq.cat)
ggsurvplot(
  fit_1,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for point estimaes of survival curves.
  xlim = c(0,2000),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 500,      # break X axis in time intervals by 500.
  risk.table.y.text.col = T,  # colour risk table text annotations.
  risk.table.y.text = F,  #show bars instead of names in text annotations in legend of risk table
  xlab="Time(Days)", risk.table.col="strata",
  palette='npg')
