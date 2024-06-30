
library(survival)
library(survminer)
setwd("")    


bioSurvival=function(inputFile=null, outFile=null){
	
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by = 2,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	
	
	pdf(file=outFile, width=6, height=5, onefile=FALSE)
	print(surPlot)
	dev.off()
}


bioSurvival(inputFile="risk.all.txt", outFile="surv.all.pdf")


