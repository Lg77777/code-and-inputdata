
library(survival)
library(survminer)

riskFile="risk.all.txt"     
tmbFile="TMB.txt"            
setwd("")     


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      


sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)


res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)


bioSurvival=function(surData=null, outFile=null){
  
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	
	#????????????
	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=F,
			           cumevents=F,
			           risk.table.height=.25)

	pdf(file=outFile, width=5.5, height=4.8, onefile=FALSE)
	print(surPlot)
	dev.off()
}


data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

a$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")


####