
library(limma)
library(survival)
library(survminer)
immuneFile="CIBERSORT-Results.txt"    
surFile="time.txt"                    
pFilter=0.05                           
setwd("")     


immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])


group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)


surTime=read.table(surFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
rt=cbind(surTime, data)
rt$futime=rt$futime/365


outTab=data.frame()
for(immune in colnames(rt)[3:ncol(rt)]){
	if(sd(rt[,immune])<0.01){next}
	data=rt[,c("futime", "fustat", immune)]
	colnames(data)=c("futime", "fustat", "immune")
	

	res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("immune"))
	res.cat=surv_categorize(res.cut)


	diff=survdiff(Surv(futime, fustat) ~ immune, data =res.cat)
	pValue=1-pchisq(diff$chisq, df=1)
	if(pValue<0.05){
		outVector=cbind(immune, res.cut$cutpoint[1], pValue)
		outTab=rbind(outTab,outVector)
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		
		
		fit=survfit(Surv(futime, fustat) ~ immune, data = res.cat)
		surPlot=ggsurvplot(fit,
						   data=res.cat,
						   pval=pValue,
						   pval.size=6,
						   legend.title=immune,
						   legend.labs=c("high","low"),
						   xlab="Time(years)",
						   break.time.by=1,
						   palette=c("red", "blue"),
						   conf.int=F,
						   surv.median.line="hv",
						   risk.table=F,
						   risk.table.title="",
						   risk.table.height=.25)
		
		pdf(file=paste0("Sur.", immune, ".pdf"), onefile = FALSE, width = 6, height = 5)
		print(surPlot)
		dev.off()
	}
}

write.table(outTab,file="survival.result.txt",sep="\t",row.names=F,quote=F)




