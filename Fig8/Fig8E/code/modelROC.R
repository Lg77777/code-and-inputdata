library(limma)
library(survival)
library(survminer)
library(timeROC)

expFile="symbol.txt"       
geneFile="TISgene.txt"       
riskFile="risk.TCGA.txt"    
tideFile="TIDE.txt"         
setwd("")     


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

eRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
geneExp=data[sameGene,]
logData=log2(geneExp+1)
TISscore=colMeans(logData)

#??ÈEscore=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(TIDEscore), names(TISscore))
score=cbind(TIDEscore[sameSample,,drop=F], TIS=TISscore[sameSample])
score=score[,c("TIDE","TIS")]

#È¥?up=sapply(strsplit(row.names(score),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
score=score[group==0,]
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)

#??Èk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(risk))
score=score[sameSample, , drop=F]
risk=risk[sameSample, c("futime","fustat","riskScore"), drop=F]
rt=cbind(risk, score)

####Col=rainbow(ncol(rt)-2, s=0.9, v=0.9)
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
	           marker=risk$riskScore,cause=1,
	           weighting='aalen',
	           times=c(1,2,3),ROC=TRUE)
pdf(file="ROC.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()

####dictTime=3     #???Text=c()
pdf(file="multiROC.pdf", width=5.5, height=5.5)
#???
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#??T(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#???end("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()


####