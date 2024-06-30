library(limma)
library(ggpubr)
setwd("")      


tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	

risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
	

sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)
	

data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	

boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
			      xlab="",
			      ylab="Tumor tmbation burden (log2)",
			      legend.title="",
			      palette = c("#0066FF","#FF0000"),
			      add = "boxplot", add.params = list(fill="white"))+ 
	stat_compare_means(comparisons = my_comparisons)
	

pdf(file="riskTMB.pdf", width=5, height=4.5)
print(boxplot)
dev.off()


