library(limma)
library(ggpubr)
tideFile="TIDE.csv"        
riskFile="risk.all.txt"     
setwd("")     


tide=read.csv(tideFile, header=T, sep=",", check.names=F, row.names=1)
tide=tide[,"TIDE",drop=F]     


group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	

sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)
	

data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
	         xlab="", ylab="TIDE",
	         palette=c("#0066FF","#FF0000"),
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")


pdf(file="TIDE.pdf", width=6, height=5)
print(gg1)
dev.off()



