library(limma)
library(ggpubr)
File1="riskscore.txt"        
File2="group.txt"     
setwd("C:\\Users\\Administrator\\Desktop\\新建文件夹 (4)\\Fig8\\Fig8J\\code")   


a=read.table(File1, header=T, sep="\t", check.names=F, row.names=1)
a=a[,"Risk score",drop=F]    


b=read.table(File2, header=T, sep="\t", check.names=F, row.names=1)
	
sameSample=intersect(row.names(a), row.names(b))
a=a[sameSample, , drop=F]
b=b[sameSample, "Responser", drop=F]
data=cbind(a, b)
	
data$Responser=ifelse(data$Responser=="high", "Responser", "non-Responser")
group=levels(factor(data$Responser))
data$Responser=factor(data$Responser, levels=c("non-Responser", "Responser"))
group=levels(factor(data$Responser))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


gg1=ggviolin(data, x="Responser", y="Risk score", fill = "Responser", 
	         xlab="Responser", ylab="Risk score",
	         palette=c("#0066FF","#FF0000"),
	         legend.title="Risk score",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

	
pdf(file="Riskscore.pdf", width=6, height=5)
print(gg1)
dev.off()




