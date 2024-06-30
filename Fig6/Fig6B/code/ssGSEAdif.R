library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)

expFile="symbol.txt"        
gmtFile="immune.gmt"       
riskFile="risk.all.txt"     
setwd("")   


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
	

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}

data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
	

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(data, risk)


data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))


p=ggboxplot(data, x="Type", y="Score", fill = "Risk",
     notch=T, outlier.shape = NA,
     xlab="",ylab="Score",add = "none", palette=c("blue","red") )
p=p+rotate_x_text(50)
p1=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")


pdf(file="immFunction.pdf", width=7.5, height=5.5)
print(p1)
dev.off()



