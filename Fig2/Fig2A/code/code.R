library(limma)         
library(pheatmap)
expFile="mRNA.txt"       
geneFile="gene.txt"      
setwd("")    
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
gread.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]


write.table(outTab, file="disuinput", sep="\t", quote=F, col.names=F)



expFile="input.txt"
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))

sigVec=c()
outTab=data.frame()
for(i in rownames(data)){
  if(sd(data[i,])<0.001){next}
  wilcoxTest=wilcox.test(data[i,] ~ sampleType)
  pvalue=wilcoxTest$p.value
  if(pvalue<0.05){
    Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    sigVec=c(sigVec, paste0(i, Sig))
    conGeneMeans=mean(data[i,1:conNum])
    treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}

write.table(outTab, file="diff.xls", sep="\t", row.names=F, quote=F)
write.table(outTab, file="diff.txt", sep="\t", row.names=F, quote=F)


exp=data[as.vector(outTab[,1]),]
expOut=rbind(ID=colnames(exp), exp)
write.table(expOut, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)


exp=log2(exp+0.1)
row.names(exp)=sigVec
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=9, height=6)
pheatmap(exp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()
