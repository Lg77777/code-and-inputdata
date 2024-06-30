

library(limma)
corFilter=0.4          
pvalueFilter=0.001     
setwd("")    


rt=read.table("MAIAGs.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.2,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])      
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=read.table("disulfidptosisExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
disulfidptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
disulfidptosis=avereps(disulfidptosis)


group=sapply(strsplit(colnames(disulfidptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
disulfidptosis=disulfidptosis[,group==0]


outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.2){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(disulfidptosis)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(disulfidptosis[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Disulfidptosis=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(Disulfidptosis=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}


write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)


disulfidptosisLncRNA=unique(as.vector(outTab[,"lncRNA"]))
disulfidptosisLncRNAexp=data[disulfidptosisLncRNA,]
disulfidptosisLncRNAexp=rbind(ID=colnames(disulfidptosisLncRNAexp), disulfidptosisLncRNAexp)
write.table(disulfidptosisLncRNAexp,file="disulfidptosisLncExp.txt",sep="\t",quote=F,col.names=F)


