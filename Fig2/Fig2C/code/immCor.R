######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

riskFile="risk.all.txt"             #风险文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
setwd("C:\\biowolf\\DRG\\33.immCor")      #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"risk",drop=F])
rt=rt[order(rt$risk, decreasing=T),]
conNum=nrow(rt[rt$risk=="low",])
treatNum=nrow(rt[rt$risk=="high",])

##########绘制柱状图##########
data=t(rt[,-ncol(rt)])
pdf("barplot.pdf", width=16, height=9)
#定义图形的颜色
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=rt
data=melt(data, id.vars=c("risk"))
colnames(data)=c("Risk", "Immune", "Expression")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))

#定义箱线图的颜色
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

#绘制箱线图
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Risk",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图形
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

