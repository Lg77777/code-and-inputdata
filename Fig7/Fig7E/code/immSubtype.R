library(RColorBrewer)       
riskFile="risk.TCGA.txt"                         
subtypeFile="Subtype_Immune_Model_Based.txt"      
setwd("")    


ChischisqTest <- function(tabledata){
	ct = chisq.test(tabledata)
	ct.pvalue = ct$p.value
	ct.pvalue = ifelse(ct.pvalue<0.001,0.001,round(ct.pvalue,3))
}


rect2text <- function(x1,y1,x2,y2,text,rect.col,text.cex,text.col="white",...){
  rect(x1,y1,x2,y2,col=rect.col,border=NA,...)
  text((x1+x2)/2,(y1+y2)/2,text,cex=text.cex,col=text.col,...)
}


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))


sameSample=intersect(row.names(subtype), row.names(risk))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
risk=risk[sameSample,]
data=cbind(subtype, as.data.frame(risk))
typeTab=table(data$subtype)
typeName=names(typeTab[typeTab>6])
merge.data=data[which(data[,"subtype"] %in% typeName),]
cliName=colnames(merge.data)[1]
colnames(merge.data)[1]="Clinical"


groups = sort(unique(merge.data$Clinical))
groups.table = table(merge.data$Clinical)
groups.num = length(groups)
samples.num = nrow(merge.data)
groups.lownum = sum(merge.data$Risk=="low")
groups.highnum = sum(merge.data$Risk=="high")
pergrouptable = table(merge.data$Clinical, merge.data$Risk)
ct.pvalue = ChischisqTest(pergrouptable)


xlim1=0; xlim2=samples.num; ylim1=1; ylim2=10
space = dist.up = 0.1
space.x = 1/100*xlim2
groups.col = col = colorRampPalette(brewer.pal(9, "Paired"))(groups.num)
high_low.col = c(rgb(235/255,116/255,106/255),rgb(30/255,181/255,184/255))
bottom.unitwidth = xlim2 / (groups.num+1+1.3)
bottom.left.width = bottom.unitwidth*1.3


pdf(file="immSubtype.pdf", width=10, height=7)
plot(1,type="n",xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),axes=F,xlab="",ylab="")


header.text.cex = 1.5
header.text = gettextf("%s TCGA patients",samples.num)
header.text.width = strwidth(header.text,cex=header.text.cex)
arrows.len = (xlim2 - header.text.width - space.x*2)/2
text(xlim2/2,9.5,header.text,cex=header.text.cex,font=2)
arrows(0,9.5,arrows.len,9.5,code=1,lwd=3,length=0.2)
arrows(xlim2-arrows.len,9.5,xlim2,9.5,code=2,lwd=3,length=0.2)
rect2text(x1=0,y1=4.5+space,x2=bottom.left.width,y2=7-space,text.cex=1.2,text = 
	"IRGPI\n groups",rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)

bottom.header.text = paste0(cliName, " group (n=" , samples.num, ")")
rect2text(x1=bottom.left.width+space.x,y1=6+space/2,x2=xlim2,y2=7-space,text.cex=1.2,
	text= bottom.header.text,rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)


bottom.left2.text = gettextf("IRGPI-low\n(n=%s)",groups.lownum)
rect2text(x1=0,y1=3+space/2,x2=bottom.left.width,y2=4.5-space/2,text.cex=1.2,text= 
	bottom.left2.text,rect.col=high_low.col[2],text.col="white",font=2)

bottom.left3.text = gettextf("IRGPI-high\n(n=%s)",groups.highnum)
rect2text(x1=0,y1=1.5+space/2,x2=bottom.left.width,y2=3-space/2,text.cex=1.2,text= 
	bottom.left3.text,rect.col=high_low.col[1],text.col="white",font=2)


start = 0
for(i in 1:length(groups.table)){
	end = groups.table[i]+start
	namesi = names(groups.table)[i]
	
	rect2text(x1=start,y1=8+space,x2=end,y2=9-space,text.cex=1.1,text= namesi, 
	  	rect.col=groups.col[i],text.col="white")
	merge.datai = merge.data[merge.data$Clinical==namesi,,drop=F]
	high.num = sum(merge.datai$Risk=="high")
	low.num = sum(merge.datai$Risk=="low")

	rect(start,7+dist.up,start+low.num,8-dist.up,col=high_low.col[2],border=NA)
	
	rect(start+low.num,7+dist.up,end,8-dist.up,col=high_low.col[1],border=NA)
	
	bottom.starti = bottom.left.width+(i-1)* bottom.unitwidth
	bottom.endi = bottom.starti+bottom.unitwidth
	subheader.text = gettextf("%s\n(n=%s, %s)",namesi,nrow(merge.datai),paste0(round(nrow(merge.datai)/samples.num*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=4.5+space,x2=bottom.endi,y2=6-space,text.cex=1.1,text= subheader.text,rect.col=groups.col[i],text.col="white")

	low.texti = gettextf("%s(%s)",low.num,paste0(round(low.num/groups.lownum*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=3+space/2,x2=bottom.endi,y2=4.5-space/2,text.cex=1.1,text= low.texti,rect.col="grey90",text.col="black")
	
	high.texti = gettextf("%s(%s)",high.num,paste0(round(high.num/groups.highnum*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=1.5+space/2,x2=bottom.endi,y2=3-space/2,text.cex=1.1,text= high.texti,rect.col="grey70",text.col="black")
	start = end
}


rect2text(x1=bottom.endi+space.x,y1=4.5+space,x2=xlim2,y2=6-space,text.cex=1.1,text= "P-value",rect.col="grey70",text.col="black",font=3)
rect2text(x1=bottom.endi+space.x,y1=1.5+space/2,x2=xlim2,y2=4.5-space/2,text.cex=1.1,text= ct.pvalue,rect.col="grey70",text.col="black")


rect(0,8+space,xlim2,9-space,lwd=2,border="grey30")
rect(0,7+space,xlim2,8-space,lwd=2,border="grey30")
rect(0,1.5+space/2,xlim2,7-space,lwd=2,border="grey30")


legend("bottom",legend=c("IRGPI-low","IRGPI-high"),col=rev(high_low.col),bty="n",ncol=2,pch=15,pt.cex=1.3,cex=1.3)
dev.off()


