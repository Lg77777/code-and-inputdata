
library(survival)
library(regplot)
library(rms)
library(survcomp)

riskFile="risk.all.txt"      
cliFile="clinical.txt"       
setwd("")     


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)


res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),     
             clickable=F,
             title="",             
             points=TRUE,           
             droplines=TRUE,        
             observation=rt[1,],    
             rank="sd",
             failtime = c(1,3,5),   
             prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")


nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


pdf(file="calibration.pdf", width=5, height=5)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)

f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('topleft', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')

cindex=concordance.index(x=nomoRisk, surv.time=rt$futime, surv.event=rt$fustat, method= "noether")
c_index=sprintf("%.03f", cindex$c.index)
c_index.ci_low=sprintf("%.03f", cindex$lower)
c_index.ci_high=sprintf("%.03f", cindex$upper)
cindexLabel=paste0(c_index, "(95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")
text(0.5, 0.1, "C-index:")
text(0.7, 0.03, cindexLabel)
dev.off()


