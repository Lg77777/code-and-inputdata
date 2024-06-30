library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
setwd("") # Set working directory
# Read input file
rt = read.table("expTime.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt$futime[rt$futime <= 0] = 1
rt$futime = rt$futime / 365
rt[, 3:ncol(rt)] = log2(rt[, 3:ncol(rt)] + 1)

# Define forest plot function
bioForest = function(coxFile = NULL, forestFile = NULL, forestCol = NULL) {
  rt <- read.table(coxFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  gene <- rownames(rt)
  hr <- sprintf("%.3f", rt$"HR")
  hrLow <- sprintf("%.3f", rt$"HR.95L")
  hrHigh <- sprintf("%.3f", rt$"HR.95H")
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  # Output plot
  pdf(file = forestFile, width = 7, height = 6)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
  
  xlim = c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  text.cex = 0.8
  text(0, n:1, gene, adj = 0, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, adj = 1)
  text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(3, n + 1, 'Hazard ratio', cex = text.cex, adj = 1)
  
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  LOGindex = 10 
  hrLow = log(as.numeric(hrLow), LOGindex)
  hrHigh = log(as.numeric(hrHigh), LOGindex)
  hr = log(as.numeric(hr), LOGindex)
  xlim = c(floor(min(hrLow, hrHigh)), ceiling(max(hrLow, hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
  abline(v = log(1, LOGindex), col = "black", lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > log(1, LOGindex), forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3)
  a1 = axis(1, labels = FALSE, tick = FALSE)
  axis(1, a1, 10^a1)
  dev.off()
}

# Build prognosis model
n = 1 # Number of random groups
coxPfilter = 0.05 # Significance threshold for Cox regression

for (i in 1:n) {
  # Split data into train and test sets
  inTrain <- createDataPartition(y = rt[, 2], p = 0.5, list = FALSE)
  train <- rt[inTrain, ]
  test <- rt[-inTrain, ]
  trainOut = cbind(id = row.names(train), train)
  testOut = cbind(id = row.names(test), test)
  
  # Univariate Cox regression analysis
  outUniTab = data.frame()
  sigGenes = c("futime", "fustat")
  for (i in colnames(train[, 3:ncol(train)])) {
    if (sd(train[, i]) > 0.1) {
      cox <- coxph(Surv(futime, fustat) ~ train[, i], data = train)
      coxSummary = summary(cox)
      coxP = coxSummary$coefficients[, "Pr(>|z|)"]
      
      # Retain significant genes
      if (coxP < coxPfilter) {
        sigGenes = c(sigGenes, i)
        outUniTab = rbind(outUniTab,
                          cbind(id = i,
                                HR = coxSummary$conf.int[, "exp(coef)"],
                                HR.95L = coxSummary$conf.int[, "lower .95"],
                                HR.95H = coxSummary$conf.int[, "upper .95"],
                                pvalue = coxSummary$coefficients[, "Pr(>|z|)"])
        )
      }
    }	
  }
  uniSigExp = train[, sigGenes]
  uniSigExpOut = cbind(id = row.names(uniSigExp), uniSigExp)
  if (ncol(uniSigExp) < 6) { next }
  
  # LASSO regression
  x = as.matrix(uniSigExp[, 3:ncol(uniSigExp)])
  y = data.matrix(Surv(uniSigExp$futime, uniSigExp$fustat))
  fit = glmnet(x, y, family = "cox", maxit = 1000)
  cvfit = cv.glmnet(x, y, family = "cox", maxit = 1000)
  coef = coef(fit, s = cvfit$lambda.min)
  index = which(coef != 0)
  actCoef = coef[index]
  lassoGene = row.names(coef)[index]
  lassoSigExp = uniSigExp[, c("futime", "fustat", lassoGene)]
  lassoSigExpOut = cbind(id = row.names(lassoSigExp), lassoSigExp)
  geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
  if (nrow(geneCoef) < 2) { next }
  
  # Build Cox model
  multiCox = coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
  multiCox = step(multiCox, direction = "both")
  multiCoxSum = summary(multiCox)
  
  # Output model formula
  outMultiTab = data.frame()
  outMultiTab = cbind(
    coef = multiCoxSum$coefficients[, "coef"],
    HR = multiCoxSum$conf.int[, "exp(coef)"],
    HR.95L = multiCoxSum$conf.int[, "lower .95"],
    HR.95H = multiCoxSum$conf.int[, "upper .95"],
    pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"]
  )
  outMultiTab = cbind(id = row.names(outMultiTab), outMultiTab)
  outMultiTab = outMultiTab[, 1:2]
  
  # Divide train set into high and low risk groups based on median risk score
  riskScore = predict(multiCox, type = "risk", newdata = train)
  coxGene = rownames(multiCoxSum$coefficients)
  coxGene = gsub("`", "", coxGene)
  outCol = c("futime", "fustat", coxGene)
  medianTrainRisk = median(riskScore)
  risk = as.vector(ifelse(riskScore > medianTrainRisk, "high", "low"))
  trainRiskOut = cbind(id = rownames(cbind(train[, outCol], riskScore, risk)), cbind(train[, outCol], riskScore, risk))
  
  # Divide test set into high and low risk groups based on median risk score
  riskScoreTest = predict(multiCox, type = "risk", newdata = test)
  riskTest = as.vector(ifelse(riskScoreTest > medianTrainRisk, "high", "low"))
  testRiskOut = cbind(id = rownames(cbind(test[, outCol], riskScoreTest, riskTest)), cbind(test[, outCol], riskScore = riskScoreTest, risk = riskTest))
  
  # Compare survival differences between high and low risk groups, obtain P value
  diff = survdiff(Surv(futime, fustat) ~ risk, data = train)
  pValue = 1 - pchisq(diff$chisq, df = 1)
  diffTest = survdiff(Surv(futime, fustat) ~ riskTest, data = test)
  pValueTest = 1 - pchisq(diffTest$chisq, df = 1)
  
  write.table(trainOut, file = "data.train.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(testOut, file = "data.test.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Output univariate Cox regression results
  write.table(outUniTab, file = "uni.trainCox.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(uniSigExpOut, file = "uni.SigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  bioForest(coxFile = "uni.trainCox.txt", forestFile = "uni.foreast.pdf", forestCol = c("red", "green"))
  
  # Output LASSO regression results
  write.table(lassoSigExpOut, file = "lasso.SigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  pdf("lasso.lambda.pdf")
  plot(fit, xvar = "lambda", label = TRUE)
  dev.off()
  pdf("lasso.cvfit.pdf")
  plot(cvfit)
  abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
  dev.off()
  
  # Output multivariate Cox regression results
  write.table(outMultiTab, file = "multiCox.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(trainRiskOut, file = "risk.train.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(testRiskOut, file = "risk.test.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  allRiskOut = rbind(trainRiskOut, testRiskOut)
  write.table(allRiskOut, file = "risk.all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  break
}

