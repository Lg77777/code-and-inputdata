
setwd("")
rm(list=ls())
# 加载必要的库
# 加载必要的库
library(ggplot2)
library(reshape2)
library(corrplot)
library(gplots)
library(PerformanceAnalytics)

# 读取数据
cibersort_results <- read.table('CIBERSORT-Results.txt', sep = '\t', header = TRUE, row.names = 1)
disulfidptosis_exp <- read.table('disulfidptosisExp.txt', sep = '\t', header = TRUE, row.names = 1)

# 确保两个数据集的样本名对齐
common_samples <- intersect(colnames(cibersort_results), colnames(disulfidptosis_exp))
cibersort_common <- cibersort_results[, common_samples]
disulfidptosis_common <- disulfidptosis_exp[, common_samples]

# 初始化相关性矩阵和P值矩阵
correlation_matrix <- matrix(NA, nrow(cibersort_common), nrow(disulfidptosis_common))
p_values <- matrix(NA, nrow(cibersort_common), nrow(disulfidptosis_common))
rownames(correlation_matrix) <- rownames(cibersort_common)
colnames(correlation_matrix) <- rownames(disulfidptosis_common)
rownames(p_values) <- rownames(correlation_matrix)
colnames(p_values) <- colnames(correlation_matrix)

# 计算每对基因的相关性和P值
for (i in 1:nrow(cibersort_common)) {
  for (j in 1:nrow(disulfidptosis_common)) {
    x <- as.numeric(cibersort_common[i, ])
    y <- as.numeric(disulfidptosis_common[j, ])
    if (all(!is.na(x)) & all(!is.na(y))) {
      cor_test <- cor.test(x, y, method = "pearson")
      correlation_matrix[i, j] <- cor_test$estimate
      p_values[i, j] <- cor_test$p.value
    }
  }
}

# 准备数据供corrplot使用
M <- correlation_matrix
res <- list(p = p_values)

# 表现形式10：多种颜色的热图
pdf(file="corplot10.pdf", width=8, height=8)
corrplot(M, method = "color", order = "original", type = "upper",col = colorRampPalette(c("blue", "white", "orange"))(50))
dev.off()
