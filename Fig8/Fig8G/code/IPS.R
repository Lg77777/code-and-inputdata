library(reshape2)
library(ggpubr)

# Set working directory (adjust the path as necessary)
setwd("")

# Read the input files
risk_data <- read.table("risk.all.txt", header = TRUE, sep = "\t", check.names = FALSE)
ips_data <- read.table("TCIA-ClinicalData.LUAD.tsv", header = TRUE, sep = "\t", check.names = FALSE)

# Ensure the IPS data has a common identifier column named 'id' for merging
ips_data$id <- ips_data$barcode

# Merge data by common identifier (assume column "id" in both data frames)
merged_data <- merge(risk_data, ips_data, by = "id")

# Initialize a list to store p-values and plots
p_values <- list()
plot_list <- list()

# Loop through each column in ips_data (excluding the first column, which is 'id')
for (column in colnames(ips_data)[-1]) {
  if (is.numeric(merged_data[[column]])) {
    # Perform Wilcoxon test only if the column is numeric
    test_result <- wilcox.test(merged_data[[column]] ~ merged_data$risk)
    p_values[[column]] <- test_result$p.value
    
    # Create boxplot
    p <- ggboxplot(
      merged_data, x = "risk", y = column, 
      color = "risk", palette = "jco",
      add = "none", ylab = column, xlab = "Risk"
    ) + stat_summary(fun = median, geom = "point", size = 2, color = "black") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        legend.position = "none"
      ) + stat_compare_means(method = "wilcox.test")
    
    # Add significance annotation
    p_value <- test_result$p.value
    if (p_value < 0.001) {
      signif <- "***"
    } else if (p_value < 0.01) {
      signif <- "**"
    } else if (p_value < 0.05) {
      signif <- "*"
    } else {
      signif <- "ns"
    }
    
    max_y <- max(merged_data[[column]], na.rm = TRUE)
    p <- p + annotate("text", x = 1.5, y = max_y, label = signif, color = "red", size = 5, vjust = -0.5)
    
    plot_list[[column]] <- p
  } else {
    p_values[[column]] <- NA
  }
}

# Arrange plots in one row
ggarrange(plotlist = plot_list, ncol = 4, nrow = 1)
