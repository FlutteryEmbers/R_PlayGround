setwd(".")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
library("pheatmap")
list.files()
sample_metadata = read.table('./scripts/SRP187485_summary.csv', header = TRUE, sep = ",")
sample_metadata <- sample_metadata[,c(1,2,5,6)]

sample_metadata$sampleID <- factor(sample_metadata$sampleID)
sample_metadata$runID <- factor(sample_metadata$runID)
sample_metadata$treatment <- factor(sample_metadata$treatment)
sample_metadata$time <- factor(sample_metadata$time)


print(sample_metadata$sampleID)
print(sample_metadata$runID)
print(sample_metadata$treatment)
print(sample_metadata$time)

# idx = which(sample_metadata$treatment == "Intact")
# sample_metadata[idx,]
# summary(sample_metadata$treatment)