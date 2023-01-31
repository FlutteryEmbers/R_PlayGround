.libPaths()
# Add if not present "/scratch/user/r.kapoor/my_R_libs" 
#.libPaths(c(.libPaths(),"/scratch/user/r.kapoor/my_R_libs/"))
setwd("/scratch/user/r.kapoor/tiller_bud")


gene_counts <- read.table("SRP187485_genes.featurecounts.tsv",header=TRUE,row.names = 1)
sample_metadata <- read.table("SRP187485_summary.csv",header=TRUE,sep=',')
# For NOT including tech replicate, ignore o/w
sample_metadata <- sample_metadata[,c(1,2,5,6)]
# Both cases
sample_metadata$sampleID <- factor(sample_metadata$sampleID)
sample_metadata$runID <- factor(sample_metadata$runID)
sample_metadata$treatment <- factor(sample_metadata$treatment)
sample_metadata$time <- factor(sample_metadata$time)


print(sample_metadata$sampleID)
print(sample_metadata$runID)
print(sample_metadata$treatment)
print(sample_metadata$time)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = sample_metadata,
                              design = ~ 1 + treatment + time + treatment:time)


ddsColl <- collapseReplicates(dds, sample_metadata$sampleID, sample_metadata$runID)

keep <- rowSums(counts(ddsColl) > 10) >= 6
ddsColl <- ddsColl[keep,]

library("BiocParallel")
register(MulticoreParam(8))

ddsColl <- DESeq(ddsColl,betaPrior=FALSE, parallel=TRUE)
resultsNames(ddsColl)

par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsColl)[["cooks"]]), range=0, las=2)

vsd <- vst(ddsColl)
plotPCA(vsd, intgroup=c("treatment", "time"))

rlog_data<- rlog(ddsColl, blind=FALSE)
saveRDS(rlog_data, file = "rlog_data_coll.rds")

plotPCA(rlog_data, intgroup=c("treatment", "time"))


#############
gene_counts <- read.table("SRP187485_genes.featurecounts.tsv",header=TRUE,row.names = 1)
sample_metadata <- read.table("SRP187485_summary.csv",header=TRUE,sep=',')
# For including tech replicate
sample_metadata <- sample_metadata[,c(1,2,4,5,6)]
# Both cases
sample_metadata$sampleID <- factor(sample_metadata$sampleID)
sample_metadata$runID <- factor(sample_metadata$runID)
sample_metadata$treatment <- factor(sample_metadata$treatment)
sample_metadata$time <- factor(sample_metadata$time)
# For including tech replicate
sample_metadata$techReplicate <- factor(sample_metadata$techReplicate)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = sample_metadata,
                              design = ~ 1 + treatment + time + treatment:time)

keep <- rowSums(counts(dds) > 10) >= 6
dds <- dds[keep,]

library("BiocParallel")
register(MulticoreParam(8))

dds <- DESeq(dds,betaPrior=FALSE, parallel=TRUE)
resultsNames(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

vsd <- vst(dds)
plotPCA(vsd, intgroup=c("treatment", "time"))

rlog_data<- rlog(dds, blind=FALSE)
saveRDS(rlog_data, file = "rlog_data.rds")

plotPCA(rlog_data, intgroup=c("treatment", "time"))

