setwd("/scratch/user/r.kapoor/TAM111_112/second_try/scripts/")

load("JGL/transformed_data_new.RData")

setwd("/scratch/user/r.kapoor/wheat/wgcna_grain/grain_RNAseq")
library(WGCNA,"/mnt/scratch/user/r.kapoor/R/library-4.0.3/00LOCK-WGCNA/00new/WGCNA/libs")
options(stringsAsFactors = FALSE);



net.out <- readRDS("HCgrainNoTOM.9.signed hybrid.it.500.rds.parts.rds")

common_genes <- names(net.out[[71]])
datExpr <- t(transformed_dds_data[rownames(transformed_dds_data) %in% common_genes,])


modules_upTAM111_A <- c(33,61,93,18,40,13,81,64,75,56)
modules_upTAM111_B <- c(6,7,16,43,82,60)
modules_upTAM111_C <- c(27,24,42,62,63,96,58,12,50,71,2,10,70)

sum((net.out[[71]] %in% modules_upTAM111_A))
sum((net.out[[71]] %in% modules_upTAM111_B))
sum((net.out[[71]] %in% modules_upTAM111_C))



setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts")

allTFnames <- read.table("TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allcofactors <- read.table("cofactors_v1.1.csv",sep=',',header=FALSE, stringsAsFactors = FALSE)

# Group A

modules_upTAM111_genes <- names(net.out[[71]][net.out[[71]] %in% modules_upTAM111_A])
modules_upTAM111_genes <- modules_upTAM111_genes[modules_upTAM111_genes %in% colnames(datExpr)]

all_genes_JGL <- unique(c(allTFnames$gene,allcofactors$V1,modules_upTAM111_genes))
all_genes_JGL <- all_genes_JGL[all_genes_JGL %in% colnames(datExpr)]

datExpr_TAM111up <- datExpr[,all_genes_JGL]


## JGL 

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

#coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111",]
#coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112",]

#coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111" & TAM_metadata$day=="14",]
#coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112" & TAM_metadata$day=="14",]

coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111" & TAM_metadata$day=="5",]
coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112" & TAM_metadata$day=="5",]

gene_data_TAM111.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM111),])
gene_data_TAM112.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM112),])

#library(JGL,lib.loc = "/scratch/user/r.kapoor/my_R_libs/")

sen_JGL_in <- list(gene_data_TAM111.npn,gene_data_TAM112.npn)
fgl.screen <- screen.fgl(Y=sen_JGL_in, lambda1=1, lambda2=0, weights = "equal")
sum(fgl.screen)

#12353

#sum(all_genes_JGL[fgl.screen] %in% allTFnames$gene)
#1748
sum(all_genes_JGL[fgl.screen] %in% allcofactors$V1)

sum(all_genes_JGL[fgl.screen] %in% unique(allTFnames$gene))

screened_gene_data_TAM111.npn<- (sen_JGL_in[[1]])[,fgl.screen]
screened_gene_data_TAM112.npn <- (sen_JGL_in[[2]])[,fgl.screen]
screened_sen_JGL_in <- list(screened_gene_data_TAM111.npn,screened_gene_data_TAM112.npn)

#saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_A.rds")
# saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_A_day14.rds")
# screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_A_day14.rds")
# fgl.results.TAM111up.A.day14 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
#                             lambda2=0,return.whole.theta = TRUE)
# 
# saveRDS(fgl.results.TAM111up.A.day14, "fgl.results.TAM111up.A.day14.rds")
# print(fgl.results.TAM111up.A.day14)

saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_A_day5.rds")
screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_A_day5.rds")
fgl.results.TAM111up.A.day5 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
                                    lambda2=0,return.whole.theta = TRUE)

saveRDS(fgl.results.TAM111up.A.day5, "fgl.results.TAM111up.A.day5.rds")
print(fgl.results.TAM111up.A.day5)

# Group B
rm(list=ls())
modules_upTAM111_genes <- names(net.out[[71]][net.out[[71]] %in% modules_upTAM111_B])
modules_upTAM111_genes <- modules_upTAM111_genes[modules_upTAM111_genes %in% colnames(datExpr)]

all_genes_JGL <- unique(c(allTFnames$gene,allcofactors$V1,modules_upTAM111_genes))
all_genes_JGL <- all_genes_JGL[all_genes_JGL %in% colnames(datExpr)]

datExpr_TAM111up <- datExpr[,all_genes_JGL]


## JGL 

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

#coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111",]
#coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112",]

coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111" & TAM_metadata$day=="14",]
coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112" & TAM_metadata$day=="14" & TAM_metadata$rep !="2",]

gene_data_TAM111.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM111),])
gene_data_TAM112.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM112),])

#library(JGL,lib.loc = "/scratch/user/r.kapoor/my_R_libs/")

sen_JGL_in <- list(gene_data_TAM111.npn,gene_data_TAM112.npn)
fgl.screen <- screen.fgl(Y=sen_JGL_in, lambda1=1, lambda2=0, weights = "equal")
sum(fgl.screen)

#12353

#sum(all_genes_JGL[fgl.screen] %in% allTFnames$gene)
#1748
sum(all_genes_JGL[fgl.screen] %in% allcofactors$V1)

sum(all_genes_JGL[fgl.screen] %in% unique(allTFnames$gene))

screened_gene_data_TAM111.npn<- (sen_JGL_in[[1]])[,fgl.screen]
screened_gene_data_TAM112.npn <- (sen_JGL_in[[2]])[,fgl.screen]
screened_sen_JGL_in <- list(screened_gene_data_TAM111.npn,screened_gene_data_TAM112.npn)

#saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_B.rds")

# saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_B_day14.rds")
# screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_B_day14.rds")
# fgl.results.TAM111up.B.day14 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
#                                   lambda2=0,return.whole.theta = TRUE)
# 
# saveRDS(fgl.results.TAM111up.B.day14, "fgl.results.TAM111up.B.day14.rds")
# print(fgl.results.TAM111up.day14)

# saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_B_day14.rem9.rds")
# screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_B_day14.rem9.rds")
# fgl.results.TAM111up.B.day14.rem9 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
#                                     lambda2=0,return.whole.theta = TRUE)
# 
# saveRDS(fgl.results.TAM111up.B.day14.rem9, "fgl.results.TAM111up.B.day14.rem9.rds")
# print(fgl.results.TAM111up.B.day14.rem9)

# Group C
rm(list=ls())
modules_upTAM111_genes <- names(net.out[[71]][net.out[[71]] %in% modules_upTAM111_C])
modules_upTAM111_genes <- modules_upTAM111_genes[modules_upTAM111_genes %in% colnames(datExpr)]

all_genes_JGL <- unique(c(allTFnames$gene,allcofactors$V1,modules_upTAM111_genes))
all_genes_JGL <- all_genes_JGL[all_genes_JGL %in% colnames(datExpr)]

datExpr_TAM111up <- datExpr[,all_genes_JGL]


## JGL 

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

#coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111",]
#coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112",]

#coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111" & TAM_metadata$day=="14",]
#coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112" & TAM_metadata$day=="14",]

coldata_TAM111 <- TAM_metadata[TAM_metadata$geno=="TAM111" & TAM_metadata$day=="14",]
coldata_TAM112 <- TAM_metadata[TAM_metadata$geno=="TAM112" & TAM_metadata$day=="14" & TAM_metadata$rep !="2",]

gene_data_TAM111.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM111),])
gene_data_TAM112.npn<-as.matrix(datExpr_TAM111up[rownames(coldata_TAM112),])

#library(JGL,lib.loc = "/scratch/user/r.kapoor/my_R_libs/")

sen_JGL_in <- list(gene_data_TAM111.npn,gene_data_TAM112.npn)
fgl.screen <- screen.fgl(Y=sen_JGL_in, lambda1=1, lambda2=0, weights = "equal")
sum(fgl.screen)

#12353

#sum(all_genes_JGL[fgl.screen] %in% allTFnames$gene)
#1748
sum(all_genes_JGL[fgl.screen] %in% allcofactors$V1)

sum(all_genes_JGL[fgl.screen] %in% unique(allTFnames$gene))

screened_gene_data_TAM111.npn<- (sen_JGL_in[[1]])[,fgl.screen]
screened_gene_data_TAM112.npn <- (sen_JGL_in[[2]])[,fgl.screen]
screened_sen_JGL_in <- list(screened_gene_data_TAM111.npn,screened_gene_data_TAM112.npn)

#saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_C.rds")

# saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_C_day14.rds")
# screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_C_day14.rds")
# fgl.results.TAM111up.C.day14 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
#                                     lambda2=0,return.whole.theta = TRUE)
# 
# saveRDS(fgl.results.TAM111up.C.day14, "fgl.results.TAM111up.C.day14.rds")
# print(fgl.results.TAM111up.C.day14)

saveRDS(screened_sen_JGL_in,file="screened_sen_JGL_in_C_day14.rem9.rds")
screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_C_day14.rem9.rds")
fgl.results.TAM111up.C.day14.rem9 <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
                                         lambda2=0,return.whole.theta = TRUE)

saveRDS(fgl.results.TAM111up.C.day14.rem9, "fgl.results.TAM111up.C.day14.rem9.rds")
print(fgl.results.TAM111up.C.day14.rem9)







