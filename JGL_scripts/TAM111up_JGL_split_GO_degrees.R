setwd("/scratch/user/r.kapoor/wheat/wgcna_grain/grain_RNAseq")

OntologiesForGenes <- readRDS("/scratch/user/r.kapoor/wheat/wgcna_grain/grain_RNAseq/OntologiesForGenes.rds")

OntologiesForGenes <- OntologiesForGenes[OntologiesForGenes$ontology == "GO",]

#amyloplast_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0009501",1], as.character))))
#amyloplast_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:2001070",1], as.character))))
amyloplast_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0019252",1], as.character))))


##Group A
setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts")

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

fgl.results.TAM111up_A <- readRDS("fgl.results.TAM111up_A.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- amyloplast_Genes
allTFnames <- allTFnames[allTFnames %in% colnames(fgl.results.TAM111up_A$theta[[1]])]

allTFnames <- names(degree.TAM111up_A[[1]])

degree.GroupA <- data.frame('Gene_name'=allTFnames,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames]),
                            'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames]))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,2]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

View(degree.GroupA)
#write.table(file="GO_0009501.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
#write.table(file="GO_2001070.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="GO_0019252.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

write.table(file="allgenes.degree.TAM111GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

##Grp B
setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts")
fgl.results.TAM111up_A <- readRDS("fgl.results.TAM111up_B.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- amyloplast_Genes
allTFnames <- allTFnames[allTFnames %in% colnames(fgl.results.TAM111up_A$theta[[1]])]

allTFnames <- names(degree.TAM111up_A[[1]])

degree.GroupA <- data.frame('Gene_name'=allTFnames,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames]),
                            'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames]))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,2]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

View(degree.GroupA)
#write.table(file="GO_0009501.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
#write.table(file="GO_2001070.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="GO_0019252.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

write.table(file="allgenes.degree.TAM111GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

##Grp C
setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts")
fgl.results.TAM111up_A <- readRDS("fgl.results.TAM111up_C.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- amyloplast_Genes
allTFnames <- allTFnames[allTFnames %in% colnames(fgl.results.TAM111up_A$theta[[1]])]

allTFnames <- names(degree.TAM111up_A[[1]])

degree.GroupA <- data.frame('Gene_name'=allTFnames,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames]),
                            'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames]))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,2]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

View(degree.GroupA)
#write.table(file="GO_0009501.degree.GroupC.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
#write.table(file="GO_2001070.degree.GroupC.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="GO_0019252.degree.GroupC.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')


write.table(file="allgenes.degree.TAM111GroupC.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')


