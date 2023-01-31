setwd("/scratch/user/r.kapoor/wheat/wgcna_grain/grain_RNAseq")

OntologiesForGenes <- readRDS("/scratch/user/r.kapoor/wheat/wgcna_grain/grain_RNAseq/OntologiesForGenes.rds")

OntologiesForGenes <- OntologiesForGenes[OntologiesForGenes$ontology == "GO",]

#amyloplast_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0009501",1], as.character))))
#amyloplast_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:2001070",1], as.character))))
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0071705",1], as.character))))
#protein transport
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0015031",1], as.character))))
#protein localization: GO:0008104
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0008104",1], as.character))))
#introcellular protein transport: GO:0006886
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0006886",1], as.character))))
####GO:0046488 
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0046488",1], as.character))))
#phospholipid biosynthetic process - GO:0008654
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0008654",1], as.character))))
# Glycerophospholipid metabolic process: GO:0006650
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[OntologiesForGenes$ID=="GO:0006650",1], as.character))))
# GO:0006644 phospholipid metabolic
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[
  OntologiesForGenes$ID=="GO:0006644",1], as.character))))
# GO:0034613: cellular protein loca
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[
  OntologiesForGenes$ID=="GO:0034613",1], as.character))))
# GO:0051603
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[
  OntologiesForGenes$ID=="GO:0051603",1], as.character))))
#GO:0033365
nitrogenCompTransport_Genes <- unique(gsub('01G','02G',unlist(lapply(OntologiesForGenes[
  OntologiesForGenes$ID=="GO:0033365",1], as.character))))


###
nitrogenCompTransport_Genes <- read.table("GO_terms/GO-0046488.txt",header=TRUE,colClasses=c("character"))
nitrogenCompTransport_Genes <- nitrogenCompTransport_Genes$Gene_stable_ID


###
nitrogenCompTransport_Genes <- read.table("GO_terms/GO-0034613.txt",header=TRUE,colClasses=c("character"))
nitrogenCompTransport_Genes <- nitrogenCompTransport_Genes$Gene_stable_ID

###
###
nitrogenCompTransport_Genes <- read.table("GO_terms/GO-0015031.txt",header=TRUE,colClasses=c("character"))
nitrogenCompTransport_Genes <- nitrogenCompTransport_Genes$Gene_stable_ID
###
###
nitrogenCompTransport_Genes <- read.table("GO_terms/GO-0006644.txt",header=TRUE,colClasses=c("character"))
nitrogenCompTransport_Genes <- nitrogenCompTransport_Genes$Gene_stable_ID
###
nitrogenCompTransport_Genes <- read.table("GO_terms/GO-0051603.txt",header=TRUE,colClasses=c("character"))
nitrogenCompTransport_Genes <- nitrogenCompTransport_Genes$Gene_stable_ID

##Group A
setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts/TAM112up/")

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

fgl.results.TAM111up_A <- readRDS("fgl.results.TAM112up_A.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- nitrogenCompTransport_Genes
allTFnames <- allTFnames[allTFnames %in% colnames(fgl.results.TAM111up_A$theta[[1]])]

allTFnames <- names(degree.TAM111up_A[[1]])


degree.GroupA <- data.frame('Gene_name'=allTFnames,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames]),
                            'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames]))



degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]))

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

View(degree.GroupA)

write.table(file="allgenes.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')


#write.table(file="GO_0009501.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
#write.table(file="GO_2001070.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="GO_terms/GO-0046488.degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

##Grp B
setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts/TAM112up/")
fgl.results.TAM111up_A <- readRDS("fgl.results.TAM112up_B_III.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- nitrogenCompTransport_Genes
allTFnames <- allTFnames[allTFnames %in% colnames(fgl.results.TAM111up_A$theta[[1]])]

allTFnames <- names(degree.TAM111up_A[[1]])

degree.GroupA <- data.frame('Gene_name'=allTFnames,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames]),
                            'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames]))



degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,2]) - as.numeric(degree.GroupA[,3]))))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]))

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

View(degree.GroupA)
#write.table(file="GO_0009501.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
#write.table(file="GO_2001070.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
##write.table(file="GO-0051603.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="GO_terms/GO-0046488.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

write.table(file="allgenes.degree.GroupB.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

