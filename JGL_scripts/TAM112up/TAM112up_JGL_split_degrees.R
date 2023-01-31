setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts/TAM112up/")

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

fgl.results.TAM111up_A <- readRDS("fgl.results.TAM112up_A.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- read.table("../TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allTFnames <- allTFnames[allTFnames$gene %in% colnames(fgl.results.TAM111up_A$theta[[1]]),]

degree.GroupA <- data.frame('TF_name'=allTFnames$gene,'TF_family'=allTFnames$TF_family,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames$gene]),
           'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames$gene]))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,4]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,3]) - as.numeric(degree.GroupA[,4]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,3]) - as.numeric(degree.GroupA[,4]))))
View(degree.GroupA)
write.table(file="degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')



fgl.results.TAM111up_B <- readRDS("fgl.results.TAM112up_B_III.rds")

degree.TAM111up_B <- net.degree(fgl.results.TAM111up_B$theta)

allTFnames <- read.table("../TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allTFnames <- allTFnames[allTFnames$gene %in% colnames(fgl.results.TAM111up_B$theta[[1]]),]

degree.GroupB <- data.frame('TF_name'=allTFnames$gene,'TF_family'=allTFnames$TF_family,'TAM111'=as.numeric(degree.TAM111up_B[[1]][allTFnames$gene]),
                            'TAM112'=as.numeric(degree.TAM111up_B[[2]][allTFnames$gene]))

order.degree.GroupB <- order(as.numeric(degree.GroupB[,4]),decreasing = TRUE)

degree.GroupB <- degree.GroupB[order.degree.GroupB,]

degree.GroupB <- cbind(degree.GroupB,'diff'= as.numeric(degree.GroupB[,3]) - as.numeric(degree.GroupB[,4]))
degree.GroupB <- cbind(degree.GroupB,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupB[,3]) - as.numeric(degree.GroupB[,4]))))

order.degree.GroupB <- order(as.numeric(degree.GroupB[,5]))

degree.GroupB <- degree.GroupB[order.degree.GroupB,]

View(degree.GroupB)
write.table(file="degree.GroupB.csv",degree.GroupB,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

