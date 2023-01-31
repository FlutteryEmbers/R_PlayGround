setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts")

myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

fgl.results.TAM111up_A <- readRDS("fgl.results.TAM111up_A.rds")

degree.TAM111up_A <- net.degree(fgl.results.TAM111up_A$theta)

allTFnames <- read.table("TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allTFnames <- allTFnames[allTFnames$gene %in% colnames(fgl.results.TAM111up_A$theta[[1]]),]

degree.GroupA <- data.frame('TF_name'=allTFnames$gene,'TF_family'=allTFnames$TF_family,'TAM111'=as.numeric(degree.TAM111up_A[[1]][allTFnames$gene]),
           'TAM112'=as.numeric(degree.TAM111up_A[[2]][allTFnames$gene]))

order.degree.GroupA <- order(as.numeric(degree.GroupA[,3]),decreasing = TRUE)

degree.GroupA <- degree.GroupA[order.degree.GroupA,]

degree.GroupA <- cbind(degree.GroupA,'diff'= as.numeric(degree.GroupA[,3]) - as.numeric(degree.GroupA[,4]))
degree.GroupA <- cbind(degree.GroupA,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupA[,3]) - as.numeric(degree.GroupA[,4]))))
View(degree.GroupA)
write.table(file="degree.GroupA.csv",degree.GroupA,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')



fgl.results.TAM111up_B <- readRDS("fgl.results.TAM111up_B.rds")

degree.TAM111up_B <- net.degree(fgl.results.TAM111up_B$theta)

allTFnames <- read.table("TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allTFnames <- allTFnames[allTFnames$gene %in% colnames(fgl.results.TAM111up_B$theta[[1]]),]

degree.GroupB <- data.frame('TF_name'=allTFnames$gene,'TF_family'=allTFnames$TF_family,'TAM111'=as.numeric(degree.TAM111up_B[[1]][allTFnames$gene]),
                            'TAM112'=as.numeric(degree.TAM111up_B[[2]][allTFnames$gene]))

order.degree.GroupB <- order(as.numeric(degree.GroupB[,3]),decreasing = TRUE)

degree.GroupB <- degree.GroupB[order.degree.GroupB,]

degree.GroupB <- cbind(degree.GroupB,'diff'= as.numeric(degree.GroupB[,3]) - as.numeric(degree.GroupB[,4]))
degree.GroupB <- cbind(degree.GroupB,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupB[,3]) - as.numeric(degree.GroupB[,4]))))
View(degree.GroupB)
write.table(file="degree.GroupB.csv",degree.GroupB,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

#Group C

fgl.results.TAM111up_C <- readRDS("fgl.results.TAM111up_C.rds")

degree.TAM111up_C <- net.degree(fgl.results.TAM111up_C$theta)

allTFnames <- read.table("TFs_v1.1.csv",sep=',',header=TRUE, stringsAsFactors = FALSE)
allTFnames <- allTFnames[allTFnames$gene %in% colnames(fgl.results.TAM111up_C$theta[[1]]),]

degree.GroupC <- data.frame('TF_name'=allTFnames$gene,'TF_family'=allTFnames$TF_family,'TAM111'=as.numeric(degree.TAM111up_C[[1]][allTFnames$gene]),
                            'TAM112'=as.numeric(degree.TAM111up_C[[2]][allTFnames$gene]))

order.degree.GroupC <- order(as.numeric(degree.GroupC[,3]),decreasing = TRUE)

degree.GroupC <- degree.GroupC[order.degree.GroupC,]

degree.GroupC <- cbind(degree.GroupC,'diff'= as.numeric(degree.GroupC[,3]) - as.numeric(degree.GroupC[,4]))
degree.GroupC <- cbind(degree.GroupC,'abs_diff'= as.numeric(abs(as.numeric(degree.GroupC[,3]) - as.numeric(degree.GroupC[,4]))))
View(degree.GroupC)
write.table(file="degree.GroupC.csv",degree.GroupC,row.names = FALSE,col.names =TRUE,quote=FALSE,sep=',')

#######
#day14 only and day5 only

###day14 only: Group A
fgl.results.TAM111up.A.day14 <- readRDS("fgl.results.TAM111up.A.day14.rds")
degrees.A.day14 <- net.degree(fgl.results.TAM111up.A.day14$theta)
write.table(file="TAM111hubs.day14.GroupA.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[1]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[1]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="TAM112hubs.day14.GroupA.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[2]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[2]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')

###day14 only: Group B
fgl.results.TAM111up.A.day14 <- readRDS("fgl.results.TAM111up.B.day14.rds")
degrees.A.day14 <- net.degree(fgl.results.TAM111up.A.day14$theta)
write.table(file="TAM111hubs.day14.GroupB.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[1]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[1]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="TAM112hubs.day14.GroupB.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[2]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[2]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')

###day14 only: Group C
fgl.results.TAM111up.A.day14 <- readRDS("fgl.results.TAM111up.C.day14.rds")
degrees.A.day14 <- net.degree(fgl.results.TAM111up.A.day14$theta)
write.table(file="TAM111hubs.day14.GroupC.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[1]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[1]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="TAM112hubs.day14.GroupC.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[2]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[2]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')


###day14 only: Group B rem 9
fgl.results.TAM111up.A.day14 <- readRDS("fgl.results.TAM111up.B.day14.rem9.rds")
degrees.A.day14 <- net.degree(fgl.results.TAM111up.A.day14$theta)
write.table(file="TAM111hubs.day14.GroupB.rem9.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[1]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[1]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')
write.table(file="TAM112hubs.day14.GroupB.rem9.csv",cbind('TAM111'=degrees.A.day14[[1]][degrees.A.day14[[2]]>=10],
                                                     'TAM112'=degrees.A.day14[[2]][degrees.A.day14[[2]]>=10]),
            row.names = TRUE,col.names =TRUE,quote=FALSE,sep=',')



