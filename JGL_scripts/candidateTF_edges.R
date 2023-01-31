
myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

fgl.results.TAM111up_A <- readRDS("fgl.results.TAM111up_A.rds")

options(stringsAsFactors = FALSE)
candidate_TF <- read.table("../ARACNe/candidate_TFs.txt")

for (TFid in 1:length(candidate_TF$V1)){
write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_A$theta)[[classid]][candidate_TF$V1[TFid]]),
net.neighbors(fgl.results.TAM111up_A$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupA.TAM111edges.csv",
quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}

classid <- 1;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_A$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_A$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupA.TAM111edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}
classid <- 2;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_A$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_A$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupA.TAM112edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}

fgl.results.TAM111up_B <- readRDS("fgl.results.TAM111up_B.rds")
classid <- 1;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_B$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_B$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupB.TAM111edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}
classid <- 2;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_B$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_B$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupB.TAM112edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}

fgl.results.TAM111up_C <- readRDS("fgl.results.TAM111up_C.rds")
classid <- 1;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_C$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_C$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupC.TAM111edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}
classid <- 2;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM111up_C$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM111up_C$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM111up.GroupC.TAM112edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}

fgl.results.TAM112up_A <- readRDS("TAM112up/fgl.results.TAM112up_A.rds")
classid <- 1;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM112up_A$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM112up_A$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM112up.GroupA.TAM111edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}
classid <- 2;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM112up_A$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM112up_A$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM112up.GroupA.TAM112edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}

fgl.results.TAM112up_B <- readRDS("TAM112up/fgl.results.TAM112up_B_III.rds")
classid <- 1;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM112up_B$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM112up_B$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM112up.GroupB.TAM111edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}
classid <- 2;
for (TFid in 1:length(candidate_TF$V1)){write.table(print(cbind(rep(candidate_TF$V1[TFid],net.degree(fgl.results.TAM112up_B$theta)[[classid]][candidate_TF$V1[TFid]]),net.neighbors(fgl.results.TAM112up_B$theta,candidate_TF$V1[TFid])[[classid]])), file="candidateTF.TAM112up.GroupB.TAM112edges.csv",quote=FALSE,col.names=FALSE,row.names=FALSE,sep=',',append=TRUE)}