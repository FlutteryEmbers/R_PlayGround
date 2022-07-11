getwd();
workingDir = ".";
setwd(workingDir);
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads();
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6;
adjacency = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1 - TOM;

geneTree = hclust(as.dist(dissTOM), method = 'average');
sizeGrWindow(12, 9)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilartiy", labels = FALSE, hang = 0.04);

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, disM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = False, minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamics"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

ModuleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(ModuleColors, colorOrder) - 1;
MEs = mergedMEs;
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")



