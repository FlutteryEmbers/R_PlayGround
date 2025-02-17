myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

setwd("/scratch/user/r.kapoor/TAM111_112/JGL_scripts/TAM112up/")

#screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_112B.rds")

screened_sen_JGL_in <- readRDS("screened_sen_JGL_in_112B_III.rds")

# Verify
# fgl.screen <- screen.fgl(Y=screened_sen_JGL_in, lambda1=1, lambda2=0, weights = "equal")
# sum(fgl.screen)

fgl.results.TAM111up <- JGL(Y=screened_sen_JGL_in,penalty="fused",lambda1=1,
                            lambda2=0,return.whole.theta = TRUE)

saveRDS(fgl.results.TAM111up, "fgl.results.TAM112up_B_III_2.rds")
print(fgl.results.TAM111up)
