myPaths <- .libPaths()
myPaths<-c("/scratch/user/r.kapoor/Sorghum/DESeq_analysis/featureCounts",myPaths)
.libPaths(myPaths)
library(JGL)

setwd("/scratch/user/r.kapoor/Sorghum/tiller_bud")


# Verify
screened_tillerBud_JGL_in <- readRDS("screened_tillerBud_JGL_in.rds")
head(rownames(t(screened_tillerBud_JGL_in[[1]])))
#fgl.results.tillerBud <- JGL(Y=screened_tillerBud_JGL_in,penalty="fused",lambda1=1,lambda2=0,return.whole.theta = TRUE)
#saveRDS(fgl.results.tillerBud, "fgl.results.tillerBud.rds")
#print(fgl.results.tillerBud)
