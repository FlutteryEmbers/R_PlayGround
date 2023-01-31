
.libPaths(c(.libPaths(),"/scratch/user/r.kapoor/my_R_libs/"))
setwd("/scratch/user/r.kapoor/tiller_bud")

library(DESeq2)
rlog_data <- readRDS("/scratch/user/r.kapoor/tiller_bud/rlog_data.rds")


coldata_Intact <- rlog_data$treatment=="Intact"
coldata_LeafRemoved <- rlog_data$treatment=="LeafRemoved"

gene_data_Intact.npn<-t(as.matrix(assay(rlog_data)[,coldata_Intact]))
gene_data_LeafRemoved.npn<-t(as.matrix(assay(rlog_data)[,coldata_LeafRemoved]))

#library(JGL,lib.loc = "/scratch/user/r.kapoor/my_R_libs/")
library(JGL)
tillerBud_JGL_in <- list(gene_data_Intact.npn,gene_data_LeafRemoved.npn)
fgl.screen <- screen.fgl(Y=tillerBud_JGL_in, lambda1=1, lambda2=0, weights = "equal")
sum(fgl.screen)

screened_gene_data_Intact.npn<- (tillerBud_JGL_in[[1]])[,fgl.screen]
screened_gene_data_LeafRemoved.npn <- (tillerBud_JGL_in[[2]])[,fgl.screen]
screened_tillerBud_JGL_in <- list(screened_gene_data_Intact.npn,screened_gene_data_LeafRemoved.npn)

saveRDS(screened_tillerBud_JGL_in,file="screened_tillerBud_JGL_in.rds")

fgl.results.tillerBud <- JGL(Y=screened_tillerBud_JGL_in,penalty="fused",lambda1=1,
                            lambda2=0,return.whole.theta = TRUE)

saveRDS(fgl.results.tillerBud, "fgl.results.tillerBud.rds")
print(fgl.results.tillerBud)
