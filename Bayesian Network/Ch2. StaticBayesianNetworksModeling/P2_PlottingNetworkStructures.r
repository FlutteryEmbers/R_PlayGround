# install.packages("BiocManager")
# BiocManager::install("Rgraphviz")
library(bnlearn)
load("Graphs.RData")

hl2 = list(arcs = vstructs(dag2, arcs=TRUE), lwd=4, col="black")
hl3 = list(arcs = vstructs(dag3, arcs=TRUE), lwd=4, col="black")
graphviz.plot(dag2, highlight=hl2, layout="fdp", main="dag2")
graphviz.plot(dag3, highlight=hl3, layout="fdp", main="dag3")
graphviz.plot(cpdag(dag2), highlight=hl2, layout="fdp", main="cpdag(dag2)")
graphviz.plot(cpdag(dag3), highlight=hl3, layout="fdp", main="cpdag(dag3)")

save.image(file = "Graphs.RData")