library(bnlearn)
library(pcalg)

load("Graphs.RData")
bn.gs = gs(marks)
bn.gs

all.equal(bn.gs, iamb(marks))
all.equal(bn.gs, inter.iamb(marks))
all.euqal(bn.gs, iamb(marks, text = "mc-cor"))

# suff.stat = list(C=cor(marks), n = nrow(marks))
# pc.fit = pc(suff.stat, indepTest = gaussCItest, p = ncol(marks), alpha = 0.05)
# pc.fit

gs.graph = as.graphAM(bn.gs)
# compareGraphs(pc.fit@graph, gs.graph)


