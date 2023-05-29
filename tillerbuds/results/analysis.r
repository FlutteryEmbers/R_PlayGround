library(JGL)
dat <- readRDS("fgl.results.tillerBud.rds")
net.neighbors(dat$theta[2], 'SORBI_3007G165701')
edges <- net.edges(dat$theta[2])

net.hubs(dat$theta[2], nhub = 10)
net.hubs(dat$theta[1], nhub = 10)

LeafIntact_network <- dat$theta[1]
LeafRemoved_network <- dat$theta[2]
LeafRemoved_network[[1]]['SORBI_3007G165701', 'SORBI_3001G293800']



neighbor1= c()
neighbor2 = c()
neighbor_diff = c()
abs_neighbor_diff = c()

gene_names <- names(dat2[[1]][,1])
for (gene_name in gene_names){
  neighbor1 = c(neighbor1, lengths(net.neighbors(dat$theta[1], gene_name)))
  neighbor2 = c(neighbor2, lengths(net.neighbors(dat$theta[2], gene_name)))
  difference = lengths(net.neighbors(dat$theta[1], gene_name)) - lengths(net.neighbors(dat$theta[2], gene_name))
  neighbor_diff = c(neighbor_diff, difference)
  abs_neighbor_diff = c(abs_neighbor_diff, abs(difference))
}

table <- data.frame(
  gene_name = gene_names,
  neighbor1 = neighbor1,
  neighbor2 = neighbor2,
  neighbor_diff = neighbor_diff,
  abs_neighbor_diff = abs_neighbor_diff
)
write.csv(table, "./results.csv", row.names=TRUE)

table[order(neighbor1, decreasing=TRUE), ][1:10, ]

gene_interest = 'SORBI_3009G247900'
neighbor_strength_1 = c()
neighbor_strength_2 = c()
strength_diffs = c()
abs_strenghth_diffs = c()
for (gene_name in gene_names){
  neighbor_strength_1 = c(neighbor_strength_1, LeafIntact_network[[1]][gene_interest, gene_name])
  neighbor_strength_2 = c(neighbor_strength_2, LeafRemoved_network[[1]][gene_interest, gene_name])
  strength_diff = LeafIntact_network[[1]][gene_interest, gene_name] - LeafRemoved_network[[1]][gene_interest, gene_name]
  strength_diffs = c(strength_diffs, strength_diff)
  abs_strenghth_diffs = c(abs_strenghth_diffs, abs(strength_diff))
}

stress_table <- data.frame(
  gene_name = gene_names,
  neighbor_strength_1 = neighbor_strength_1,
  neighbor_strength_2 = neighbor_strength_2,
  strength_diffs = strength_diffs,
  abs_strenghth_diffs = abs_strenghth_diffs
)
write.csv(table, "./stress_results.csv", row.names=TRUE)

   
