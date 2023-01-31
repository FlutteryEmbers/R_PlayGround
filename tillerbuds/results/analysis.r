library(JGL)
dat <- readRDS("fgl.results.tillerBud.rds")
net.neighbors(dat$theta[2], 'SORBI_3007G165701')
edges <- net.edges(dat$theta[2])

net.hubs(dat$theta[2], nhub = 10)
net.hubs(dat$theta[1], nhub = 10)

dat2 = dat$theta[1]


neighbor1= c()
neighbor2 = c()
deg_diff = c()
abs_deg_diff = c()

gene_names <- names(dat2[[1]][,1])
for (gene_name in gene_names){
  neighbor1 = c(neighbor1, lengths(net.neighbors(dat$theta[1], gene_name)))
  neighbor2 = c(neighbor2, lengths(net.neighbors(dat$theta[2], gene_name)))
  difference = lengths(net.neighbors(dat$theta[1], gene_name)) - lengths(net.neighbors(dat$theta[2], gene_name))
  deg_diff = c(deg_diff, difference)
  abs_deg_diff = c(abs_deg_diff, abs(difference))
}

table <- data.frame(
  gene_name = gene_names,
  neighbor1 = neighbor1,
  neighbor2 = neighbor2,
  deg_diff = deg_diff,
  abs_deg_diff = abs_deg_diff
)

table[order(neighbor1, decreasing=TRUE), ][1:10, ]

   
