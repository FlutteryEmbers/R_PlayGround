library(bnlearn)
data(marks)
str(marks)

# Assign Edges
ug = empty.graph(names(marks))
arcs(ug, check.cycles = TRUE) = matrix(
  c("MECH", "VECT", "MECH", "ALG",  "VECT", "MECH",
    "VECT", "ALG",  "ALG",  "MECH", "ALG",  "VECT",
    "ALG",  "ANL", "ALG",   "STAT", "ANL",  "ALG",
    "ANL",   "STAT", "STAT", "ALG", "STAT", "ANL"), 
  ncol = 2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))
ug

dag = empty.graph(names(marks))
arcs(dag) = matrix(
  c("VECT", "MECH", "ALG", "MECH", "ALG", "VECT",
    "ANL", "ALG", "STAT", "ALG", "STAT", "ANL"),
  ncol = 2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))
dag

# Create Graph By Adjacency Matrix
mat = matrix(
  c(0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 
  nrow=5, 
  dimnames = list(nodes(dag), nodes(dag)))
mat
dag2 = empty.graph(nodes(dag))
amat(dag2) = mat
all.equal(dag, dag2)

dag3 = empty.graph(nodes(dag))
dag3 = set.arc(dag3, "VECT", "MECH")
dag3 = set.arc(dag3, "ALG", "MECH")
dag3 = set.arc(dag3, "ALG", "VECT")
dag3 = set.arc(dag3, "ANL", "ALG")
dag3 = set.arc(dag3, "STAT", "ALG")
dag3 = set.arc(dag3, "STAT", "ANL")
all.equal(dag, dag3)

# moral of a graph
all.equal(ug, moral(dag))

node.ordering(dag)

# provide a synthetic description of the local dependence structure around that node.
nbr(dag, "ANL")
mb(dag, "ANL")

# show that both sets describe symmetric relationships
"ANL" %in% mb(dag, "ALG")
"ALG" %in% mb(dag, "ANL")

chld = children(dag, "VECT")
par = parents(dag, "VECT")
o.par = sapply(chld, parents, x = dag) # Find parents of a given node
unique(c(chld, par, o.par[o.par != "VECT"]))
mb(dag, "VECT")

# log-likelihood score of a Bayesian network
score(dag, data = marks, type = "loglik-g")
dag.eq = reverse.arc(dag, "STAT", "ANL")
score(dag.eq, data = marks, type = "loglik-g")

# V-structure
vstructs(dag)
vstructs(dag.eq)

dag2 = drop.arc(dag, from="STAT", to="ANL")
dag3 = drop.arc(dag, from="ALG", to="VECT")
vstructs(dag2)
vstructs(dag3)
all.equal(cpdag(dag2), cpdag(dag3))

all.equal(moral(dag2), moral(dag3))
all.equal(moral(dag2), moral(dag))
all.equal(moral(dag3), moral(dag))

save.image(file = "Graphs.RData")
