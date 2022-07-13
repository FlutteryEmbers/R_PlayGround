# install.packages("bnlearn")
library(bnlearn)
data("lizards")
# lizards = read.table("lizards.txt", header = TRUE)
str(lizards)
summary(lizards)

levels(lizards[, "Species"])
levels(lizards$Species)

levels(lizards[, "Height"])
levels(lizards$Height)

levels(lizards[, "Diameter"])
levels(lizards$Diameter)

table(lizards[, c(3, 2, 1)])

Sagrei.lizards = lizards[lizards$Species == 'Sagrei', ]
Distichus.lizards = lizards[lizards$Species == 'Distichus', ]
par(mfrow = c(2, 2))
plot(Sagrei.lizards$Height)
plot(Distichus.lizards$Height)
plot(Sagrei.lizards$Diameter)
plot(Distichus.lizards$Diameter)

diam = numeric(length = nrow(lizards))
narrow = (lizards$Diameter == "narrow")
wide = (lizards$Diameter == "wide")
diam[narrow] = runif(n = 252, min = 2, max = 4)
diam[wide] = runif(n = 157, min = 4, max = 6)
new.data = data.frame(Species = lizards$Species, Sim.Diameter = diam)

summary(new.data$Sim.Diameter)

is.sagrei = (new.data$Species == "Sagrei")
summary(new.data[is.sagrei, "Sim.Diameter"])
summary(new.data[!is.sagrei, "Sim.Diameter"])
var(new.data[is.sagrei, "Sim.Diameter"])
var(new.data[!is.sagrei, "Sim.Diameter"])

boxplot(Sim.Diameter ~ Species, data = new.data, ylab = "Diameter(inches)")
abline(h  = 4, lty = "dashed")


