
load("/Volumes/@Single_cell/SingleCellCourse_2021/schnapps/Tcells_d2_10X.schnapps.RData")
library(phateR)

plot(prcomp(t(assays(scEx)[[1]])$x))

system.time(
tree.phate <- phate(t(assays(scEx)[[1]])$x))
)
summary(tree.phate)


palette(rainbow(10))
plot(tree.phate)

tree.phate <- phate(tree.data$data, gamma=0, t=120, init=tree.phate)
# plot embedding
palette(rainbow(10))
plot(tree.phate)


data(tree.data)
tree.data$branches

