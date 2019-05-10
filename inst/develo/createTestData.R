load("data/scEx.RData")
dat = rowData(scEx)

rownames(dat) = dat$symbol
colnames (dat) = c("Descr1","ChrName1","Gene.St1","Gene.E1","Gene.Name1","biotype1","symbol","id" )
head(dat)
write.file.csv(dat, row.names=TRUE, file="data/scExGenesRowNames.csv" )


load("data/scEx.RData")
dat = rowData(scEx)
