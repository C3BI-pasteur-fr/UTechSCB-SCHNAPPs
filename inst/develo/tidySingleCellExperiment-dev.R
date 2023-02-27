library(tidySingleCellExperiment)
cp = load(file = "scEx1.RData")
scEx
colData(scEx)
tsce = scEx %>% group_by(dbCluster) %>% slice_sample(n=50) %>% ungroup()

class(tsce)
scEx[,tsce$.cell] %>% class()
subSce = scEx[,tsce$.cell] 
colData(subSce)

all(colnames(subSce) == tsce$.cell)
all(rownames(colData(subSce))== tsce$.cell)

