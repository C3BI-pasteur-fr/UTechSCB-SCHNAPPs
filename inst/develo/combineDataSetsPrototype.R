cp = load("/Volumes/LaCie2022/scData/scEx.51e3d27df0e E13,5.RData")
cp
scEx2 = scEx$scEx
cp = load("/Volumes/LaCie2022/scData/scEx.51e7382b945 E16,5.RData")
cp
scEx= scEx$scEx
allCols = union(colnames(colData(scEx2)), colnames(colData(scEx)))
colData(scEx2)[,setdiff(allCols,colnames(colData(scEx2)) )] = NA
colData(scEx)[,setdiff(allCols,colnames(colData(scEx)) )] = NA
newcDat = DataFrame(rbind(colData(scEx2),colData(scEx)))

commGenes = intersect(rownames(scEx), rownames(scEx2))
mat =  cbind(as.data.frame(assays(scEx)[["counts"]])[commGenes,], as.data.frame(assays(scEx2)[["counts"]])[commGenes,]) %>% as.matrix %>% as("CsparseMatrix")
mat

allCols = union(colnames(rowData(scEx2)), colnames(rowData(scEx)))
rowData(scEx2)[,setdiff(allCols,colnames(rowData(scEx2)) )] = NA
rowData(scEx)[,setdiff(allCols,colnames(rowData(scEx)) )] = NA
newrDat = DataFrame(rbind(rowData(scEx2),rowData(scEx)))

newrDat
combscEx <- SingleCellExperiment(
  assay = list(counts = mat),
  colData = newcDat[colnames(mat),],
  rowData = newrDat[commGenes,]
)

save(file = "/Volumes/LaCie2022/scData/smd.e16.RData", list = c("combscEx"))
