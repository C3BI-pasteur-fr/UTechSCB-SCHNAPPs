# pairwise correlation
# take way too long!
corrP = correlatePairs(scEx_log)

# bootstrap cluster
library(pheatmap)
library(bluster)
ass.prob <- bootstrapStability(scEx_log, FUN=function(x) {
  g <- buildSNNGraph(x)
  igraph::cluster_walktrap(g)$membership
}, clusters=scEx_log@colData$dbCluster)

pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=colorRampPalette(c("white", "blue"))(100))

InteractiveHeatmap

fps = dir(path = "/Volumes/LaCie2022/RStudio_history/katja/hist_2022-Jul-10.13.04/", full.names = T, pattern = "*.RData")
for (fp in fps){
  cp = load(fp)
  if("inp" %in% cp){
    if(is.null(inp)){
      inp = inputS
      cat(file = stderr(), fp)
      save(file = fp, list = cp)
    }
  }
}
cp = load("/Volumes/LaCie2022/RStudio_history/katja/hist_2022-Jul-10.13.04/scEx_log.c91b68a4e307.RData")