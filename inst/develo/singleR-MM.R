library(SingleR)
hpca.se <- MouseRNAseqData()
library(scater)
load("~/Downloads/puceat-mm-Oct2019-avx.sml.RData")
scEx <- normalize(scEx)
common <- intersect(rownames(scEx), rownames(hpca.se))
hpca.se <- hpca.se[common, ]
scEx <- scEx[common, ]
scEx <- logNormCounts(scEx)

pred.hpca <- SingleR(test = scEx, ref = hpca.se, labels = hpca.se$label.main)
pred.hpca
table(pred.hpca$labels)


plotScoreHeatmap(pred.hpca, show.labels = TRUE)

to.remove <- pruneScores(pred.hpca)
summary(to.remove)

plotScoreDistribution(pred.hpca, show = "delta.med", ncol = 3, show.nmads = 3)

new.pruned <- pred.hpca$labels
new.pruned[pruneScores(pred.hpca, nmads = 5)] <- NA
table(new.pruned, useNA = "always")

library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
