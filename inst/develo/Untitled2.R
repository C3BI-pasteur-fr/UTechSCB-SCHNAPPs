# scatr qc

cellStats <- perCellQCMetrics(scEx)

colnames(colData(scEx))

colnames(rowData(scEx))


plotRowData(scEx)

plotColData(scEx,
  x = "total_features_by_counts",
  y = "pct_counts_feature_control", colour = "Mutation_Status"
) +
  theme(legend.position = "top") +
  stat_smooth(method = "lm", se = FALSE, size = 2, fullrange = TRUE)

plotScater(scEx,
  block1 = "total_features_by_counts", block2 = "total_counts",
  colour_by = "pct_counts_in_top_50_features", nfeatures = 300, exprs_values = "counts"
)
