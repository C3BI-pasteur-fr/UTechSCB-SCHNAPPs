library(monocle3)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(scater)
# scEx object
cp <- load("~/paper1.RData")
scEx1 <- scEx
cp <- load("~/paper2.RData")
common <- intersect(rownames(scEx), rownames(scEx1))
which(colnames(scEx1) %in% colnames(scEx))
colnames(scEx) <- paste0(colnames(scEx), "2")

scEx <- cbind(scEx1[common, ], scEx[common, ])

scEx <- applySimpleR(scEx, HumanPrimaryCellAtlasData(), "humanPrimaryCellAtlas")



colData(scEx)
featureMeta <- rowData(scEx)
cellMeta <- colData(scEx)
cellMeta$cell_names <- cellMeta$barcode
featureMeta$gene_short_name <- featureMeta$symbol
cds <- new_cell_data_set(assays(scEx)[[1]],
  cell_metadata = colData(scEx),
  gene_metadata = featureMeta
)
cds <- preprocess_cds(cds, num_dim = 100)
# batch correction:
cds <- align_cds(cds, alignment_group = "sampleNames")

# dimension reduction
cds <- reduce_dimension(cds, reduction_method = "LSI")

# LSI clustering
cds <- cluster_cells(cds)

# key function
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by = "sampleNames")


cp <- load("~/Downloads/project.2019-10-08.RData")
cds <- new_cell_data_set(assays(scEx)[[1]],
  cell_metadata = colData(scEx),
  gene_metadata = featureMeta
)
for (rdn in reducedDimNames(scEx)) {
  rd <- reducedDim(scEx, type = rdn, withDimnames = TRUE)
  reducedDim(cds, type = rdn) <- rd
}


# reducedDim(cds) <- reducedDim(scEx)
if (!"UMAP" %in% reducedDimNames(cds)) {
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
}

# not working
# monocle3::partitions(cds) <- scEx$dbCluster

cds@clusters[[1]] <- list("UMAP" = list())
cds@clusters[["UMAP"]]$clusters <- scEx$dbCluster
cds@clusters[["UMAP"]]$partitions <- scEx$dbCluster
cds <- monocle3::learn_graph(cds, use_partition = F, verbose = F)

# > names(cds@principal_graph_aux[["UMAP"]])
# [1] "stree"                             "Q"
# [3] "R"                                 "objective_vals"
# [5] "history"                           "pr_graph_cell_proj_closest_vertex"
# [7] "dp_mst"                            "pr_graph_cell_proj_tree"
# [9] "pr_graph_cell_proj_dist"

cds@principal_graph_aux[["UMAP"]]$stree

order_cells(cds1)


names(cds@principal_graph_aux[["UMAP"]])
cds@principal_graph_aux[["UMAP"]]
# clustering
# UMAP clustering
# creates a complex data structure including an igraph and statistics
# cds1 <- cluster_cells(cds, reduction_method = "UMAP", k=20)
# class(cds@clusters[[1]][[1]][[1]]) = igraph
# > names(cds@clusters)
# [1] "UMAP" "PCA"
# > names(cds@clusters[[1]])
# [[1]][[1]] = not needed
# > names(cds@clusters[[1]][[1]])
# [1] "g"          "relations"  "distMatrix"
# "coord"      "edge_links" "optim_res"
# [1] "cluster_result" "partitions"     "clusters"
# partitions: cluster-of-cluster names, factor
# cds@clusters[[1]][[2]] (1-5)  => needed for learn_graph
# partitions: cluster names, factor
# cds@clusters[[1]][[3]] (1-8)

length(cds1@clusters[["UMAP"]]$clusters)


cds1 <- monocle3::learn_graph(cds1, use_partition = F)

cds <- monocle3::learn_graph(cds, use_partition = F)

cds1 <- order_cells(cds1, root_cells = "AGGTCATCATGAGCGA-1")

plot_cells(cds1,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
)
