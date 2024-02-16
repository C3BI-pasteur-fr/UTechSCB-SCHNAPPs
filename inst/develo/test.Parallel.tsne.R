library(SCHNAPPs)
library(Rtsne)
library(SingleCellExperiment)
cp = load(file='~/SCHNAPPsDebug/tsne.RData')
np=15
start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 0

)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 1
)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 2
)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 4
)
end_time <- Sys.time()
end_time - start_time



start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 8
)
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
t=Rtsne::Rtsne(
  pca$x[, 1:np],
  pca = FALSE, dims = gQC_tsneDim,
  perplexity = gQC_tsnePerplexity,
  theta = gQC_tsneTheta,
  check_duplicates = FALSE, num_threads = 16
)
end_time <- Sys.time()
end_time - start_time
