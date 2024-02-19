require(Seurat)
library(rlang)
library(future)
library(future.apply)
library(future.callr)
# save(file = "parallel.Seurat.FindMarkers.RData", list = ls())
cp = load("parallel.Seurat.FindMarkers.RData")
cp
pryr::object_size(seurDat[1:6400,], 
                  ident.1 = cells.1,
                  ident.2 = cells.2,
                  min.pct = 0.00000000001,
                  test.use = test,
                  logfc.threshold = 0.000001,
                  slot = Layers(seurDat)[1])

pryr::object_size(seurDat, 
                  ident.1 = cells.1,
                  ident.2 = cells.2,
                  min.pct = 0.00000000001,
                  test.use = test,
                  logfc.threshold = 0.000001,
                  slot = Layers(seurDat)[1])

options("future.globals.maxSize")
options(future.globals.maxSize= 2024^3)
options("future.globals.maxSize")

plan(callr, workers = 12) # works
sEnv = rlang::new_environment(parent = parent.frame())
sEnv$seurDat = seurDat
# seurDat[1:6400,]
ms = multisession( workers = 12, rscript_libs=c("Seurat","future"))
plan(multisession( workers = 8, rscript_libs=c("Seurat","future"))) # Time difference of 4.150328 mins
plan(callr, workers = 4) # Time difference of 4.78646 mins# needs about 4 times the memory of object
plan(callr, workers = 8) # Time difference of 3.717518 mins # needs about 4 times the memory of object
plan(sequential) # Time difference of 15.4785 mins

# full data
plan(callr, workers = 8) # Time difference of 48.82352 mins
plan(callr, workers = 4) # Time difference of 1.090113 hours
plan(sequential) # Time difference of 4.003956 hours
plan(multisession( workers = 8, rscript_libs=c("Seurat","future"))) # they need about 5 GB of RAM each, Time difference of 46.19285 mins
plan(multisession( workers = 4, rscript_libs=c("Seurat","future")))
plan()

{
  start.time = Sys.time()
  markers <- tryCatch.W.E(
    # parallel
    # plan("multisession", workers = 4)
    # plan("multicore", workers = 6)
    # register(MulticoreParam(6))
    
    Seurat::FindMarkers(seurDat, 
                        ident.1 = cells.1,
                        ident.2 = cells.2,
                        min.pct = 0.00000000001,
                        test.use = test,
                        logfc.threshold = 0.000001,
                        slot = Layers(seurDat)[1]
                        # test.use = "wilcox" # p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "bimod"  # p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "roc"    # myAUC   avg_diff power pct.1 pct.2
                        # test.use = "t"      # p_val avg_logFC pct.1 pct.2 p_val_adj
                        # test.use = "negbinom" # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "poisson"  # needs UMI; p_val  avg_logFC pct.1 pct.2    p_val_adj
                        # test.use = "LR"     # p_val  avg_logFC pct.1 pct.2 p_val_adj
                        # test.use = "MAST" # not working: Assay in position 1, with name et is unlogged. Set `check_sanity = FALSE` to override and then proceed with caution.
                        # test.use = "DESeq2" # needs UMI # done separately because the estimating process isn't working with 0s
    )
  )
  Sys.time() - start.time
}
