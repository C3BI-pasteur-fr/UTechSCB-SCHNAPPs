update.packages()

if (!require("devtools"))
  install.packages("devtools")
# devtools::install_github("mul118/shinyMCE")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# update bioconductor packages if required
BiocManager::install()
# BiocManager::install("BiocSingular")
BiocManager::install("InteractiveComplexHeatmap")
BiocManager::install("SingleR") # needed to prepare data.
devtools::install_github("haowulab/Wind")
devtools::install_github("nghiavtr/BPSC")
BiocManager::install("DEsingle")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("MAST")
BiocManager::install("monocle")
BiocManager::install("monocle3")
BiocManager::install("scDD")
BiocManager::install("limma")
BiocManager::install("Seurat")
devtools::install_github("statOmics/zingeR")
BiocManager::install("scater")
install.packages("aggregation")
devtools::install_github("cysouw/qlcMatrix")

devtools::install_github("Zhangxf-ccnu/scDEA")

BiocManager::install("InteractiveComplexHeatmap")
install.packages("bookdown")
BiocManager::install ("GSVA")
BiocManager::install("GSEABase")

devtools::install_github("BaderLab/Tempora")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github('saezlab/liana')
Sys.setenv(LD_LIBRARY_PATH=paste("/opt/local/lib/", "/opt/local/", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
Sys.setenv(HDF5_CFLAGS="-I/opt/local/include/")

# these might be set in sources if compilation fails
HDF5_CFLAGS="-I/opt/local/include/"
HDF5_LIBS="-L/opt/local/lib -lhdf5"

remotes::install_github("bnprks/BPCells", 
                        configure.vars=c(CC="clang -arch x86_64 -I /opt/local/include/"),
                        build_opts = c( "--no-resave-data", "--no-manual", "--no-build-vignettes"))


setRepositories(ind=1:6)
options(repos="http://cran.rstudio.com/")
if(!require(devtools)) { install.packages("devtools") }
library(devtools) 

install_github("cBioPortal/cgdsr")
devtools::install_github("RausellLab/CelliD", ref = "legacy")
devtools::install_github("rcannood/SCORPIUS")

devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE)
