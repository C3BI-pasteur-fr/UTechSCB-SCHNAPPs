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
Sys.setenv(LD_LIBRARY_PATH=paste("/usr/lib/", "/opt/local/lib/", "/opt/local/", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
Sys.setenv(HDF5_CFLAGS="-I/opt/local/include/")
Sys.setenv(HDF5_LIBS="-L/opt/local/lib -lhdf5 -Wl,-rpath,/usr/lib/")
Sys.setenv(HDF5_LIBS="-Wl,-rpath,/usr/lib/")
#Sys.setenv(PATH = paste("/opt/homebrew/bin/", Sys.getenv("PATH"), sep = .Platform$path.sep))
# these might be set in sources if compilation fails
HDF5_CFLAGS="-I/opt/local/include/"
HDF5_LIBS="-L/opt/local/lib -lhdf5"

#Warning: unknown option ‘--configure-vars=PKG_LIBS=-Wl,-rpath,/usr/lib/’

remotes::install_github("bnprks/BPCells/r", 
                        configure.vars=c("PKG_CXXFLAGS='-Wl,-rpath,/usr/lib/'", CC="clang -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5 -Wl,-rpath,/usr/lib/"),
                        build_opts = c( "--no-resave-data", "--no-manual", "--no-build-vignettes"), 
                        verbose=T,dependencies = T, INSTALL_opts = c("--with-keep.source", "--install-tests"))


setRepositories(ind=1:6)
options(repos="http://cran.rstudio.com/")
if(!require(devtools)) { install.packages("devtools") }
library(devtools) 

install_github("cBioPortal/cgdsr")
devtools::install_github("RausellLab/CelliD", ref = "legacy")
devtools::install_github("rcannood/SCORPIUS")

devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE)




