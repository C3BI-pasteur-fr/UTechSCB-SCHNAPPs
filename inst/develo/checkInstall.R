# not sure this is needed
# BiocManager::install("diffcyt")
# 
# BiocManager::install("HDCytoData")
# BiocManager::install("CATALYST")
# 
# devtools::install_github("hrbrmstr/dtupdate")
# library(dtupdate)
# dtupdate::github_update()

# brew install poppler librsvg rust leptonica



install.packages("devtools", build = T,type = "source", upgrade = "always")
install.packages("BiocManager")

# this needs FC defined in Makevars (nothing else works)
install.packages("mnormt", 
                 dependencies  = TRUE, 
                 build = T,
                 type = "source", 
                 upgrade = "always",
                 configure.args  = "FC=/usr/local/bin/gfortran",
                 verbose = T
)

# FC=/usr/local/bin/gfortran-13 
# F77=/usr/local/bin/gfortran-13 

install.packages("psychTools", 
                 dependencies  = TRUE, 
                 build = T,
                 type = "source", 
                 upgrade = "always",
                 configure.args  = "FC=/usr/local/bin/gfortran",
                 verbose = T
)

remotes::install_github("satijalab/seurat", dependencies  = TRUE, build = T,type = "source", upgrade = "always")
remotes::install_github("satijalab/seurat-data", quiet = TRUE, dependencies  = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install("scran", dependencies  = TRUE, build = T, update = TRUE, ask = F, type = "source", upgrade = "always")
devtools::install_github("kassambara/ggpubr", dependencies  = TRUE, build = T,type = "source", upgrade = "always")

# ‘psychTools’, ‘Seurat’, ‘scran’, ‘ggpubr’ are not available for package ‘SCHNAPPs’
# Warning messages:
#   1: packages ‘multtest’, ‘limma’, ‘monocle’, ‘rtracklayer’, ‘MAST’ are not available for this version of R


devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE, build = T,type = "source", upgrade = "always")


# 
# 
# PATH=Sys.getenv("PATH")
# PATH = paste0("/usr/local/bin:",PATH)
# Sys.setenv(PATH = PATH)
# PKG_CONFIG_PATH = "/usr/local/lib/pkgconfig/tesseract.pc"
# Sys.setenv(PKG_CONFIG_PATH = PKG_CONFIG_PATH)
# Sys.getenv("PKG_CONFIG_PATH")

# R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'
# sudo port install leptonica
install.packages("stringi", configure.args="--disable-pkg-config")

# sudo apt install libtesseract-dev libleptonica-dev tesseract-ocr-eng gnome-mahjongg librsvg2-dev libmagick++-dev
install.packages("tesseract", build = T,type = "source", upgrade = "always")
install.packages("rsvg", build = T,type = "source", upgrade = "always")
install.packages("gifski", build = T,type = "source", upgrade = "always")
install.packages("magick", build = T,type = "source", upgrade = "always")
install.packages("pdftools", build = T,type = "source", upgrade = "always")
install.packages("anndata", build = T,type = "source", upgrade = "always")

install.packages("future.callr", build = T,type = "source", upgrade = "always")
install.packages("RcppArmadillo", build = T,type = "source", upgrade = "always")

install.packages("RcppEigen", build = T,type = "source", upgrade = "always")

BiocManager::install('BiocParallel', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('DEsingle', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('DESeq2', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('edgeR', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('MAST', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('monocle', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('limma', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install("Rhdf5lib", dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F,force = TRUE)
BiocManager::install("rhdf5", dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('monocle', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
Sys.setenv(BPCELLS_DEBUG_INSTALL="true")
remotes::install_github("bnprks/BPCells", build = T,type = "source", upgrade = "always")

# remotes::install_github("satijalab/seurat", "seurat5", dependencies  = TRUE)
# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
remotes::install_github("satijalab/azimuth", quiet = TRUE, dependencies  = TRUE, build = T,type = "source", upgrade = "always")
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE, dependencies  = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install('biovizBase')
BiocManager::install('motifmatchr')
BiocManager::install('chromVAR')
BiocManager::install('ggseqlogo')


remotes::install_github("stuart-lab/signac", quiet = F, dependencies  = TRUE, build = T,type = "source", upgrade = "always")
install.packages("psychTools", dependencies  = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install("scran", dependencies  = TRUE, build = T, update = TRUE, ask = F, type = "source", upgrade = "always")

devtools::install_github("kassambara/ggpubr", dependencies  = TRUE, build = T,type = "source", upgrade = "always")

devtools::install_github('nghiavtr/BPSC', dependencies  = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install('Seurat', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
devtools::install_github('statOmics/zingeR', dependencies  = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install('SingleCellExperiment', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install('scater', dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
devtools::install_github('Zhangxf-ccnu/scDEA', dependencies  = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies  = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("haowulab/Wind", dependencies  = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("renozao/xbioc", dependencies  = TRUE, build = T,type = "source", upgrade = "always")

checkIfLibrariesInstalled <- function() {
  pkgList = c("archivist", "Biobase", "BiocGenerics", "BiocManager", "BiocParallel", 
              "BiocSingular", "callr", "colourpicker", "ComplexHeatmap", 
              "cowplot", "crayon", "debugme", "dendsort", "DESeq2",  "devtools", 
              "digest", "doParallel", "dplyr", "DT", "edgeR", "evaluate", "future",
              "GenomeInfoDb", "GenomicRanges", "ggalluvial", "ggnetwork", "ggplot2",
              "ggplotify", "ggpubr", "glue",  "gtools", "hdf5r",
              "heatmaply", "Hmisc", "hms", "igraph",  
              "IRanges", "irlba", "kableExtra", "knitr", "kohonen", 
              "limma", "magrittr", "manhattanly", "MASS", "Matrix", "mclust",
              "network", "orca", "parallel", "pdftools", 
              "pheatmap", "plotly", "plyr", "profvis", "pryr", "psychTools", 
              "RColorBrewer", "reactlog", "reactlog", "reshape2", "rintrojs", 
              "rmarkdown", "Rsomoclu", "Rtsne", "S4Vectors", "scater",
               "SCHNAPPs", "scran", "Seurat", "shiny", "shinyBS", 
              "shinycssloaders", "shinydashboard", "shinydashboardPlus", "shinyjqui",
              "shinyjs",  "shinytest", "shinythemes", "shinyTree",
              "shinyWidgets", "SIMLR", "SingleCellExperiment",  "spatstat", 
              "stringr", "SummarizedExperiment",  "threejs", "tibble", 
              "tidyr", "tidySingleCellExperiment", "tidyverse", "tools", "uwot"
              ) 
  

   
  missing = pkgList[!pkgList %in% installed.packages()]
  if(length(missing) == 0) {
    print("all should be good")
  } else {
    print(paste("The following packages are missing:" , missing))
    install.packages(missing, dependencies = NULL, build = T,type = "source", upgrade = "always")
  }
  for (pg in pkgList){
    library(pg,character.only = TRUE)
  }
}


checkIfLibrariesInstalled()
# new R version:
# 

devtools::install_github("cysouw/qlcMatrix")
devtools::install_github("Zhangxf-ccnu/scDEA")


biocList = c('BiocSingular', 'SingleR', "GSEABase", "GSVA","MAST", "monocle", "multtest", "rtracklayer", "scDEA", "SingleR", "Tempora", "Wind", "xbioc",
             'multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer', 
             'IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics',
             'DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment',
             'S4Vectors', 'BiocParallel', 'GSEABase', 'GSVA', "InteractiveComplexHeatmap"
             )
BiocManager::install(biocList, dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install("GSVA", dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)
BiocManager::install("InteractiveComplexHeatmap", dependencies  = TRUE, build = T,type = "source", upgrade = "always", update = TRUE, ask = F)

gitList = c('briatte/ggnetwork', 'mul118/shinyMCE',
            'RausellLab/CelliD',
            'BaderLab/Tempora', 'C3BI-pasteur-fr/UTechSCB-SCHNAPPs',
            'Albluca/distutils', 'Albluca/ElPiGraph.R'
            )
devtools::install_github(gitList, dependencies = T)
devtools::install_github('RausellLab/CelliD', dependencies = T)

# devtools::install_version("Matrix",version = "1.6.1.1")

gList = c('ggnetwork', 'shinyMCE',
            'CellID',
            'Tempora', 'SCHNAPPs',
            'distutils', 'ElPiGraph.R'
)

for (pg in gList) {
  require (`pg`,character.only = T)
}

BiocManager::install("GSEABase")
BiocManager::install("GSVA")

instList = c('BiocManager', 'pdftools',
             'shinycssloaders', 'network', 'igraph', 'mclust', 'shinyTree', 'shinydashboard', 'hdf5r',
             'forcats', 'kohonen', 'SCORPIUS', 'shinyBS', 'threejs', 'DT', 'shinythemes'
  
)
# setRepositories(ind = c(1,2,3))
# 
# devtools::install_version('spatstat', version = '1.64-1', repos = 'http://cran.us.r-project.org')
# update / reinstall macports https://trac.macports.org/wiki/Migration

install.packages("RJSONIO")

devtools::install_github("briatte/ggnetwork", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("mul118/shinyMCE",  dependencies = TRUE, build = T,type = "source", upgrade = "always")
# locally install destiny
# BiocManager::install("destiny")
# devtools::install_github("theislab/destiny", ref = "legacy", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("RausellLab/CelliD", ref = "legacy", dependencies = TRUE, build = T,type = "source", upgrade = "always")
# devtools::install_github("C3BI-pasteur-fr/TemporaFork", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("Albluca/distutils", dependencies = TRUE, build = T,type = "source", upgrade = "always") 
devtools::install_github("Albluca/ElPiGraph.R", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("rcannood/SCORPIUS", build_vignettes = TRUE, dependencies = TRUE, build = T,type = "source", upgrade = "always")

for (pg in instList) {
  require (`pg`,character.only = T)
}


for (pg in biocList) {
  require (`pg`,character.only = T)
}

for (pg in gList) {
  require (`pg`,character.only = T)
}

install.packages("igraph")


# the following package use openmp and require a special Makevars file

# CC=/usr/local/bin/gcc-13  #or path to clang-4.0 if that's your compiler..
# CXX=/usr/local/bin/g++-13
# CXX1X=/usr/local/bin/g++-13
# SHLIB_CXXLD=/usr/local/bin/g++-13
# FC=/usr/local/bin/gfortran-13 
# F77=/usr/local/bin/gfortran-13 
# MAKE=make -j2
# CXX11=/usr/local/bin/g++-13
# 
# SHLIB_OPENMP_CFLAGS=-fopenmp
# SHLIB_OPENMP_CXXFLAGS=-fopenmp
# SHLIB_OPENMP_FCFLAGS=-fopenmp
# SHLIB_OPENMP_FFLAGS=-fopenmp
# 
# ## -- compiling for OpenMP -> https://stackoverflow.com/a/5008957/271775
# PKG_CXXFLAGS=-I../inst/include -fopenmp
# ## -- linking for OpenMP
# ## PKG_LIBS= -fopenmp -lgomp
# PKG_LIBS= -L/usr/local/lib/gcc/ -lgomp `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"`
# LDFLAGS = -L/usr/local/lib/gcc/ -L/usr/local/lib/
#   CPPFLAGS+=-I/usr/local/include
# R_OPENMP_CFLAGS = -fopenmp
# R_OPENMP_FFLAGS = -fopenmp

# LDFLAGS="-L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++"
# echo 'export PATH="/usr/local/opt/llvm/bin:$PATH"' >> /Users/bernd/.bash_profile
# export LDFLAGS="-L/usr/local/opt/llvm/lib"
# export CPPFLAGS="-I/usr/local/opt/llvm/include"
# echo 'export PATH="/usr/local/opt/cython/bin:$PATH"' >> /Users/bernd/.bash_profile

openmpPacks = c("RcppEigen", "RcppML", "BH", "Rcpp", "RcppArmadillo", "RcppEigen", "RcppProgress",
                "RhpcBLASctl", "data.table", "dqrng", "magick", "mgcv", "Matrix", "Rtsne",
                "packrat", "renv", "sitmo")

install.packages(
  openmpPacks, 
  configure.vars = "CC=/usr/local/opt/llvm/bin/clang CXX=/usr/local/opt/llvm/bin/clang++",
  configure.args = "--with-openmp", build = T,type = "source", upgrade = "always",verbose = T
)

install.packages(openmpPacks, build = T,type = "source", upgrade = "always")
problems = c( "Rtsne")
install.packages(problems, build = T,type = "source", upgrade = "always",verbose = T,keep_outputs = T)

devtools::install_github("zdebruine/RcppML", build_vignettes = TRUE, dependencies = TRUE, build = T,type = "source", upgrade = "always", force = T) 
devtools::install_github("zdebruine/singlet", build_vignettes = TRUE, dependencies = TRUE, build = T,type = "source", upgrade = "always", force = T)

devtools::install_github("jkrijthe/Rtsne", build_vignettes = TRUE, dependencies = TRUE, build = T,type = "source", upgrade = "always", force = T,
                         configure.vars = "CC=/usr/local/opt/llvm/bin/clang CXX=/usr/local/opt/llvm/bin/clang++",
                         configure.args = "--with-openmp",verbose = T)
packs = c("archivist", "BiocParallel", "BiocSingular", "callr", "colourpicker", "ComplexHeatmap", 
          "cowplot", "crayon", "debugme", "digest", "dplyr", "DT", "edgeR", "evaluate", "ggalluvial", "ggplot2", 
          "ggpubr", "glue", "gtools", "heatmaply", "Hmisc", "hms", "rintrojs", "irlba", "kableExtra", "knitr",
          "kohonen", "magrittr", "manhattanly", "Matrix", "mclust", "orca", "parallel", "pdftools", "pheatmap",
          "plotly", "pryr", "psychTools", "RColorBrewer", "reactlog", "reshape2", "rintrojs", 
          "rmarkdown", "Rsomoclu", "Rtsne", "scater", "scDEA", "scran", "Seurat", "shiny", "shinyBS", 
          "shinydashboard", "shinydashboardPlus", "shinyjqui", "shinyjs", "shinyMCE", "shinytest", 
          "shinythemes", "shinyTree", "shinyWidgets", "SIMLR", "SingleCellExperiment", "SingleR", "stringr", 
          "threejs", "tibble", "tidyr", "tidySingleCellExperiment", "tidyverse", "uwot", "Wind", "xbioc", 
          "InteractiveComplexHeatmap", "plyr", "shinycssloaders", "profvis", "future", "tools", "doParallel", "reactlog", "MASS", "dendsort")

packs[!packs %in% installed.packages()]

