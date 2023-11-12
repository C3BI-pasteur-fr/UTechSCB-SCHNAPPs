# not sure this is needed
# BiocManager::install("diffcyt")
# 
# BiocManager::install("HDCytoData")
# BiocManager::install("CATALYST")
# 
# devtools::install_github("hrbrmstr/dtupdate")
# library(dtupdate)
# dtupdate::github_update()

install.packages("RcppArmadillo")
install.packages("RcppEigen")
remotes::install_github("satijalab/seurat", "seurat5", dependencies  = TRUE)
install.packages("psychTools")
BiocManager::install("scran", dependencies  = TRUE)
devtools::install_github("kassambara/ggpubr", dependencies  = TRUE)

devtools::install_github('nghiavtr/BPSC')
BiocManager::install('DEsingle')
devtools::install_github('nghiavtr/BPSC')
BiocManager::install('DESeq2')
BiocManager::install('edgeR')
BiocManager::install('MAST')
BiocManager::install('monocle')
BiocManager::install('limma')
BiocManager::install('Seurat')
devtools::install_github('statOmics/zingeR')
BiocManager::install('SingleCellExperiment')
BiocManager::install('scater')
devtools::install_github('Zhangxf-ccnu/scDEA')


checkIfLibrariesInstalled <- function() {
  pkgList = c("archivist", "Biobase", "BiocGenerics", "BiocManager", "BiocParallel", 
              "BiocSingular", "callr", "colourpicker", "ComplexHeatmap", 
              "cowplot", "crayon", "debugme", "dendsort", "DESeq2",  "devtools", 
              "digest", "doParallel", "dplyr", "DT", "edgeR", "evaluate", "future",
              "GenomeInfoDb", "GenomicRanges", "ggalluvial", "ggnetwork", "ggplot2",
              "ggplotify", "ggpubr", "glue", "GSEABase", "GSVA", "gtools", "hdf5r",
              "heatmaply", "Hmisc", "hms", "igraph",  
              "IRanges", "irlba", "kableExtra", "knitr", "kohonen", 
              "limma", "magrittr", "manhattanly", "MASS", "MAST", "Matrix", "mclust",
              "monocle", "multtest", "network", "orca", "parallel", "pdftools", 
              "pheatmap", "plotly", "plyr", "profvis", "pryr", "psychTools", 
              "RColorBrewer", "reactlog", "reactlog", "reshape2", "rintrojs", 
              "rmarkdown", "Rsomoclu", "rtracklayer", "Rtsne", "S4Vectors", "scater",
              "scDEA", "SCHNAPPs", "scran", "Seurat", "shiny", "shinyBS", 
              "shinycssloaders", "shinydashboard", "shinydashboardPlus", "shinyjqui",
              "shinyjs",  "shinytest", "shinythemes", "shinyTree",
              "shinyWidgets", "SIMLR", "SingleCellExperiment", "SingleR", "spatstat", 
              "stringr", "SummarizedExperiment", "Tempora", "threejs", "tibble", 
              "tidyr", "tidySingleCellExperiment", "tidyverse", "tools", "uwot", 
              "Wind", "xbioc") 
  

  
  missing = pkgList[!pkgList %in% installed.packages()]
  if(length(missing) == 0) {
    print("all should be good")
  } else {
    print(paste("The following packages are missing:" , missing))
    install.packages(missing)
  }
  for (pg in pkgList){
    library(pg,character.only = TRUE)
  }
}


checkIfLibrariesInstalled()
# new R version:
# 

devtools::install_github("Zhangxf-ccnu/scDEA")


biocList = c('BiocSingular', 'SingleR', 
             'multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer', 
             'IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics',
             'DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment',
             'S4Vectors', 'BiocParallel', 'GSEABase', 'GSVA', "InteractiveComplexHeatmap"
             )
BiocManager::install(biocList)
gitList = c('briatte/ggnetwork', 'mul118/shinyMCE',
            'RausellLab/CelliD',
            'C3BI-pasteur-fr/TemporaFork', 'C3BI-pasteur-fr/UTechSCB-SCHNAPPs',
            'Albluca/distutils', 'Albluca/ElPiGraph.R'
            )
devtools::install_github(gitList, dependencies = T)

devtools::install_version("Matrix",version = "1.6.1.1")

gList = c('ggnetwork', 'shinyMCE',
            'CellID',
            'Tempora', 'SCHNAPPs',
            'distutils', 'ElPiGraph.R'
)


instList = c('BiocManager', 'pdftools',
             'shinycssloaders', 'network', 'igraph', 'mclust', 'shinyTree', 'shinydashboard', 'hdf5r',
             'forcats', 'kohonen', 'SCORPIUS', 'shinyBS', 'threejs', 'DT', 'shinythemes'
  
)
# setRepositories(ind = c(1,2,3))
# 
# devtools::install_version('spatstat', version = '1.64-1', repos = 'http://cran.us.r-project.org')
# update / reinstall macports https://trac.macports.org/wiki/Migration


devtools::install_github("briatte/ggnetwork", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("mul118/shinyMCE",  dependencies = TRUE, build = T,type = "source", upgrade = "always")
# locally install destiny
# BiocManager::install("destiny")
# devtools::install_github("theislab/destiny", ref = "legacy", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("RausellLab/CelliD", ref = "legacy", dependencies = TRUE, build = T,type = "source", upgrade = "always")
BiocManager::install("GSVA")
devtools::install_github("C3BI-pasteur-fr/TemporaFork", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE, build = T,type = "source", upgrade = "always")
devtools::install_github("Albluca/distutils", dependencies = TRUE, build = T,type = "source", upgrade = "always") 
devtools::install_github("Albluca/ElPiGraph.R", dependencies = TRUE, build = T,type = "source", upgrade = "always")

for (pg in instList) {
  require (`pg`,character.only = T)
}


for (pg in biocList) {
  require (`pg`,character.only = T)
}

for (pg in gList) {
  require (`pg`,character.only = T)
}

devtools::install_github('chris-mcginnis-ucsf/DoubletFinder', dependencies = TRUE, build = T,type = "source", upgrade = "always")

install.packages("igraph")
