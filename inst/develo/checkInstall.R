# new R version:
# 


biocList = c('BiocSingular', 'SingleR', 
             'multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer', 
             'IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics',
             'DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment',
             'S4Vectors', 'BiocParallel', 'GSEABase', 'GSVA'
             )
gitList = c('briatte/ggnetwork', 'mul118/shinyMCE',
            'RausellLab/CelliD',
            'C3BI-pasteur-fr/TemporaFork', 'C3BI-pasteur-fr/UTechSCB-SCHNAPPs',
            'Albluca/distutils', 'Albluca/ElPiGraph.R'
            )
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

