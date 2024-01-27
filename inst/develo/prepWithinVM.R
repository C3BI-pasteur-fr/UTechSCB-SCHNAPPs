# VBoxManage modifymedium kafkaworkshoposajava-disk1.vdi –resize 20000
# disablee apport:
# sudo vii /etc/default/apport
# set enabled=0


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
devtools::install_github("Zhangxf-ccnu/scDEA")
BiocManager::install("glmGamPoi")

BiocManager::install("InteractiveComplexHeatmap")
install.packages("bookdown")
install.packages("future.callr")
BiocManager::install ("GSVA")
BiocManager::install("GSEABase")

# copy /home/bernd/Rstudio/SCORPIUS to data directory and install from within VM


devtools::install_github("BaderLab/Tempora")

BiocManager::install("CelliD")

devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs", dependencies = TRUE)
