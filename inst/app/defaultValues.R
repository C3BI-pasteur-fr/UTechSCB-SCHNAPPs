# defaultValueSingleGene='PF3D7_0413000, PF3D7_0420800, PF3D7_0421400, PF3D7_0711600, PF3D7_0711800, PF3D7_0712100, PF3D7_0809000,
# PF3D7_1240500, PF3D7_1241000, PF3D7_0412500, PF3D7_0412800, PF3D7_0421000, PF3D7_0712700, PF3D7_0808500, PF3D7_1240800'
# defaultValueMultiGenes='PF3D7_0413000, PF3D7_0420800, PF3D7_0421400, PF3D7_0711600, PF3D7_0711800, PF3D7_0712100, PF3D7_0809000,
# PF3D7_1240500, PF3D7_1241000, PF3D7_0412500, PF3D7_0412800, PF3D7_0421000, PF3D7_0712700, PF3D7_0808500, PF3D7_1240800'
allowedColors <- unique(c(
   "#2D96FA", "#8D8889", "#E8E0E0", "#FF2F7D", "#60A948", "#732BA3", "#FFC720", "#E96B0C", "#22ABA4", 
   "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e","#000000", "#c51b7d", "#e9a3c9", "#fde0ef", 
   "#e0e0e0", "#FAF4F5", "#999999", "#414144", "#E3E9EB", "#D4070F", "#1C1AAF", "#4d4d4d", "#e6f5d0", "#a1d76a", "#4d9221", "#762a83","white", 
   "#af8dc3", "#e7d4e8", "#d9f0d3", "#7fbf7b", "#1b7837", "#b35806", "#f1a340", "#fee0b6", "#d8daeb", "#998ec3", 
   "#542788",  "#fddbc7", "#d1e5f0", "#ef8a62", "#b2182b", "#67a9cf", "#2166ac"))
names(allowedColors) = make.names(1:length(allowedColors))

# .schnappsEnv$defaultValues[["minGenesGS"]] = 100
# defaultValueSingleGene <- "CD52"
# defaultValueMultiGenes <- "CD52, S100A4, S100A9, S100A8"
# defaultValueRegExGene <- "" # tip: '^CD7$|^KIT$; genes with min expression
# DEBUG <- TRUE
# DEBUGSAVE <- FALSE
