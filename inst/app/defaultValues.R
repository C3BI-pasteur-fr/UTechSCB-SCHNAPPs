# defaultValueSingleGene='PF3D7_0413000, PF3D7_0420800, PF3D7_0421400, PF3D7_0711600, PF3D7_0711800, PF3D7_0712100, PF3D7_0809000,
# PF3D7_1240500, PF3D7_1241000, PF3D7_0412500, PF3D7_0412800, PF3D7_0421000, PF3D7_0712700, PF3D7_0808500, PF3D7_1240800'
# defaultValueMultiGenes='PF3D7_0413000, PF3D7_0420800, PF3D7_0421400, PF3D7_0711600, PF3D7_0711800, PF3D7_0712100, PF3D7_0809000,
# PF3D7_1240500, PF3D7_1241000, PF3D7_0412500, PF3D7_0412800, PF3D7_0421000, PF3D7_0712700, PF3D7_0808500, PF3D7_1240800'
allowedColors <- unique(c("#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e","#c51b7d","#e9a3c9",
                         "#fde0ef","#e6f5d0","#a1d76a","#4d9221","#762a83","#af8dc3","#e7d4e8","#d9f0d3",
                         "#7fbf7b","#1b7837","#b35806","#f1a340","#fee0b6","#d8daeb","#998ec3","#542788",
                         "#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac","#b2182b","#ef8a62",
                         "#fddbc7","#e0e0e0","#999999","#4d4d4d"))

defaultValueSingleGene <- "CD52"
defaultValueMultiGenes <- "CD52, S100A4, S100A9, S100A8"
defaultValueRegExGene <- "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG <- TRUE
DEBUGSAVE <- FALSE
