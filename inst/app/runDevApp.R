#' this is used to run the app without installing it.
#'
#'
library(reactlog)
# if (!exists(".schnappsEnv")) {
.schnappsEnv <- new.env(parent=emptyenv())
# }

localContributionDir = "~/Rstudio/SCHNAPPsContributions-dont/"
# localContributionDir = ""
defaultValueSingleGene = "itgae" # CD52
defaultValueMultiGenes = "CD52, S100A9, S100A4" # itgae, cd69, itga1" # CD52, S100A9, S100A4
# defaultValueMultiGenes = "prf1, Gzmb, IFNG, PDCD1, HAVCR2, LAG3, TSC22D3,ZFP36L2"
defaultValueRegExGene = "" # tip: '^CD7$|^KIT$; genes with min expression
DEBUG = T
DEBUGSAVE = F
historyPath = "~/Rstudio/Schnapps/history"
historyPath = NULL

assign(".SCHNAPPs_locContributionDir", localContributionDir, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign(".SCHNAPPs_defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUG", DEBUG, envir = .schnappsEnv)
assign(".SCHNAPPs_DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("localContributionDir", localContributionDir, envir = .schnappsEnv)
assign("defaultValueSingleGene", defaultValueSingleGene, envir = .schnappsEnv)
assign("defaultValueMultiGenes", defaultValueMultiGenes, envir = .schnappsEnv)
assign("defaultValueRegExGene", defaultValueRegExGene, envir = .schnappsEnv)
assign("DEBUG", DEBUG, envir = .schnappsEnv)
assign("DEBUGSAVE", DEBUGSAVE, envir = .schnappsEnv)
assign("historyPath", historyPath, envir = .schnappsEnv)
ls(.schnappsEnv)


# # Scran parameters
# defaultValues = list()
# defaultValues[["coEtgMinExpr"]] = 100
# defaultValues[["pcaScale"]] = FALSE
# defaultValues[["cellsFiltersOut"]] = "D2ex_54, D2ex_64, D2ex_93, D2ex_95, D2ex_96, D3en1_8, D3en1_9, D3en1_12, D3en1_14, D3en1_15, D3en1_18, D3en1_19, D3en1_20, D3en1_21, D3en1_22, D3en1_23, D3en1_25, D3en1_31, D3en1_32, D3en1_34, D3en1_35, D3en1_36, D3en1_39, D3en1_43, D3en1_44, D3en1_47, D3en1_51, D3en1_52, D3en1_53, D3en1_57, D3en1_58, D3en1_59, D3en1_61, D3en1_62, D3en1_64, D3en1_66, D3en1_70, D3en1_72, D3en1_74, D3en1_76, D3en1_77, D3en1_78, D3en1_80, D3en1_82, D3en1_86, D3en1_91, D3en1_92, D3en1_95, D3en1_96, D3en2_1, D3en2_4, D3en2_8, D3en2_9, D3en2_10, D3en2_11, D3en2_14, D3en2_15, D3en2_18, D3en2_19, D3en2_20, D3en2_28, D3en2_29, D3en2_34, D3en2_35, D3en2_36, D3en2_40, D3en2_41, D3en2_43, D3en2_44, D3en2_46, D3en2_47, D3en2_49, D3en2_53, D3en2_54, D3en2_59, D3en2_65, D3en2_67, D3en2_69, D3en2_70, D3en2_75, D3en2_76, D3en2_77, D3en2_78, D3en2_79, D3en2_80, D3en2_82, D3en2_84, D3en2_86, D3en2_87, D3en2_88, D3en2_89, D3en2_95, D3en2_96, D3en3_3, D3en3_4, D3en3_12, D3en3_13, D3en3_14, D3en3_15, D3en3_17, D3en3_19, D3en3_20, D3en3_22, D3en3_23, D3en3_24, D3en3_25, D3en3_26, D3en3_31, D3en3_33, D3en3_35, D3en3_37, D3en3_38, D3en3_42, D3en3_43, D3en3_45, D3en3_47, D3en3_49, D3en3_56, D3en3_57, D3en3_58, D3en3_60, D3en3_61, D3en3_62, D3en3_63, D3en3_66, D3en3_69, D3en3_70, D3en3_71, D3en3_76, D3en3_79, D3en3_88, D3en3_90, D3en3_91, D3en3_93, D3en3_96, D3en4_1, D3en4_5, D3en4_10, D3en4_11, D3en4_13, D3en4_15, D3en4_17, D3en4_18, D3en4_21, D3en4_24, D3en4_27, D3en4_28, D3en4_29, D3en4_34, D3en4_36, D3en4_37, D3en4_40, D3en4_44, D3en4_45, D3en4_49, D3en4_50, D3en4_51, D3en4_52, D3en4_53, D3en4_55, D3en4_57, D3en4_59, D3en4_62, D3en4_66, D3en4_67, D3en4_70, D3en4_75, D3en4_76, D3en4_77, D3en4_78, D3en4_79, D3en4_80, D3en4_88, D3en4_93, D3en4_94, D3en4_95, D3ex_1, D3ex_2, D3ex_7, D3ex_8, D3ex_10, D3ex_11, D3ex_13, D3ex_14, D3ex_15, D3ex_17, D3ex_18, D3ex_19, D3ex_20, D3ex_21, D3ex_22, D3ex_23, D3ex_25, D3ex_26, D3ex_27, D3ex_28, D3ex_29, D3ex_30, D3ex_31, D3ex_32, D3ex_33, D3ex_35, D3ex_37, D3ex_38, D3ex_42, D3ex_43, D3ex_44, D3ex_45, D3ex_46, D3ex_47, D3ex_49, D3ex_50, D3ex_52, D3ex_53, D3ex_56, D3ex_57, D3ex_58, D3ex_59, D3ex_60, D3ex_61, D3ex_62, D3ex_63, D3ex_64, D3ex_65, D3ex_67, D3ex_69, D3ex_70, D3ex_72, D3ex_73, D3ex_75, D3ex_78, D3ex_79, D3ex_80, D3ex_81, D3ex_82, D3ex_84, D3ex_85, D3ex_86, D3ex_89, D3ex_93, D3ex_94, D3ex_95, D3ex_96, D74_4, D74_11, D74_19, D74_30, D74_38, D74_50, D74_72, D74_80, D71_13, D71_24, D71_25, D71_46, D71_50, D71_74, D71_77, D71_96, D72_6, D72_17, D72_26, D72_51, D72_56, D72_60, D73_53, D73_66, D73_74, D73_75, D73_77, D73_80, D73_88, D10631_1, D10631_6, D10631_13, D10631_16, D10631_28, D10631_30, D10631_31, D10631_32, D10631_35, D10631_49, D10631_59, D10631_77, D10631_86, D101_2, D101_3, D101_4, D101_6, D101_9, D101_12, D101_15, D101_18, D101_19, D101_20, D101_23, D101_26, D101_30, D101_32, D101_33, D101_38, D101_39, D101_40, D101_41, D101_47, D101_49, D101_54, D101_55, D101_60, D101_61, D101_62, D101_63, D101_64, D101_65, D101_67, D101_68, D101_70, D101_71, D101_75, D101_76, D101_78, D101_79, D101_80, D101_85, D101_86, D101_87, D101_90, D101_92, D101_94, D101_96, D102_2, D102_5, D102_8, D102_9, D102_10, D102_12, D102_13, D102_15, D102_16, D102_19, D102_20, D102_22, D102_23, D102_27, D102_28, D102_29, D102_30, D102_31, D102_32, D102_33, D102_34, D102_35, D102_36, D102_37, D102_38, D102_41, D102_42, D102_43, D102_44, D102_46, D102_48, D102_49, D102_52, D102_53, D102_56, D102_60, D102_62, D102_65, D102_67, D102_70, D102_72, D102_74, D102_77, D102_84, D102_85, D102_88, D102_89, D102_90, D102_92, D102_95, D17All1_13, D17All1_94, D17All2_5, D17All2_50, D17All2_89, D17TGFB_60, D17TGFB_62, D2ex_94, D3en1_40, D3en2_51, D3en3_18, D3en3_40, D3en3_77, D3en3_81, D3en3_83, D3en4_90, D3en4_96, D74_48, D74_84, D74_96, D71_89, D72_94, D73_47, D10631_47, D10631_75, D10631_78, D10631_79, D10631_80, D10631_83, D10631_84, D10631_94, D101_8, D101_11, D101_46, D102_47, D102_80, D102_81, D102_83, D172444_94, D172444_95, D172444_96, D17All1_14, D17All1_15, D17All1_17, D17All1_85, D17All2_96, D17TGFB_96"
# defaultValues[["selectIds"]] = ""
# defaultValues[["sampleInput"]] =FALSE
# defaultValues[["whichscLog"]] = "calcLog"
# defaultValues[["normalizationRadioButton"]] = "DE_scaterNormalization"
# defaultValues[["pcaN"]] = 951
# defaultValues[["tabsetCluster"]] = "snnGraph"
# defaultValues[["minGenesGS"]] = 1
# defaultValues[["pcaRank"]] = 14
# defaultValues[["minGenes"]] = 1
# defaultValues[["alluiv1"]] = "cluster"
# defaultValues[["alluiv2"]] = "dbCluster"
# defaultValues[["pcaCenter"]] =FALSE


# Seurat parameters
defaultValues = list()
defaultValues[["cellsFiltersOut"]] = "AAAGATCTGGGCAA-1, AAAGCAGAAGCCAT-1, AACGCCCTGCTTAG-1, AAGGTCTGGTATGC-1, AATGTAACGTTTGG-1, AATTACGAGTAGCT-1, ACACAGACACCTGA-1, ACATGGTGCGTTGA-1, ACCTGGCTGTCTTT-1, ACTTAAGACCACAA-1, ACTTGTACCCGAAT-1, ACTTTGTGCGATAC-1, AGAGGTCTACAGCT-1, ATCACGGATTGCTT-1, ATTACCTGGGCATT-1, CACGCTACTTGACG-1, CAGTGTGAACACGT-1, CCAATGGAACAGCT-1, CCAGTCTGCGGAGA-1, CGACCTTGGCAAGG-1, CGAGCCGACGACAT-1, CGGAATTGCACTAG-1, CGTAACGAATCAGC-1, CGTACCACGCTACA-1, CGTACCTGGACGAG-1, CTAGTTTGAGTACC-1, CTCAGCTGTTTCTG-1, CTCATTGATTGCTT-1, CTGGCACTGGACAG-1, CTTAACACGAGCTT-1, CTTAAGCTTCCTCG-1, GAAAGATGTTTGCT-1, GAACGTTGACGGAG-1, GAATGGCTAAGATG-1, GACCATGACTCTCG-1, GACTGAACAACCGT-1, GCCACTACCTACTT-1, GCGAAGGAGAGCTT-1, GCTACAGATCTTAC-1, GGCACGTGTGAGAA-1, GTCAACGATCAGGT-1, GTGAACACAGATCC-1, GTGTCAGAATGCTG-1, GTTAAAACTTCGCC-1, TAAGATACCCACAA-1, TACGCAGACGTCTC-1, TACGCGCTCTTCTA-1, TACGGCCTGTCCTC-1, TATCACTGACTGTG-1, TCCCGATGCTGTGA-1, TCGCACACCATCAG-1, TCGTGAGAACTGTG-1, TGAAGCTGAGACTC-1, TGAGACACTGTGCA-1, TGAGCTGAGCGAGA-1, TGGAGACTGAAACA-1, TGGATGTGATGTCG-1, TGGCAATGGAGGGT-1, TGGTCAGACCGTTC-1, TGTTAAGATTGGCA-1, TTACTCGAACGTTG-1, TTCAAGCTTCCAAG-1"
defaultValues[["selectIds"]] = ""
defaultValues[["pcaN"]] = 2000
defaultValues[["pcaScale"]] = TRUE
defaultValues[["sampleInput"]] =FALSE
defaultValues[["whichscLog"]] = "calcLog"
defaultValues[["normalizationRadioButton"]] = "DE_seuratLogNorm"
defaultValues[["hvgSelection"]] = "vst"
defaultValues[["alluiv1"]] = "seurartCluster"
defaultValues[["alluiv2"]] = "dbCluster"
defaultValues[["tabsetCluster"]] = "seurat_Clustering"
defaultValues[["minGenesGS"]] = 1
defaultValues[["minGenes"]] = 1
defaultValues[["seurClustDims"]] = 10
defaultValues[["seurClustk.param"]] = 20
defaultValues[["useSeuratPCA"]] = TRUE


assign("defaultValues", defaultValues, envir = .schnappsEnv)




devscShinyApp = TRUE
packagePath <<- "inst/app"
source(paste0(packagePath,  "/ui.R"))
source(paste0(packagePath,  "/server.R"))

app <- shinyApp(ui = scShinyUI, server = scShinyServer, enableBookmarking = "server")
options(shiny.reactlog=TRUE)
runApp(app)

# schnapps(Ï
# defaultValueMultiGenes = "IL7R, CCR7,CD14, LYZ ,IL7R, S100A4,MS4A1 ,CD8A,FCGR3A, MS4A7 ,GNLY, NKG7,FCER1A, CST3,PPBP",
# defaultValueSingleGene = "MS4A1", DEBUG=TRUE
# )


# sctkEx = SCtkExperiment(assays=list(counts=as.matrix(assays(scEx)[['counts']]), 
#                                     logcounts = as.matrix(assays(scEx)[['logcounts']])),
#                         colData = colData(scEx),
#                         rowData = rowData(scEx))
# singleCellTK(inSCE = sctkEx)




