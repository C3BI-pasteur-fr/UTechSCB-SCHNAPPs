options(shiny.sanitize.errors = FALSE)
library(SCHNAPPs)
#
#schnappsLite(data = "epdc.rn-sham2.v2.lite.RData", DEBUG = T, historyPath = "history")
#schnappsLite(data = "rnAllEPDC.lite.RData", DEBUG = T, historyPath = "history")
defaultValueMultiGenes = "wt1,pdgfrb,agtr1a,col3a1,col1a1,postn,tbx18,scx,npr1,vegfa,npr2,notch2,tagln,acta2,ctgf, Rack1, Bmp4,  Eef2, col3a1,col1a1,tgfb3,postn"
defaultValueSingleGene = "wt1"


defaultValues = list()
defaultValues[["coEtgMinExpr"]] = 100


schnapps(defaultValues = defaultValues, 
         localContributionDir = "/SCHNAPPsContributions/", DEBUG = T, historyPath = "history", port = 3838)

