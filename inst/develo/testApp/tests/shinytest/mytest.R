app <- ShinyDriver$new("/Users/bernd/Rstudio/UTechSCB-SCHNAPPs/inst/develo/testApp/", loadTimeout = 1e+06, seed = 1,debug = "all")
step = 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshotInit("mytest")
Sys.sleep(time = 60)
app$setInputs(sideBarID = "parameters")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sidebarItemExpanded = "GeneralQC")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "gQC_umiHist")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "gQC_sampleHist")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "gQC_variancePC")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "DE_scaterQC")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sidebarItemExpanded = "Co-expression")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "coexpressionAll")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "coexpressionSelected")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()

#not working?
# app$setInputs(sideBarID = "CoExpressionViolin")
# Sys.sleep(time = 10)
# step = step + 1
# cat(file = stderr(), paste("step ", step), "\n")
# app$snapshot()

app$setInputs(sideBarID = "alluvialTab", wait_ = FALSE, values_ = FALSE)
Sys.sleep(time = 120)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
# app$setInputs(sidebarItemExpanded = "DataExploration")
# Sys.sleep(time = 10)
# step = step + 1
# cat(file = stderr()paste(, "step ", st)ep, "\n")
# app$snapshot()
# app$setInputs(sideBarID = "DE_expression")
# Sys.sleep(time = 10)
step = step + 1
# cat(file = stderr()paste(, "step ", st)ep, "\n")
# app$snapshot()
app$setInputs(sideBarID = "DE_panelPlot", wait_ = FALSE, values_ = FALSE)
Sys.sleep(time = 60)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sidebarItemExpanded = "Subclusteranalysis")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "sCA_dge")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(updateDGEParameters = "click")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "parameters")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "Intro")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$snapshot(list(output = "introRMD"))
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "parameters")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sidebarItemExpanded = "GeneralQC")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "gQC_umiHist")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$snapshot(list(output = "gQC_plotUmiHist"))
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$setInputs(sideBarID = "gQC_sampleHist")
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
app$snapshot(list(output = "gQC_plotSampleHist"))
app$setInputs(sideBarID = "gQC_variancePC")
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot(list(output = "gQC_variancePCA"))
app$setInputs(sideBarID = "DE_scaterQC")
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot(list(output = "DE_scaterQC"))
app$setInputs(sidebarItemExpanded = "Co-expression")
app$setInputs(sideBarID = "coexpressionAll")
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot(list(output = "coExpHeatmapModule-pHeatMapPlot"))
app$setInputs(sideBarID = "coexpressionSelected")
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot(list(output = "coE_selected-clusterPlot"))
Sys.sleep(time = 10)
step = step + 1
cat(file = stderr(), paste("step ", step), "\n")
app$snapshot()
