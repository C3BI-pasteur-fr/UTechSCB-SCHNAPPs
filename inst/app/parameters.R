
if (DEBUG) cat(file = stderr(), paste("parameters:", length(allTabs), " ", "\n"))

normaliztionChoices <- list(rawNormalization = "rawNormalization")
# parameterContributions = list()
# localContributionDir <- .SCHNAPPs_locContributionDir
parFiles <- dir(path = c(paste0(packagePath,  "/contributions"), localContributionDir), pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
for (fp in parFiles) {
  if (DEBUG) {
    cat(file = stderr(), paste(fp, "\n"))
  }

  myNormalizationChoices <- c()
  source(fp, local = TRUE)
  if (length(myNormalizationChoices) > 0) {
    for (li in 1:length(myNormalizationChoices)) {
      liVal <- myNormalizationChoices[[li]]
      if (length(liVal) > 0) {
        # if (DEBUG) cat(file = stderr(), paste("normalization Choice: ", liVal, "\n"))
        oldNames <- names(normaliztionChoices)
        normaliztionChoices[[length(normaliztionChoices) + 1]] <- liVal
        names(normaliztionChoices) <- c(oldNames, names(myNormalizationChoices)[li])
      }
    }
  }
  if (DEBUG) {
    cat(file = stderr(), paste("end:", fp, "\n"))
    cat(file = stderr(), paste("end:", normaliztionChoices, "\n"))
  }
}

# here we add content to the page on the rigth (main visualization window)
allTabs[[length(allTabs) + 1]] <- list(
  shinydashboard::tabItem(
    "normalizations",
    list(
      tags$h3("Parameters for normalization to be used"),
      tags$p("SCHNAPPs uses normalized data for plots and calculations (unless stated otherwise). Here, the specific method can be set. rawNormalization means that no normalization will be performed."),
      tags$p("A table containing the first 20 cells and the normalized values is shown."),
      fluidRow(column(
        10,
        radioButtons(
          inputId = "normalizationRadioButton",
          label = "Normalization to use",
          choices = normaliztionChoices,
          selected = "DE_logNormalization",
          width = "100%"
        )
       )),
      fluidRow(column(10, verbatimTextOutput("normalizationRadioButtonValue"))),
      wellPanel(
        # This outputs the dynamic UI component
        uiOutput("normalizationsParametersDynamic")
      )
    ),
    tableSelectionUi("normalizationResult")
  )
)
if (DEBUG) {
  cat(file = stderr(), paste("end: parameters.R\n"))
}
