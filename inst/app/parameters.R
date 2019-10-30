
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
      fluidRow(div(tags$h3("Parameters for normalization to be used"), align = "center")),
      tags$p("SCHNAPPs uses normalized data for plots and calculations (unless stated otherwise). Here, the specific method can be set. rawNormalization means that no normalization will be performed."),
      tags$p("A table containing the first 20 cells and the normalized values is shown."),
      box(title = "Normalization method to use", solidHeader = TRUE, width = 12, status = 'primary', 
          fluidRow(
            column(width = 12,
                   radioButtons(
                     inputId = "normalizationRadioButton",
                     label = "choose a normalization method",
                     choices = normaliztionChoices,
                     selected = "DE_logNormalization",
                     width = "100%"
                   )
            )
          ),
          fluidRow(column(width = 10, 
                          verbatimTextOutput("normalizationRadioButtonValue"))),
          fluidRow(column(width = 12,
                          wellPanel(
                            # This outputs the dynamic UI component
                            uiOutput("normalizationsParametersDynamic")
                          ))),
          fluidRow(column(width = 12, offset = 1,
                          actionButton("updateNormalization", "apply changes", width = '80%', 
                                       style = "color: #fff; background-color: #A00272; border-color: #2e6da4")
          ))
      ),
      checkbsTT("normalizationRadioButton"),
      checkbsTT("normalizationRadioButtonValue"),
      checkbsTT("updateNormalization")
    ),
    box(
      title = "Sample of normalized values", solidHeader = TRUE, width = 12, status = 'primary',
      collapsible = TRUE, collapsed = TRUE,
      fluidRow(
        column( width = 12,
                tableSelectionUi("normalizationResult")
        )
      )
    )
    
  )
)
if (DEBUG) {
  cat(file = stderr(), paste("end: parameters.R\n"))
}
