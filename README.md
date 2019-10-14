<img src="inst/www/images/schnappsLogo.png" align="right" alt="" width="120" />

# SCHNAPPs - Single Cell sHiNy APP(s)

## Overview

Shiny app for the exploration and analysis of single cell RNAseq data as it comes from 10X or MARSseq technologies. It is currently being developed based on user requests of the Cytometry and Biomarkers UTechS at the Pasteur Institute, Paris. The goal is to enable the users of our platform to explore their data, select cells they would like to work with and then perform the final analysis together with the bioinformatics support at Pasteur.


## Installation

```
update.packages()
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("mul118/shinyMCE")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# update bioconductor packages if required
BiocManager::install()
BiocManager::install("BiocSingular")
devtools::install_github("C3BI-pasteur-fr/UTechSCB-SCHNAPPs")
```

### history functionality

To take advantage of the history functionality orca needs to be installed:
(https://github.com/plotly/orca#installation)

orca is part of ploty so nothing to be done for R. But the pdftools are required:

```
install.packages("pdftools")

```

Notes:

1. once history check box is checked only the following plots will be recorded, the current plot will not be saved.
2. any selection in the plot will not be visible. It is not a screen-shot.
3. The date of recording is added to the title.
4. Creating of the history will take some time.
5. once activated, any changes of parameters will trigger re-plot and also a save to the history. This "feature" can also be used to create of the current plot by adding a "," in one of the fields.

## create sample data set

Load a small set of 200 PBMC cells and save to a file in the local directory. This file can be uploaded using the app.

```
data("scEx", package = "SCHNAPPs")
save(file = "scEx.Rdata", list = "scEx")
```


## Running schnapps

To start the app:

```
library(SCHNAPPs)
schnapps()
```

### history functionality

Plots can be automatically stored in a PDF file as they are created. This allows to somehow track what is being done. To enable this, the parameter historyFile has to be set to a file (probably non-existing, otherwise it will be appended to).


### generate data files

A singleCellExperiment object is required, saved in a file RData object using 

```
save(file = "filename.RData", "singleCellExperiementObject")
```



Please see [GitHub](https://c3bi-pasteur-fr.github.io/UTechSCB-SCHNAPPs/) for further documentation on how to use schnapps.

Please see [GitHub Contributions](https://github.com/baj12/SCHNAPPsContributions) for additional tools, not directly part of the SCHNAPPs pacakage.



## Extending SCHAPPs

See [SCHNAPPsContribution](https://github.com/baj12/SCHNAPPsContributions) for examples and dummy contributions on how add functionality.

## Credits

The original version of this app is based on CellView (https://github.com/mohanbolisetty/CellView), but it was substantially modified. It helped me get started.

We are also greatful to all the members of the single cell working group at Pasteur Paris.
