#https://github.com/jcrodriguez1989/chatgpt.git

Sys.setenv(OPENAI_API_KEY = "")
library(chatgpt)
library(tidyverse)
Sys.setenv("OPENAI_MODEL" = "gpt-3.5-turbo")
Sys.setenv("OPENAI_MODEL" = "gpt-4-0314")
f = deparse(appendAnnotation) %>% paste(collapse = "\n")
addInfo = "the following function code requires two input variables. scEx is a SingleCellExpression object and annFile that is a file object of a csv file that has either the same row names as scEx or the row names are the column names of scEx. Include tests for incorrect input and multiple annotation file.\n\n"

chatgpt::comment_code(f)
chatgpt::create_unit_tests(paste(addInfo, f) )


install.packages("gptstudio")
install.packages("pak")
# Enable repository from jameshwade
options(repos = c(
  jameshwade = "https://jameshwade.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))
# Download and install gpttools in R
install.packages("gpttools")
library(gpttools)
# Browse the gpttools manual pages
help(package = "gpttools")
Sys.setenv(OPENAI_API_KEY = "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")


