# cellchat
# 
devtools::install_github("sqjin/CellChat")

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

showDatabaseCategory(CellChatDB)
