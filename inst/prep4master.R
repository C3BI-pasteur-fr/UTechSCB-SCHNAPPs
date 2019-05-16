## copy dev to prep4master
# git checkout prep4master
# git merge dev
getwd()
unlink("inst/develo", recursive = TRUE)
#### test
## copy to master
# git checkout master
# git merge prep4master