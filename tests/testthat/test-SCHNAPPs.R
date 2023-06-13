test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
.schnappsEnv$DEBUG = FALSE
sfFile = system.file("app/","serverFunctions.R", package = "SCHNAPPs" )
source(sfFile)

test_that("printTimeEnd", {
  expect_snapshot(printTimeEnd(start.time = as.POSIXct("2001-03-20 14:12:22 GMT"),messtr = "test"))
})
