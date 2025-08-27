test_that("errors are returned if the object is not of class DEprot", {
  expect_error(compare.imp.methods(DEprot.object = "not DEprot object", percentage.test = 100, sample.group.column = "combined.id"))
})


##########################################################################################

test_that("no error is returned when comparing the imputations", {
  expect_no_failure(suppressWarnings(compare.imp.methods(DEprot.object = DEprot::test.toolbox$dpo.norm, percentage.test = 100, sample.group.column = "combined.id", missForest.cores = 1)))
})
