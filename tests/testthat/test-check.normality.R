######################################    ERRORS    ######################################
test_that("errors are returned if the object is not of class DEprot", {
  expect_error(check.normality(DEprot.object = DEprot::sample.config))
})


test_that("errors are returned if the data required are not correct", {
  expect_error(check.normality(DEprot.object = DEprot::test.toolbox$dpo.imp, which.data = "wrong data"))
})


##########################################################################################

test_that("no errors are returned if the object class DEprot", {
  expect_no_error(check.normality(DEprot.object = DEprot::test.toolbox$dpo.imp))
})
