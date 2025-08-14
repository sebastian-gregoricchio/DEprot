######################################    ERRORS    ######################################
test_that("errors are returned if the object is not of class DEprot.analyses", {
  expect_error(check.pvalues(DEprot::sample.config))
})


test_that("errors are returned if the contrast indicated is not available", {
  expect_error(check.pvalues(DEprot::test.toolbox$diff.exp.limma,contrast = 100))
})


##########################################################################################

test_that("errors are returned if the object is not of class DEprot", {
  expect_no_error(check.pvalues(DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})
