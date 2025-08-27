##########################################################################################

test_that("no errors are returned if object of class DEprot.analyses is provided", {
  expect_no_error(compare.ranking(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})
