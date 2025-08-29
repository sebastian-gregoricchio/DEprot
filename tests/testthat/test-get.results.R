test_that("get.results works for DEprot.analyses objects", {
  expect_no_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, 1))
})


test_that("get.results does not work for object not of class DEprot.analyses", {
  expect_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$dpo.imp, 1))
})

test_that("get.results does not work if the contrast requested is not available", {
  expect_error(get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, 100))
})
