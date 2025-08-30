######################################    ERRORS    ######################################
test_that("errors are returned if the object is not of class DEprot or DEprot.analyses", {
  expect_error(perform.PCA(DEprot.object = DEprot::sample.config))
})


##########################################################################################

test_that("works for data of class DEprot on unimputed data", {
  expect_no_error(perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp, which.data = "normalized"))
})


test_that("works for data of class DEprot on a subset of samples", {
  expect_no_error(perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp, which.data = "normalized", sample.subset = DEprot::test.toolbox$dpo.norm@metadata$column.id[1:8]))
})
