######################################    ERRORS    ######################################
test_that("errors are returned if the object is not of class DEprot or DEprot.analyses", {
  expect_error(SAINTq(DEprot.object = DEprot::sample.config))
})

test_that("errors are returned if the object is not provided", {
  expect_error(SAINTq())
})

test_that("errors are returned if the metadata.column is not provided", {
  expect_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp))
})

test_that("errors are returned if the control is not provided", {
  expect_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                      metadata.column = "combined.id"))
})


test_that("errors are returned if the metadata column is not present", {
  expect_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                      metadata.column = "wrong_column"))
})


##########################################################################################

test_that("SAINTq function works on imputed counts", {
  expect_no_error(suppressMessages(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                          control = "BCa_6h.DMSO",
                                          metadata.column = "combined.id")))
})


test_that("SAINTq function works on raw counts", {
  expect_no_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                         control = "BCa_6h.DMSO",
                         metadata.column = "combined.id",
                         which.data = "raw",
                         verbose = FALSE))
})


test_that("SAINTq function works on normalized counts", {
  expect_no_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                         control = "BCa_6h.DMSO",
                         metadata.column = "combined.id",
                         which.data = "normalized",
                         verbose = FALSE))
})


test_that("SAINTq function works on randomized counts", {
  expect_no_error(SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
                         control = "BCa_6h.DMSO",
                         metadata.column = "combined.id",
                         which.data = "randomized",
                         verbose = FALSE))
})


###################################### methods ###########################################

test_that("SAINTq show-method works", {
  res = SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
               control = "BCa_6h.DMSO",
               metadata.column = "combined.id",
               which.data = "randomized",
               verbose = FALSE)

  expect_no_error(show(res))
})


test_that("SAINTq summary-method works", {
  res = SAINTq(DEprot.object = DEprot::test.toolbox$dpo.imp,
               control = "BCa_6h.DMSO",
               metadata.column = "combined.id",
               which.data = "randomized",
               verbose = FALSE)

  expect_no_error(summary(res))
})



