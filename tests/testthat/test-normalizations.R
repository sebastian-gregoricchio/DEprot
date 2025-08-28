test_that("harmonize.batches works", {
  dpo = DEprot::test.toolbox$dpo.raw
  dpo@metadata$batch = c(rep("A",6), rep("B",6))
  expect_no_error(harmonize.batches(DEprot.object = dpo, batch.column = "batch", cores = 1))
})

test_that("harmonize.batches returns an error when the batch column is not available", {
  expect_error(harmonize.batches(DEprot.object = DEprot::test.toolbox$dpo.raw, batch.column = "batch", cores = 1))
})


test_that("normalize.counts works", {
  expect_no_error(normalize.counts(DEprot.object = DEprot::test.toolbox$dpo.raw))
})


