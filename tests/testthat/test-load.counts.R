test_that("errors are returned if the table is not numeric", {
  cnt = DEprot::unimputed.counts
  cnt[8,1:10] = "A"
  expect_error(load.counts(counts = iris,
                           metadata = DEprot::sample.config[-5,],
                           log.base = 2))

})



test_that("errors are returned if the counts and metadata do not match", {
  expect_error(load.counts(counts = DEprot::unimputed.counts,
                           metadata = DEprot::sample.config[1:6,],
                           log.base = 2))

})



test_that("the function is working", {
  expect_no_error(load.counts(counts = DEprot::unimputed.counts,
                              metadata = DEprot::sample.config,
                              log.base = 2))

})
