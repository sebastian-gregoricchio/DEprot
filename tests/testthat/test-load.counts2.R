######################################    ERRORS    ######################################
test_that("errors are returned if the table is not numeric", {
  cnt = DEprot::unimputed.counts
  cnt[8,1:10] = "A"
  expect_error(load.counts2(counts = iris,
                            metadata = DEprot::sample.config[-5,],
                            data.type = "norm",
                            log.base = 2))

})



test_that("errors are returned if the counts and metadata do not match", {
  expect_error(load.counts2(counts = DEprot::unimputed.counts,
                            metadata = DEprot::sample.config[1:6,],
                            data.type = "norm",
                            log.base = 2))

})



test_that("errors are returned if the column.id idnidcated is not present in the metadata", {
  expect_error(load.counts2(counts = DEprot::unimputed.counts,
                            metadata = DEprot::sample.config[1:6,],
                            data.type = "norm",
                            log.base = 2,
                            column.id = "not a column"))

})


##########################################################################################

test_that("the function load.counts.2 is working with normalized data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "norm",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with raw data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "raw",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with imputed data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "imputed",
                               log.base = 2))

})


test_that("the function load.counts.2 is working with randomized data", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "randomized",
                               log.base = 2))

})



test_that("the function load.counts.2 is working with linear data", {
  expect_no_error(suppressMessages(load.counts2(counts = DEprot::unimputed.counts,
                                                metadata = DEprot::sample.config,
                                                data.type = "randomized",
                                                log.base = 1)))

})
