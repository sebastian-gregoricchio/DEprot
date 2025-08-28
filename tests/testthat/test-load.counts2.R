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

test_that("the function load.counts.2 is working", {
  expect_no_error(load.counts2(counts = DEprot::unimputed.counts,
                               metadata = DEprot::sample.config,
                               data.type = "norm",
                               log.base = 2))

})
