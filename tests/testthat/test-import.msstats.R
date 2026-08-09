## ----------------------------------------------------------------------------------------
##  import.msstats()
##  MSstats is not a dependency of DEprot: the tests are run on tables having the same
##  structure as the ProteinLevelData returned by dataProcess()/proteinSummarization().
## ----------------------------------------------------------------------------------------

test_that("a label-free ProteinLevelData table is imported", {
  tb <- make.msstats.lfq()

  dpo <- suppressWarnings(import.msstats(object = tb, verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(nrow(dpo@norm.counts), 5)
  expect_equal(ncol(dpo@norm.counts), 4)
  expect_equal(dpo@log.base, 2)
})


test_that("the metadata is reconstructed from the MSstats annotation", {
  tb <- make.msstats.lfq()

  dpo <- suppressWarnings(import.msstats(object = tb, verbose = FALSE))

  expect_true("column.id" %in% colnames(dpo@metadata))
  ## GROUP and SUBJECT are constant within a run and must be carried over
  expect_true(any(c("GROUP", "SUBJECT") %in% colnames(dpo@metadata)))
  expect_equal(nrow(dpo@metadata), ncol(dpo@norm.counts))
})


test_that("the list returned by dataProcess() is accepted", {
  tb <- make.msstats.lfq()

  dpo <- suppressWarnings(import.msstats(object = list(ProteinLevelData = tb,
                                                       FeatureLevelData = data.frame()),
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
})


test_that("the legacy 'RunlevelData' element is accepted as well", {
  tb <- make.msstats.lfq()

  dpo <- suppressWarnings(import.msstats(object = list(RunlevelData = tb), verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
})


test_that("the MSstatsPTM nesting is resolved", {
  tb <- make.msstats.lfq()

  dpo <- suppressWarnings(import.msstats(object = list(PROTEIN = list(ProteinLevelData = tb)),
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
})


test_that("an .rds file is read from disk", {
  dir <- local.tmpdir()
  path <- file.path(dir, "msstats.rds")
  saveRDS(list(ProteinLevelData = make.msstats.lfq()), path)

  dpo <- suppressWarnings(import.msstats(object = path, verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
})


test_that("the data are flagged as imputed when MBimpute was used", {
  tb <- make.msstats.lfq()
  tb$NumImputedFeature <- c(1, rep(0, nrow(tb) - 1))

  dpo <- suppressWarnings(import.msstats(object = tb, verbose = FALSE))

  expect_true(dpo@imputed)
  expect_false(is.null(dpo@imputed.counts))
})


test_that("a TMT table is routed to the isobaric reader", {
  tb <- make.msstats.lfq()
  tb$Channel <- rep(c("channel.1", "channel.2"), length.out = nrow(tb))
  tb$Mixture <- "mix1"
  tb$Run <- tb$originalRUN
  tb$Abundance <- tb$LogIntensities
  tb$LogIntensities <- NULL

  dpo <- suppressWarnings(import.msstats(object = tb, verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  ## one column per run x channel combination
  expect_true(ncol(dpo@norm.counts) >= 2)
})



## ----------------------------------------------------------------------------------------
##  Error branches
## ----------------------------------------------------------------------------------------

test_that("a missing 'object' is reported", {
  expect_error(import.msstats(verbose = FALSE))
})


test_that("an object of the wrong class is rejected", {
  expect_error(import.msstats(object = 1:10, verbose = FALSE))
  expect_error(import.msstats(object = list(something = data.frame(a = 1)), verbose = FALSE))
})


test_that("a file that does not exist is reported", {
  expect_error(import.msstats(object = file.path(tempdir(), "nope.rds"), verbose = FALSE))
})


test_that("a table without the expected columns is rejected", {
  expect_error(import.msstats(object = data.frame(a = 1, b = 2), verbose = FALSE))
})


test_that("an unknown 'type' is rejected", {
  expect_error(import.msstats(object = make.msstats.lfq(), type = "not.a.type", verbose = FALSE))
})
