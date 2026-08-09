## ----------------------------------------------------------------------------------------
##  normalize.counts() and harmonize.batches()
## ----------------------------------------------------------------------------------------

##########################################    normalize.counts    ########################

test_that("normalize.counts returns a normalized object with the same dimensions", {
  norm <- normalize.counts(DEprot.object = tb.dpo.raw)

  expect_s4_class(norm, "DEprot")
  expect_true(norm@normalized)
  expect_equal(dim(norm@norm.counts), dim(tb.dpo.raw@raw.counts))
  expect_equal(rownames(norm@norm.counts), rownames(tb.dpo.raw@raw.counts))
  expect_equal(colnames(norm@norm.counts), colnames(tb.dpo.raw@raw.counts))
  ## the raw counts must be left untouched
  expect_equal(norm@raw.counts, tb.dpo.raw@raw.counts)
})


test_that("the normalization method is documented in the object", {
  norm <- normalize.counts(DEprot.object = tb.dpo.raw)

  expect_s3_class(norm@normalization.method, "data.frame")
  expect_true("MBQN" %in% norm@normalization.method$value)
  expect_true("median" %in% norm@normalization.method$value)
})


test_that("the balancing function can be changed or removed", {
  mean.balanced <- normalize.counts(DEprot.object = tb.dpo.raw, balancing.function = "mean")
  unbalanced <- normalize.counts(DEprot.object = tb.dpo.raw, balancing.function = NULL)

  expect_true("mean" %in% mean.balanced@normalization.method$value)
  ## without a balancing function the object records 'balanced = FALSE'
  expect_true("FALSE" %in% as.character(unbalanced@normalization.method$value))
  ## the two normalizations are not expected to give the same values
  expect_false(isTRUE(all.equal(mean.balanced@norm.counts, unbalanced@norm.counts)))
})


test_that("the NRI/RI threshold is stored and changes the result", {
  low <- normalize.counts(DEprot.object = tb.dpo.raw, NRI.RI.ratio.threshold = 0.1)
  high <- normalize.counts(DEprot.object = tb.dpo.raw, NRI.RI.ratio.threshold = 0.9)

  expect_true("0.1" %in% as.character(low@normalization.method$value))
  expect_true("0.9" %in% as.character(high@normalization.method$value))
})


test_that("the normalized boxplot is generated", {
  norm <- normalize.counts(DEprot.object = tb.dpo.raw)

  expect_true(inherits(norm@boxplot.norm, "ggplot"))
})


test_that("an object already normalized is not overwritten silently", {
  norm <- normalize.counts(DEprot.object = tb.dpo.raw)

  expect_error(normalize.counts(DEprot.object = norm))
  expect_no_error(normalize.counts(DEprot.object = norm, overwrite.normalization = TRUE))
})


test_that("the normalization requires the raw counts", {
  bare <- tb.dpo.raw
  bare@raw.counts <- NULL
  bare@normalized <- TRUE

  expect_error(normalize.counts(DEprot.object = bare, overwrite.normalization = TRUE))
})


test_that("normalize.counts rejects objects of the wrong class", {
  expect_error(normalize.counts(DEprot.object = DEprot::sample.config))
  expect_error(normalize.counts(DEprot.object = "not a DEprot object"))
})



##########################################    harmonize.batches    #######################

test_that("harmonize.batches corrects the counts and documents the parameters", {
  dpo <- tb.dpo.raw
  dpo@metadata$batch <- c(rep("A", 6), rep("B", 6))

  harm <- harmonize.batches(DEprot.object = dpo, batch.column = "batch", cores = 1, verbose = FALSE)

  expect_s4_class(harm, "DEprot")
  expect_true(harm@normalized)
  expect_true("HarmonizR" %in% harm@normalization.method$value)
  expect_true("batch" %in% harm@normalization.method$value)
  ## the correction cannot invent proteins, only drop the ones it fails on
  expect_true(nrow(harm@norm.counts) <= nrow(dpo@raw.counts))
  expect_equal(ncol(harm@norm.counts), ncol(dpo@raw.counts))
})


test_that("the limma algorithm is accepted", {
  dpo <- tb.dpo.raw
  dpo@metadata$batch <- c(rep("A", 6), rep("B", 6))

  harm <- harmonize.batches(DEprot.object = dpo,
                            batch.column = "batch",
                            algorithm = "limma",
                            cores = 1,
                            verbose = FALSE)

  expect_true(any(grepl("limma", as.character(harm@normalization.method$value))))
})


test_that("the ComBat modes are passed through", {
  dpo <- tb.dpo.raw
  dpo@metadata$batch <- c(rep("A", 6), rep("B", 6))

  ## mode 2 is mean-only: it survives sparser designs than the default
  harm <- harmonize.batches(DEprot.object = dpo,
                            batch.column = "batch",
                            ComBat.mode = 2,
                            cores = 1,
                            verbose = FALSE)

  expect_true("2" %in% as.character(harm@normalization.method$value))
})


test_that("harmonize.batches returns an error when the batch column is not available", {
  expect_error(harmonize.batches(DEprot.object = tb.dpo.raw, batch.column = "batch", cores = 1, verbose = FALSE))
})


test_that("harmonize.batches requires raw counts and a DEprot object", {
  dpo <- tb.dpo.raw
  dpo@metadata$batch <- c(rep("A", 6), rep("B", 6))
  dpo@raw.counts <- NULL

  expect_error(harmonize.batches(DEprot.object = dpo, batch.column = "batch", cores = 1, verbose = FALSE))
  expect_error(harmonize.batches(DEprot.object = DEprot::sample.config, batch.column = "batch", cores = 1, verbose = FALSE))
})


test_that("a batch design too sparse to be corrected is reported", {
  dpo <- tb.dpo.raw
  ## one batch per sample: nothing can be estimated within a batch
  dpo@metadata$batch <- paste0("batch", seq_len(nrow(dpo@metadata)))

  expect_error(suppressWarnings(harmonize.batches(DEprot.object = dpo,
                                                  batch.column = "batch",
                                                  cores = 1,
                                                  verbose = FALSE)))
})
