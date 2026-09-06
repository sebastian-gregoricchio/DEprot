## ----------------------------------------------------------------------------------------
##  harmonize.batches() and detect.outliers(): the remaining branches
## ----------------------------------------------------------------------------------------

#' Object carrying a usable batch design
#'
#' The batches must be blocks of samples: an alternating design leaves too few proteins
#' quantified in both batches for ComBat to estimate anything.

batched.object <-
  function(object = tb.dpo.raw,
           n.batches = 2) {
    object@metadata$batch <- rep(paste0("batch", seq_len(n.batches)),
                                 each = ceiling(nrow(object@metadata) / n.batches))[seq_len(nrow(object@metadata))]
    return(object)
  }


harmonize <-
  function(object = batched.object(), ...) {
    suppressWarnings(suppressMessages(
      harmonize.batches(DEprot.object = object, batch.column = "batch", cores = 1, verbose = FALSE, ...)))
  }



## ----------------------------------------------------------------------------------------
##  harmonize.batches
## ----------------------------------------------------------------------------------------

test_that("the corrected counts land in the normalized slot and keep the sample order", {
  object <- batched.object()
  harm <- harmonize(object)

  expect_s4_class(harm, "DEprot")
  expect_true(harm@normalized)
  expect_equal(colnames(harm@norm.counts), colnames(object@raw.counts))
  ## the raw counts are the input and must be left untouched
  expect_equal(harm@raw.counts, object@raw.counts)
  expect_true(is.matrix(harm@norm.counts))
})


test_that("the parameters of the correction are stored in the object", {
  harm <- harmonize()
  parameters <- harm@normalization.method

  expect_s3_class(parameters, "data.frame")
  expect_setequal(parameters$param,
                  c("package", "batch.column", "algorithm", "ComBat.mode", "block", "cores"))
  expect_equal(parameters$value[parameters$param == "package"], "HarmonizR")
  expect_equal(parameters$value[parameters$param == "batch.column"], "batch")
  expect_equal(parameters$value[parameters$param == "block"], "null")
})


test_that("every ComBat mode is accepted", {
  for (mode in 1:4) {
    harm <- harmonize(ComBat.mode = mode)
    expect_s4_class(harm, "DEprot")
    expect_equal(harm@normalization.method$value[harm@normalization.method$param == "ComBat.mode"],
                 as.character(mode))
  }
})


test_that("the limma algorithm is accepted and reported", {
  harm <- harmonize(algorithm = "limma")

  expect_s4_class(harm, "DEprot")
  expect_equal(harm@normalization.method$value[harm@normalization.method$param == "algorithm"], "limma")
  ## the boxplot of the normalized counts is rebuilt
  expect_true(inherits(harm@boxplot.norm, "ggplot"))
})


test_that("a blocking variable can be passed to HarmonizR", {
  harm <- harmonize(block = 2)

  expect_s4_class(harm, "DEprot")
  expect_equal(harm@normalization.method$value[harm@normalization.method$param == "block"], "2")
})


test_that("three batches are handled as well as two", {
  harm <- harmonize(object = batched.object(n.batches = 3))

  expect_s4_class(harm, "DEprot")
  expect_equal(ncol(harm@norm.counts), ncol(tb.dpo.raw@raw.counts))
})


test_that("the proteins that could not be corrected are reported", {
  ## HarmonizR drops the proteins it cannot model: the warning names how many were lost
  object <- batched.object()

  harm <- harmonize(object)

  expect_true(nrow(harm@norm.counts) <= nrow(object@raw.counts))
  expect_true(all(rownames(harm@norm.counts) %in% rownames(object@raw.counts)))
})


test_that("a design too sparse to be corrected is reported", {
  ## one batch per sample: nothing can be estimated within a batch
  object <- batched.object(n.batches = nrow(tb.dpo.raw@metadata))

  expect_error(harmonize(object))
})


test_that("the wrong inputs are rejected", {
  expect_error(harmonize.batches(DEprot.object = DEprot::sample.config,
                                 batch.column = "batch", cores = 1, verbose = FALSE))

  expect_error(harmonize.batches(DEprot.object = batched.object(),
                                 batch.column = "not.a.column", cores = 1, verbose = FALSE))

  no.raw <- batched.object()
  no.raw@raw.counts <- NULL
  expect_error(harmonize.batches(DEprot.object = no.raw,
                                 batch.column = "batch", cores = 1, verbose = FALSE))
})


test_that("the messages of HarmonizR can be printed", {
  ## 'verbose' is forwarded as the verbosity level of HarmonizR
  expect_no_error(suppressWarnings(suppressMessages(
    harmonize.batches(DEprot.object = batched.object(),
                      batch.column = "batch",
                      cores = 1,
                      verbose = TRUE))))
})



## ----------------------------------------------------------------------------------------
##  detect.outliers: the metric-availability branches
## ----------------------------------------------------------------------------------------

test_that("every 'missingness.data' keyword is accepted", {
  for (w in c("auto", "raw", "normalized", "randomized")) {
    out <- detect.outliers(DEprot.object = tb.dpo.imp,
                           missingness.data = w,
                           verbose = FALSE)

    expect_s4_class(out, "DEprot.outliers")
    expect_true(is.character(out@missingness.data.used))
  }

  expect_error(detect.outliers(DEprot.object = tb.dpo.imp,
                               missingness.data = "not.a.type",
                               verbose = FALSE))
})


test_that("the missingness metric is dropped when the counts hold no missing value", {
  ## the imputed counts are complete: every sample has a missing rate of zero, hence the
  ## metric carries no information and must be declared unavailable
  complete <- tb.dpo.imp
  complete@raw.counts <- complete@imputed.counts
  complete@norm.counts <- complete@imputed.counts
  complete@random.counts <- NULL

  out <- suppressWarnings(detect.outliers(DEprot.object = complete,
                                          missingness.data = "raw",
                                          verbose = FALSE))

  expect_s4_class(out, "DEprot.outliers")
  expect_false("missingness" %in% out@metrics.available)
})


test_that("an absolute threshold keeps a metric that the z-scores would have dropped", {
  complete <- tb.dpo.imp
  complete@raw.counts <- complete@imputed.counts
  complete@norm.counts <- complete@imputed.counts

  out <- suppressWarnings(detect.outliers(DEprot.object = complete,
                                          missingness.data = "raw",
                                          missingness.max = 0.5,
                                          min.flags = 1,
                                          verbose = FALSE))

  expect_s4_class(out, "DEprot.outliers")
})


test_that("'min.flags' is capped by the number of metrics available", {
  out <- suppressWarnings(detect.outliers(DEprot.object = tb.dpo.imp,
                                          min.flags = 10,
                                          verbose = FALSE))

  expect_true(out@parameters$min.flags <= length(out@metrics.available))
})


test_that("the number of PCs is capped by the samples available", {
  out <- suppressWarnings(detect.outliers(DEprot.object = tb.dpo.imp,
                                          n.PCs = 50,
                                          verbose = FALSE))

  expect_s4_class(out, "DEprot.outliers")
  expect_true(!is.null(out@PCA))
})


test_that("the messages are printed when requested", {
  expect_message(detect.outliers(DEprot.object = tb.dpo.imp, verbose = TRUE))
})


test_that("a sample flagged by every metric is reported as an outlier", {
  ## one sample is made deliberately aberrant: shifted in intensity and mostly missing
  object <- tb.dpo.imp
  target <- object@metadata$column.id[1]

  object@imputed.counts[, target] <- rev(object@imputed.counts[, target]) + 10
  object@norm.counts[, target] <- NA
  object@norm.counts[1:2, target] <- object@imputed.counts[1:2, target]

  out <- suppressWarnings(detect.outliers(DEprot.object = object,
                                          missingness.data = "normalized",
                                          min.flags = 1,
                                          verbose = FALSE))

  expect_true(target %in% out@outliers)
  expect_true(out@metrics$outlier[out@metrics$column.id == target])
  ## the summary printed names the outliers instead of 'none'
  expect_no_error(suppressMessages(show(out)))
})


test_that("the group column restricts the correlation to the samples of the same group", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         group.column = "condition",
                         verbose = FALSE)

  expect_equal(out@group.column, "condition")
  expect_true(all(out@metrics$group %in% as.character(tb.dpo.imp@metadata$condition)))
})


test_that("a group holding a single sample falls back on the other samples", {
  ## with one sample per group there is no within-group partner: the correlation is then
  ## computed against all the other samples rather than being left undefined
  object <- tb.dpo.imp
  object@metadata$lonely <- object@metadata$column.id

  out <- suppressWarnings(detect.outliers(DEprot.object = object,
                                          group.column = "lonely",
                                          verbose = FALSE))

  expect_s4_class(out, "DEprot.outliers")
  expect_false(all(is.na(out@metrics$median.correlation)))
})
