## ----------------------------------------------------------------------------------------
##  compare.ranking()
## ----------------------------------------------------------------------------------------

test_that("the returned object is the panel of six correlation plots", {
  cmp <- compare.ranking(DEprot.analyses.object = tb.limma, contrast = 1)

  expect_true(inherits(cmp, "patchwork") | inherits(cmp, "ggplot"))
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(cmp[[1]])))
})


test_that("every contrast can be compared", {
  for (i in seq_along(tb.limma@contrasts)) {
    expect_no_error(compare.ranking(DEprot.analyses.object = tb.limma, contrast = i))
  }
})


test_that("the contrast must be given as a number", {
  expect_error(compare.ranking(DEprot.analyses.object = tb.limma,
                               contrast = names(tb.limma@contrasts)[1]))
})


test_that("the colors can be customized", {
  expect_no_error(compare.ranking(DEprot.analyses.object = tb.limma,
                                  contrast = 1,
                                  color.up = "darkred",
                                  color.down = "navy",
                                  regression.line.color = "black"))
})


test_that("errors are returned for wrong inputs", {
  expect_error(compare.ranking(DEprot.analyses.object = tb.dpo.imp, contrast = 1))
  expect_error(compare.ranking(DEprot.analyses.object = DEprot::sample.config, contrast = 1))
  expect_error(compare.ranking(DEprot.analyses.object = tb.limma, contrast = 100))
  expect_error(compare.ranking(DEprot.analyses.object = tb.limma))
})




## ----------------------------------------------------------------------------------------
##  compare.imp.methods()
##  The full comparison runs every imputation method available: the tests below switch the
##  slow ones off and keep one or two representative methods, so that the function is
##  walked through without paying for the whole benchmark.
## ----------------------------------------------------------------------------------------

## the arguments shared by all the calls below: only the pcaMethods-based methods are run
imp.args <- list(percentage.test = 100,
                 sample.group.column = "combined.id",
                 which.data = "normalized",
                 run.missForest = FALSE,
                 run.kNN = FALSE,
                 run.tkNN = FALSE,
                 run.corkNN = FALSE,
                 run.LLS = FALSE,
                 run.SVD = TRUE,
                 run.BPCA = FALSE,
                 run.PPCA = FALSE,
                 run.RegImpute = FALSE,
                 pcaMethods.nPCs.to.test = 2,
                 seed = 1,
                 verbose = FALSE)

run.comparison <-
  function(...) {
    suppressWarnings(do.call(compare.imp.methods,
                             c(list(DEprot.object = tb.dpo.norm), utils::modifyList(imp.args, list(...)))))
  }



test_that("errors are returned if the object is not of class DEprot", {
  expect_error(compare.imp.methods(DEprot.object = "not DEprot object",
                                   percentage.test = 100,
                                   sample.group.column = "combined.id"))
})


test_that("a group column that does not exist is rejected", {
  expect_error(compare.imp.methods(DEprot.object = tb.dpo.norm,
                                   percentage.test = 100,
                                   sample.group.column = "not.a.column",
                                   which.data = "normalized"))
})


test_that("counts that are not available cannot be compared", {
  bare <- tb.dpo.raw
  bare@norm.counts <- NULL
  bare@random.counts <- NULL

  expect_error(compare.imp.methods(DEprot.object = bare,
                                   percentage.test = 100,
                                   sample.group.column = "combined.id",
                                   which.data = "randomized"))
})



##########################################################################################

test_that("the RMSE object carries the scores of the methods that were run", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison(run.LLS = TRUE)

  expect_s4_class(cmp, "DEprot.RMSE")
  expect_equal(cmp@percentage.test, 100)
  expect_equal(cmp@seed, 1)
  expect_true(length(cmp@imputed.objects) >= 1)
  expect_true(length(cmp@RMSE.tables) >= 1)
  expect_true(nrow(cmp@RMSE.scores) >= 1)
})


test_that("the test dataset and the missingness are stored", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison()

  expect_s4_class(cmp@original.DEprot.object, "DEprot")
  ## the test dataset is the subset of proteins that were fully quantified, on which the
  ## missing values are then simulated: it cannot have more rows than the original counts
  expect_true(nrow(cmp@test.dataset) <= nrow(tb.dpo.norm@norm.counts))
  expect_equal(ncol(cmp@test.dataset), ncol(tb.dpo.norm@norm.counts))

  ## with a 'sample.group.column' the missingness is reported group by group
  expect_true(all(as.numeric(unlist(cmp@fraction.missing.values[sapply(cmp@fraction.missing.values, is.numeric)])) >= 0))
})


test_that("the missingness is a single value when no group column is provided", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison(sample.group.column = NULL)

  expect_true(all(as.numeric(unlist(cmp@fraction.missing.values)) >= 0))
  expect_s4_class(cmp, "DEprot.RMSE")
})


test_that("the diagnostic plots are generated", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison(normalize.color.bar = FALSE,
                        low.residual.color = "blue",
                        zero.residual.color = "grey90",
                        high.residual.color = "red")

  expect_true(length(cmp@correlation.plots) >= 1)
  expect_true(!is.null(cmp@density.residuals))
})


test_that("the spearman correlation is accepted", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison(correlation.method = "spearman")

  expect_s4_class(cmp, "DEprot.RMSE")
})


test_that("the show and summary methods return the report", {
  skip_if_not_installed("pcaMethods")

  cmp <- run.comparison()

  ## show() assembles the correlation plots, summary() returns the RMSE scores
  expect_no_error(suppressWarnings(show(cmp)))
  expect_s3_class(summary(cmp), "data.frame")
})


test_that("the full benchmark runs with the default methods", {
  skip_on_cran()
  skip_if_not_installed("missForest")
  skip_if_not_installed("pcaMethods")

  expect_no_error(suppressWarnings(compare.imp.methods(DEprot.object = tb.dpo.norm,
                                                       percentage.test = 100,
                                                       sample.group.column = "combined.id",
                                                       missForest.cores = 1,
                                                       missForest.max.iterations = 3,
                                                       which.data = "normalized",
                                                       seed = 1,
                                                       verbose = FALSE)))
})
