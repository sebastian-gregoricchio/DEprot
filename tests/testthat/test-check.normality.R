## ----------------------------------------------------------------------------------------
##  check.normality()
## ----------------------------------------------------------------------------------------

######################################    ERRORS    ######################################

test_that("errors are returned if the object is not of class DEprot", {
  expect_error(check.normality(DEprot.object = DEprot::sample.config))
  expect_error(check.normality(DEprot.object = "not an object"))
})


test_that("errors are returned if the data required are not correct", {
  expect_error(check.normality(DEprot.object = tb.dpo.imp, which.data = "wrong data"))
})


test_that("counts that are not available cannot be checked", {
  bare <- tb.dpo.raw
  bare@norm.counts <- NULL
  bare@imputed.counts <- NULL
  bare@random.counts <- NULL

  expect_error(check.normality(DEprot.object = bare, which.data = "normalized", verbose = FALSE))
  expect_error(check.normality(DEprot.object = bare, which.data = "imputed", verbose = FALSE))
  expect_error(check.normality(DEprot.object = bare, which.data = "randomized", verbose = FALSE))
})



##########################################################################################

test_that("the normality object carries one test and two plots per sample", {
  norm <- check.normality(DEprot.object = tb.dpo.imp, verbose = FALSE)

  expect_s4_class(norm, "DEprot.normality")
  expect_true(is.logical(norm@norm.statement))
  expect_equal(length(norm@qqplots), ncol(tb.dpo.imp@imputed.counts))
  expect_equal(length(norm@densities), ncol(tb.dpo.imp@imputed.counts))
  expect_true(inherits(norm@qqplots[[1]], "ggplot"))
  expect_true(inherits(norm@densities[[1]], "ggplot"))
})


test_that("the Anderson-Darling tests are stored one per sample", {
  norm <- check.normality(DEprot.object = tb.dpo.imp, verbose = FALSE)

  ## 'norm.AD.tests' is a list of htest objects returned by nortest::ad.test()
  expect_type(norm@norm.AD.tests, "list")
  expect_equal(length(norm@norm.AD.tests), ncol(tb.dpo.imp@imputed.counts))
  expect_s3_class(norm@norm.AD.tests[[1]], "htest")

  p.values <- vapply(norm@norm.AD.tests, FUN = function(x) {x$p.value}, FUN.VALUE = numeric(1))
  expect_true(all(p.values >= 0 & p.values <= 1))
})


test_that("the p-value threshold is stored and drives the verdict", {
  permissive <- check.normality(DEprot.object = tb.dpo.imp, p.threshold = 1e-10, verbose = FALSE)
  strict <- check.normality(DEprot.object = tb.dpo.imp, p.threshold = 0.99, verbose = FALSE)

  expect_equal(permissive@p.threshold, 1e-10)
  expect_equal(strict@p.threshold, 0.99)
  ## a threshold of ~1 rejects the normality of essentially any real distribution
  expect_false(strict@norm.statement)
})


test_that("every count type available can be checked", {
  for (w in c("raw", "normalized", "imputed")) {
    norm <- check.normality(DEprot.object = tb.dpo.imp, which.data = w, verbose = FALSE)
    expect_s4_class(norm, "DEprot.normality")
  }
})


test_that("the show method reports the verdict", {
  norm <- check.normality(DEprot.object = tb.dpo.imp, verbose = FALSE)

  ## the method writes through message(), hence on the standard error stream
  expect_message(show(norm))
})


test_that("the plot method assembles the qqplots and the densities", {
  norm <- check.normality(DEprot.object = tb.dpo.imp, verbose = FALSE)

  p <- suppressWarnings(plot(norm, n.samples = 2))
  expect_true(inherits(p, "patchwork") | inherits(p, "ggplot"))
})




## ----------------------------------------------------------------------------------------
##  check.pvalues()
## ----------------------------------------------------------------------------------------

######################################    ERRORS    ######################################

test_that("errors are returned if the object is not of class DEprot.analyses", {
  expect_error(check.pvalues(DEprot::sample.config))
  expect_error(check.pvalues(tb.dpo.imp))
})


test_that("errors are returned if the contrast indicated is not available", {
  expect_error(check.pvalues(tb.limma, contrast = 100))
})


test_that("the contrast must be given as a number", {
  expect_error(check.pvalues(tb.limma, contrast = names(tb.limma@contrasts)[1]))
})



##########################################################################################

test_that("the object carries the three diagnostic plots", {
  pv <- check.pvalues(tb.limma, contrast = 1)

  expect_s4_class(pv, "DEprot.pvalues")
  expect_true(inherits(pv@pvalue.distribution, "ggplot"))
  expect_true(inherits(pv@padjusted.distribution, "ggplot"))
  expect_true(inherits(pv@pvalue.rank, "ggplot"))
})


test_that("every contrast of the object can be checked", {
  for (i in seq_along(tb.limma@contrasts)) {
    expect_s4_class(check.pvalues(tb.limma, contrast = i), "DEprot.pvalues")
  }
})


test_that("the aesthetic arguments are accepted", {
  pv <- check.pvalues(tb.limma,
                      contrast = 1,
                      histogram.binwidth = 0.1,
                      p.value.color = "forestgreen",
                      p.adjusted.color = "orange")

  expect_s4_class(pv, "DEprot.pvalues")
  ## the plots must be buildable with the custom parameters
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(pv@pvalue.distribution)))
})


test_that("the show method assembles the plots", {
  pv <- check.pvalues(tb.limma, contrast = 1)

  expect_no_error(suppressWarnings(show(pv)))
})
