## ----------------------------------------------------------------------------------------
##  detect.outliers()
##  The function combines three independent metrics (correlation, Mahalanobis distance on
##  the PCs, missingness) and flags a sample when at least 'min.flags' of them fire.
## ----------------------------------------------------------------------------------------

test_that("the object returned carries the metrics of every sample", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         which.data = "imputed",
                         verbose = FALSE)

  expect_s4_class(out, "DEprot.outliers")
  expect_equal(nrow(out@metrics), nrow(tb.dpo.imp@metadata))
  expect_true("column.id" %in% colnames(out@metrics))
  expect_equal(out@data.used, "imputed")
  expect_equal(out@correlation.method, "pearson")
})


test_that("the correlation matrix is square and its diagonal is masked", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp, verbose = FALSE)

  expect_equal(dim(out@correlation.matrix), rep(length(out@sample.subset), 2))
  ## the self-correlation is set to NA: a sample must not be its own best partner
  expect_true(all(is.na(diag(out@correlation.matrix))))
  off.diagonal <- out@correlation.matrix[lower.tri(out@correlation.matrix)]
  expect_true(all(off.diagonal >= -1 & off.diagonal <= 1, na.rm = TRUE))
})


test_that("the spearman correlation is accepted", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         correlation.method = "spearman",
                         verbose = FALSE)

  expect_equal(out@correlation.method, "spearman")
})


test_that("a group column restricts the correlation to the samples of the same group", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         group.column = "condition",
                         verbose = FALSE)

  expect_equal(out@group.column, "condition")
  expect_s4_class(out, "DEprot.outliers")
})


test_that("the analysis can be restricted to a subset of samples", {
  samples <- tb.dpo.imp@metadata$column.id[1:6]

  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         sample.subset = samples,
                         verbose = FALSE)

  expect_equal(sort(out@sample.subset), sort(samples))
  expect_equal(nrow(out@metrics), length(samples))
})


test_that("the missingness metric is computed on unimputed counts", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         missingness.data = "auto",
                         verbose = FALSE)

  expect_true(!is.null(out@missingness.data.used))
  expect_true(length(out@metrics.available) > 0)
})


test_that("the thresholds change the number of flagged samples", {
  strict <- detect.outliers(DEprot.object = tb.dpo.imp,
                            correlation.z.th = 0,
                            mahalanobis.padj.th = 1,
                            missingness.z.th = 0,
                            min.flags = 1,
                            verbose = FALSE)

  loose <- detect.outliers(DEprot.object = tb.dpo.imp,
                           correlation.z.th = -100,
                           mahalanobis.padj.th = 0,
                           missingness.z.th = 100,
                           min.flags = 3,
                           verbose = FALSE)

  expect_true(length(strict@outliers) >= length(loose@outliers))
  expect_equal(length(loose@outliers), 0)
})


test_that("the absolute thresholds are accepted alongside the z-scores", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         correlation.min = 0.99,
                         missingness.max = 0.01,
                         min.flags = 1,
                         verbose = FALSE)

  expect_s4_class(out, "DEprot.outliers")
  expect_equal(out@parameters$correlation.min, 0.99)
  expect_equal(out@parameters$missingness.max, 0.01)
})


test_that("the number of PCs and the scaling can be tuned", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         n.PCs = 2,
                         center.data = FALSE,
                         scale.data = FALSE,
                         verbose = FALSE)

  expect_s4_class(out, "DEprot.outliers")
  expect_true(!is.null(out@PCA))
})


test_that("the p-value correction method can be changed", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp,
                         padj.method = "bonferroni",
                         verbose = FALSE)

  expect_equal(out@parameters$padj.method, "bonferroni")
})


test_that("the diagnostic plots are stored in the object", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp, verbose = FALSE)

  expect_true(length(out@plot.list) > 0)
  expect_true(inherits(out@plot, "ggplot") | inherits(out@plot, "patchwork"))
})


test_that("the show method does not fail", {
  out <- detect.outliers(DEprot.object = tb.dpo.imp, verbose = FALSE)

  expect_no_error(suppressWarnings(suppressMessages(show(out))))
})



## ----------------------------------------------------------------------------------------
##  Error branches
## ----------------------------------------------------------------------------------------

test_that("an object of the wrong class is rejected", {
  expect_error(detect.outliers(DEprot.object = data.frame(a = 1), verbose = FALSE))
})


test_that("counts that are not available cannot be used", {
  dpo.bare <- tb.dpo.norm
  dpo.bare@imputed.counts <- NULL

  expect_error(detect.outliers(DEprot.object = dpo.bare, which.data = "imputed", verbose = FALSE))
  expect_error(detect.outliers(DEprot.object = dpo.bare, which.data = "not.a.type", verbose = FALSE))
})


test_that("a group column that does not exist is rejected", {
  expect_error(detect.outliers(DEprot.object = tb.dpo.imp,
                               group.column = "not.a.column",
                               verbose = FALSE))
})


test_that("a sample subset that does not match the counts is rejected", {
  expect_error(detect.outliers(DEprot.object = tb.dpo.imp,
                               sample.subset = c("not.a.sample"),
                               verbose = FALSE))
})


test_that("too few samples cannot support the multivariate metrics", {
  expect_error(detect.outliers(DEprot.object = tb.dpo.imp,
                               sample.subset = tb.dpo.imp@metadata$column.id[1:2],
                               verbose = FALSE))
})
