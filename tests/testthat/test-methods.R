test_that("method-show for class DEprot works", {
  expect_no_error(show(DEprot::test.toolbox$dpo.imp))
})


test_that("method-show for class DEprot.analyses works", {
  expect_no_error(show(DEprot::test.toolbox$diff.exp.limma))
})


test_that("method-summary for class DEprot.analyses works", {
  expect_no_error(summary(DEprot::test.toolbox$diff.exp.limma))
})


test_that("method-show for class DEprot.PCA works", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_no_error(show(pca))
})


# test_that("method-show for class DEprot.correlation works", {
#   correlation.heatmap = plot.correlation.heatmap(DEprot.object = DEprot::test.toolbox$dpo.imp)
#   expect_no_error(show(correlation.heatmap))
# })


# test_that("method-show for class DEprot.upset works", {
#   upset = plot.upset(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma)
#   expect_no_error(show(upset))
# })



# test_that("method-show for class DEprot.contrast.heatmap works", {
#   contrast.heatmap = heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5, use.uncorrected.pvalue = TRUE)
#   expect_no_error(show(contrast.heatmap))
# })


# test_that("method-show for class DEprot.counts.heatmap works", {
#   counts.heatmap = heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1)
#   expect_no_error(show(counts.heatmap))
# })


test_that("method-show for class DEprot.enrichResult works", {
  expect_no_error(show(DEprot::test.toolbox$gsea.results))
})



test_that("method-show for class DEprot.pvalues works", {
  expect_no_error(show(check.pvalues(DEprot::test.toolbox$diff.exp.limma, contrast = 1)))
})


test_that("method-show for class DEprot.normality works", {
  expect_no_error(show(suppressMessages(check.normality(DEprot.object = DEprot::test.toolbox$dpo.imp))))
})


test_that("method-show for class DEprot.RMSE works", {
  expect_no_error(show(suppressWarnings(compare.imp.methods(DEprot.object = DEprot::test.toolbox$dpo.norm, percentage.test = 100, sample.group.column = "combined.id", missForest.cores = 1, which.data = "normalized"))))
})

test_that("method-summary for class DEprot.RMSE works", {
  expect_no_error(summary(suppressWarnings(compare.imp.methods(DEprot.object = DEprot::test.toolbox$dpo.norm, percentage.test = 100, sample.group.column = "combined.id", missForest.cores = 1, which.data = "normalized"))))
})












