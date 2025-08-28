test_that("class DEprot is working", {
  expect_s4_class(DEprot::test.toolbox$dpo.raw, "DEprot")
})


test_that("class DEprot is working", {
  expect_s4_class(DEprot::test.toolbox$diff.exp.limma, "DEprot.analyses")
})


test_that("class DEprot.PCA is working", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_s4_class(pca, "DEprot.PCA")
})


test_that("class DEprot.correlation is working", {
  correlation = plot.correlation.heatmap(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_s4_class(correlation, "DEprot.correlation")
})


test_that("class DEprot.upset is working", {
  upset = plot.upset(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma)
  expect_s4_class(upset, "DEprot.upset")
})


test_that("class DEprot.contrast.heatmap is working", {
  contrast.heatmap = heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 2)
  expect_s4_class(contrast.heatmap, "DEprot.contrast.heatmap")
})



test_that("class DEprot.counts.heatmap is working", {
  counts.heatmap = heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma)
  expect_s4_class(counts.heatmap, "DEprot.counts.heatmap")
})



test_that("class DEprot.enrichResult is working", {
  expect_s4_class(DEprot::test.toolbox$gsea.results, "DEprot.enrichResult")
})



test_that("class DEprot.normality is working", {
  expect_s4_class(check.normality(DEprot.object = DEprot::test.toolbox$dpo.imp), "DEprot.normality")
})



test_that("class DEprot.pvalues is working", {
  expect_s4_class(check.pvalues(DEprot::test.toolbox$diff.exp.limma, contrast = 1), "DEprot.pvalues")
})


test_that("class DEprot.RMSE is working", {
  RMSE = suppressWarnings(compare.imp.methods(DEprot.object = DEprot::test.toolbox$dpo.norm, percentage.test = 100, sample.group.column = "combined.id", missForest.cores = 1,
                                              run.missForest = FALSE, run.BPCA = FALSE, run.kNN = FALSE, run.tkNN = FALSE, run.corkNN = FALSE, run.LLS = FALSE, run.SVD = TRUE, run.PPCA = FALSE, run.RegImpute = FALSE))
  expect_s4_class(RMSE, "DEprot.RMSE")
})







