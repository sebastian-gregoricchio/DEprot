test_that("contrast.scatter works", {
  expect_no_error(contrast.scatter(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast.x = 1, contrast.y = 2))
})


test_that("expression.boxplot works", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14"))
})


test_that("expression.boxplot works (Z-score)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14", scale.expression = TRUE))
})


test_that("heatmap.contrasts works", {
  expect_no_error(heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5))
})


test_that("heatmap.contrasts works using the uncorrected p-value", {
  expect_no_error(heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5, use.uncorrected.pvalue = TRUE))
})


test_that("plot.MA works ", {
  expect_no_error(plot.MA(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})


test_that("plot.volcano works ", {
  expect_no_error(plot.volcano(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})



test_that("plot.PC.cumulative works ", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_no_error(plot.PC.cumulative(DEprot.PCA.object = pca))
})


test_that("plot.PC.scatter works ", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_no_error(plot.PC.scatter(DEprot.PCA.object = pca, color.column = "condition", shape.column = "replicate", label.column = "replicate"))
})

test_that("plot.PC.scatter.123 works ", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_no_error(plot.PC.scatter.123(DEprot.PCA.object = pca, color.column = "condition", shape.column = "replicate", label.column = "replicate"))
})



test_that("plot.counts works", {
  expect_no_error(plot.counts(DEprot.object = DEprot::test.toolbox$dpo.imp))
})


test_that("plot.upset works", {
  expect_no_error(plot.upset(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma))
})



test_that("protein.summary works", {
  expect_no_error(protein.summary(DEprot.object = DEprot::test.toolbox$dpo.imp, group.column = "combined.id"))
})

