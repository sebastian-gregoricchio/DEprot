######## contrast.scatter
test_that("contrast.scatter works", {
  expect_no_error(contrast.scatter(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast.x = 1, contrast.y = 2))
})



######### expression.boxplot
test_that("expression.boxplot works", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14", shape.column = "replicate"))
})

test_that("expression.boxplot works (grouped)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14", group.by.metadata.column = "combined.id"))
})

test_that("expression.boxplot works (Z-score)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp, protein.id = "protein.14", scale.expression = TRUE))
})

test_that("expression.boxplot works (raw data)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp, protein.id = "protein.14", scale.expression = TRUE, which.data = "raw"))
})

test_that("expression.boxplot works (randomized data)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp, protein.id = "protein.14", scale.expression = TRUE, which.data = "random"))
})

test_that("expression.boxplot works (normalized data)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp, protein.id = "protein.14", scale.expression = TRUE, which.data = "norm"))
})

test_that("expression.boxplot works (with pairwise)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14", shape.column = "replicate", pairwise.comparisons = TRUE))
})

test_that("expression.boxplot works (with pairwise, stars)", {
  expect_no_error(expression.boxplot(DEprot.object = DEprot::test.toolbox$dpo.imp,protein.id = "protein.14", shape.column = "replicate", pairwise.comparisons = TRUE, pairwise.p.label = "stars"))
})




####### Hetamaps
test_that("heatmap.contrasts works", {
  expect_no_error(heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5))
})


test_that("heatmap.contrasts works using the uncorrected p-value", {
  expect_no_error(heatmap.contrasts(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5, use.uncorrected.pvalue = TRUE))
})



test_that("heatmap.counts works for differential analyses", {
  expect_no_error(heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})

test_that("heatmap.counts works for differential analyses (scaled rows)", {
  expect_no_error(heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, scale = "row"))
})

test_that("heatmap.counts works for differential analyses (scaled columns)", {
  expect_no_error(heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, scale = "column"))
})

test_that("heatmap.counts works for differential analyses (top.n proteins)", {
  expect_no_error(heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5, contrast = 1, scale = "column"))
})


test_that("heatmap.counts works for differential analyses (average groups)", {
  expect_no_error(heatmap.counts(DEprot.object = DEprot::test.toolbox$diff.exp.limma, top.n = 5, contrast = 1, scale = "column", group.by.metadata.column = "combined.id"))
})




### MA plot
test_that("plot.MA works ", {
  expect_no_error(plot.MA(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})

test_that("plot.MA works with uncorrected pvalue (label in box)", {
  expect_no_error(plot.MA(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, use.uncorrected.pvalue = T, dot.labels = "protein.19", labels.in.boxes = TRUE))
})

test_that("plot.MA works with uncorrected pvalue (label as text)", {
  expect_no_error(plot.MA(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, use.uncorrected.pvalue = T, dot.labels = "protein.19", labels.in.boxes = FALSE))
})




test_that("plot.volcano works ", {
  expect_no_error(plot.volcano(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1))
})

test_that("plot.volcano works with uncorrected pvalue", {
  expect_no_error(plot.volcano(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, use.uncorrected.pvalue = TRUE))
})

test_that("plot.volcano works (label in box)", {
  expect_no_error(plot.volcano(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, dot.labels = "protein.19", labels.in.boxes = TRUE))
})


test_that("plot.volcano works (label as text)", {
  expect_no_error(plot.volcano(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma, contrast = 1, dot.labels = "protein.19", labels.in.boxes = FALSE, label.top.n = 1, use.uncorrected.pvalue = TRUE))
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


test_that("plot.PC.biplot works ", {
  pca = perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
  expect_no_error(plot.PC.biplot(DEprot.PCA.object = pca, color.column = "condition", shape.column = "replicate", label.column = "replicate"))
})




test_that("plot.counts works", {
  expect_no_error(plot.counts(DEprot.object = DEprot::test.toolbox$dpo.imp))
})


test_that("plot.upset works", {
  expect_no_error(plot.upset(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma))
})



test_that("protein.summary works (frequency)", {
  expect_no_error(protein.summary(DEprot.object = DEprot::test.toolbox$dpo.imp, group.column = "combined.id", n.labels = "frequency", show.frequency = TRUE))
})

test_that("protein.summary works (percentage)", {
  expect_no_error(protein.summary(DEprot.object = DEprot::test.toolbox$dpo.imp, group.column = "combined.id", n.labels = "percentage"))
})

test_that("protein.summary works (counts)", {
  expect_no_error(protein.summary(DEprot.object = DEprot::test.toolbox$dpo.imp, group.column = "combined.id", n.labels = "counts"))
})





test_that("contrast.LFQ works (text label)", {
  expect_no_error(contrast.LFQ(DEprot::test.toolbox$diff.exp.limma, dot.labels = "protein.17", label.font.size = 3, labels.in.boxes = TRUE, show.only.differential = TRUE, show.only.significant = TRUE, log2FC.scale.min = -1.1, log2FC.scale.max = 1.1))
})

test_that("contrast.LFQ works (box label)", {
  expect_no_error(contrast.LFQ(DEprot::test.toolbox$diff.exp.limma, dot.labels = "protein.17", label.font.size = 3, labels.in.boxes = FALSE, show.only.differential = TRUE, show.only.significant = TRUE, log2FC.scale.min = -1.1, log2FC.scale.max = 1.1))
})

test_that("contrast.LFQ does not work when wrong contrast provided", {
  expect_error(contrast.LFQ(DEprot::test.toolbox$diff.exp.limma, contrast = 100))
})

test_that("contrast.LFQ does not work when the object is not of class DEprot.analyses", {
  expect_error(contrast.LFQ(DEprot::test.toolbox$dpo.raw))
})





