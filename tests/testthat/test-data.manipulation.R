test_that("filter.proteins works", {
  expect_no_error(filter.proteins(DEprot.object = DEprot::test.toolbox$diff.exp.limma, proteins = c("protein.29","protein.30"), mode = "remove"))
})



test_that("get.metadata works from DEprot", {
  expect_no_error(get.metadata(DEprot::test.toolbox$dpo.raw))
})

test_that("get.metadata works from diff.analyses", {
  expect_no_error(get.metadata(DEprot::test.toolbox$diff.exp.limma))
})

test_that("get.metadata works from PCA", {
  expect_no_error(get.metadata(DEprot::test.toolbox$diff.exp.limma@analyses.result.list$condition_FBS.vs.6h.DMSO$PCA.data))
})

test_that("get.metadata works from correlation", {
  expect_no_error(get.metadata(plot.correlation.heatmap(DEprot.object = DEprot::test.toolbox$dpo.imp)))
})

test_that("get.metadata does not work for not DEprot objects", {
  expect_error(get.metadata(DEprot.object = DEprot::test.toolbox$geneset))
})






test_that("get.results works", {
  expect_no_error(get.results(DEprot::test.toolbox$diff.exp.limma, 1))
})


test_that("randomize.missing.value works", {
  expect_no_error(suppressMessages(randomize.missing.values(DEprot.object = DEprot::test.toolbox$dpo.norm, group.column = "combined.id")))
})




test_that("remove.undetected.proteins works (on raw data)", {
  expect_no_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.imp, min.n.samples = 1, which.data = "raw"))
})

test_that("remove.undetected.proteins works (on normalized data)", {
  expect_no_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.imp, min.n.samples = 1, which.data = "normalized"))
})

test_that("remove.undetected.proteins works (on imputed data)", {
  expect_no_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.imp, min.n.samples = 1, which.data = "imputed"))
})

test_that("remove.undetected.proteins returns an error if raw data are required but not available", {
  dpo = DEprot::test.toolbox$dpo.imp
  dpo@raw.counts = NULL
  expect_error(remove.undetected.proteins(DEprot.object = dpo, min.n.samples = 1, which.data = "raw"))
})

test_that("remove.undetected.proteins returns an error if normalized data are required but not available", {
  expect_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.raw, min.n.samples = 1, which.data = "norm"))
})

test_that("remove.undetected.proteins returns an error if imputed data are required but not available", {
  expect_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.norm, min.n.samples = 1, which.data = "imp"))
})


test_that("remove.undetected.proteins returns an error if the type of data required is not recognized", {
  expect_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.norm, min.n.samples = 1, which.data = "hello"))
})





test_that("rename.samples works", {
  expect_no_error(rename.samples(DEprot.object = DEprot::test.toolbox$dpo.imp, metadata.column = "column.id"))
})


test_that("rescale.bait works", {
  expect_no_error(rescale.bait(DEprot.object = DEprot::test.toolbox$dpo.imp, bait.id = "protein.29"))
})
