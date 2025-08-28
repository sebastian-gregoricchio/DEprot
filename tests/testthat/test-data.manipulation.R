test_that("filter.proteins works", {
  expect_no_error(filter.proteins(DEprot.object = DEprot::test.toolbox$diff.exp.limma, proteins = c("protein.29","protein.30"), mode = "remove"))
})


test_that("get.metadata works", {
  expect_no_error(get.metadata(DEprot::test.toolbox$diff.exp.limma))
})


test_that("get.results works", {
  expect_no_error(get.results(DEprot::test.toolbox$diff.exp.limma, 1))
})


test_that("randomize.missing.value works", {
  expect_no_error(suppressMessages(randomize.missing.values(DEprot.object = DEprot::test.toolbox$dpo.norm, group.column = "combined.id")))
})

test_that("remove.undetected.proteins works", {
  expect_no_error(remove.undetected.proteins(DEprot.object = DEprot::test.toolbox$dpo.norm, min.n.samples = 1))
})

test_that("remove.undetected.proteins works", {
  expect_no_error(rename.samples(DEprot.object = DEprot::test.toolbox$dpo.imp, metadata.column = "column.id"))
})


test_that("rescale.bait works", {
  expect_no_error(rescale.bait(DEprot.object = DEprot::test.toolbox$dpo.imp, bait.id = "protein.29"))
})
