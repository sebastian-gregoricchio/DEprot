##########################################################################################

test_that("the function geneset.enrichment is working for GSEA analyses with rank based on corrlation", {
  expect_no_error(suppressMessages(suppressWarnings(geneset.enrichment(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
                                                                       contrast = 1,
                                                                       TERM2GENE = DEprot::test.toolbox$geneset,
                                                                       enrichment.type = "GSEA",
                                                                       gsea.rank.method = "corrlation",
                                                                       pvalueCutoff = 1,
                                                                       qvalueCutoff = 1))))

})



test_that("the function geneset.enrichment is working for GSEA analyses with rank based on foldchange", {
  expect_no_error(invisible(geneset.enrichment(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
                                               contrast = 1,
                                               TERM2GENE = DEprot::test.toolbox$geneset,
                                               enrichment.type = "GSEA",
                                               gsea.rank.method = "foldchange",
                                               pvalueCutoff = 1,
                                               qvalueCutoff = 1)))

})



test_that("the function geneset.enrichment is working for ORA analyses", {
  expect_no_error(invisible(geneset.enrichment(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
                                               contrast = 1,
                                               TERM2GENE = DEprot::test.toolbox$geneset,
                                               enrichment.type = "ORA",
                                               diff.status.category = "FBS",
                                               pvalueCutoff = 1,
                                               qvalueCutoff = 1)))

})




test_that("the function simplify is working for ORA analyses", {
  expect_no_failure(suppressMessages(simplify.enrichment(enrichment.results = DEprot::test.toolbox$ora.results)))
})



test_that("the function simplify.enrichment is working for GSEA analyses", {
  expect_no_failure(suppressMessages(simplify.enrichment(enrichment.results = DEprot::test.toolbox$gsea.results)))
})


test_that("the function simplify.enrichment is not working when the input is not of class DEprot.enrichResult", {
  expect_error(simplify.enrichment(enrichment.results = DEprot::test.toolbox$geneset))

})




test_that("the function NES.plot is working for GSEA analyses", {
  expect_no_error(NES.plot(enrichResult = DEprot::test.toolbox$gsea.results))
})


test_that("the function NES.plot is not working for ORA analyses", {
  expect_error(NES.plot(enrichResult = DEprot::test.toolbox$ora.results))
})


test_that("the function plot.GSEA is working for GSEA analyses", {
  expect_no_error(plot.GSEA(gsea.results = DEprot::test.toolbox$gsea.results, geneset.id = 1))
})

test_that("the function plot.GSEA is not working for ORA analyses", {
  expect_error(plot.GSEA(gsea.results = DEprot::test.toolbox$ora.results, geneset.id = 1))
})









