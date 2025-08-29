test_that("errors are returned if the object is not of class DEprot (prolfqua)", {
  expect_error(diff.analyses.prolfqua(DEprot.object = DEprot::test.toolbox$geneset,
                                      contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                           c("condition", "6h.10nM.E2", "6h.DMSO")),
                                      linear.FC.th = 1.2))
})



test_that("errors are returned if the object is not containing imputed data (prolfqua)", {
  expect_error(diff.analyses.prolfqua(DEprot.object = DEprot::test.toolbox$dpo.raw,
                                      contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                           c("condition", "6h.10nM.E2", "6h.DMSO")),
                                      linear.FC.th = 1.2))
})


##########################################################################################

test_that("no error is returned when performing differential analyses (prolfqua)", {
  expect_no_failure(diff.analyses.prolfqua(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                           contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                                c("condition", "6h.10nM.E2", "6h.DMSO")),
                                           linear.FC.th = 1.2))
})



