## Object used by the structural tests: the model is fitted once and reused.
## The counts must still contain the missing values, hence the normalized (not imputed) object.
dpo.proDA =
  suppressWarnings(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                       contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                            c("condition", "6h.10nM.E2", "6h.DMSO")),
                                       check.missingness = FALSE,
                                       linear.FC.th = 1.2))



test_that("errors are returned if the object is not of class DEprot (proDA)", {
  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$geneset,
                                   contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                        c("condition", "6h.10nM.E2", "6h.DMSO")),
                                   linear.FC.th = 1.2))
})



test_that("errors are returned if the object is not containing normalized data (proDA)", {
  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.raw,
                                   contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                   which.data = "normalized",
                                   linear.FC.th = 1.2))
})



test_that("errors are returned if counts without missing values are required (proDA)", {
  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.imp,
                                   contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                   which.data = "imputed",
                                   linear.FC.th = 1.2))

  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.random,
                                   contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                   which.data = "randomized",
                                   linear.FC.th = 1.2))
})



test_that("errors are returned if the replicate model is required without a replicate column (proDA)", {
  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                   contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                   include.rep.model = TRUE,
                                   linear.FC.th = 1.2))
})



test_that("errors are returned if the contrast is not available in the metadata (proDA)", {
  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                   contrast.list = list(c("not.a.column", "FBS", "6h.DMSO")),
                                   linear.FC.th = 1.2))

  expect_error(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                   contrast.list = list(c("condition", "not.a.condition", "6h.DMSO")),
                                   linear.FC.th = 1.2))
})


##########################################################################################

test_that("no error is returned when performing differential analyses (proDA)", {
  expect_no_error(suppressWarnings(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                                       contrast.list = list(c("condition", "FBS", "6h.DMSO"),
                                                                            c("condition", "6h.10nM.E2", "6h.DMSO")),
                                                       check.missingness = FALSE,
                                                       linear.FC.th = 1.2)))
})



test_that("no error is returned when including the replicate in the model (proDA)", {
  expect_no_error(suppressWarnings(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                                       contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                       include.rep.model = TRUE,
                                                       replicate.column = "replicate",
                                                       check.missingness = FALSE,
                                                       linear.FC.th = 1.2)))
})



test_that("no error is returned when checking the missingness pattern (proDA)", {
  expect_no_error(suppressWarnings(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                                       contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                       check.missingness = TRUE,
                                                       linear.FC.th = 1.2)))
})



test_that("no error is returned when subsampling the proteins used for the hyper-parameters (proDA)", {
  expect_no_error(suppressWarnings(diff.analyses.proDA(DEprot.object = DEprot::test.toolbox$dpo.norm,
                                                       contrast.list = list(c("condition", "FBS", "6h.DMSO")),
                                                       n.subsample = 25,
                                                       seed = 1234,
                                                       check.missingness = FALSE,
                                                       linear.FC.th = 1.2)))
})



test_that("no error is returned when reapplying differential analyses thresholds (proDA)", {
  expect_no_error(reapply.thresholds(DEprot.analyses.object = dpo.proDA, linear.FC = 1.45))
})


##########################################################################################

test_that("the object returned is a DEprot.analyses generated by proDA", {
  expect_s4_class(dpo.proDA, "DEprot.analyses")
  expect_equal(dpo.proDA@differential.analyses.params$stat.test, "proDA")
  expect_equal(dpo.proDA@differential.analyses.params$counts.used, "normalized")
  expect_equal(length(dpo.proDA@analyses.result.list), 2)
})



test_that("the results table keeps one row per protein and contains the proDA-specific columns", {
  res = dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$results

  expect_equal(nrow(res), nrow(DEprot::test.toolbox$dpo.norm@norm.counts))
  expect_true(all(c("n.detected.FBS", "n.detected.6h.DMSO", "n.approx") %in% colnames(res)))
  expect_true(all(c("log2.Fold_FBS.vs.6h.DMSO", "lfcSE", "diff.status") %in% colnames(res)))
  expect_false(TRUE %in% is.na(res$diff.status))
})



test_that("the number of detected values corresponds to the counts table", {
  res = dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$results
  counts = DEprot::test.toolbox$dpo.norm@norm.counts
  samples = dpo.proDA@contrasts$condition_FBS.vs.6h.DMSO$group.1

  expect_equal(res$n.detected.FBS,
               unname(rowSums(!is.na(counts[,samples, drop = FALSE]))[match(res$prot.id, rownames(counts))]))
})



test_that("the statistic corresponds to the ratio between log2(FoldChange) and lfcSE", {
  res = dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$results
  res = res[!is.na(res$statistic),]

  expect_equal(res$statistic,
               res$log2.Fold_FBS.vs.6h.DMSO / res$lfcSE,
               tolerance = 1e-6)
})



test_that("a dropout curve is estimated for each sample", {
  curves = dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$proDA.fit$dropout.curves

  expect_equal(nrow(curves), ncol(DEprot::test.toolbox$dpo.norm@norm.counts))
  expect_false(TRUE %in% is.na(curves$position))
})



test_that("contrasts sharing the same metadata column share a single fit", {
  expect_identical(dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$proDA.fit$fit,
                   dpo.proDA@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$proDA.fit$fit)
})



test_that("the contrast vector is built on the groups of the contrast", {
  cntrst = dpo.proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$proDA.fit$contrast

  expect_equal(sum(cntrst), 0)
  expect_equal(unname(cntrst[make.names("FBS")]), 1)
  expect_equal(unname(cntrst[make.names("6h.DMSO")]), -1)
})
