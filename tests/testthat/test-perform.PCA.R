## ----------------------------------------------------------------------------------------
##  perform.PCA()
##  Two engines are used depending on the data: stats::prcomp on complete matrices, and
##  pcaMethods::pca(method = "nipals") as soon as missing values are present. Only the
##  second one honours 'n.PCs'.
## ----------------------------------------------------------------------------------------

######################################    ERRORS    ######################################

test_that("errors are returned if the object is not of class DEprot or DEprot.analyses", {
  expect_error(perform.PCA(DEprot.object = DEprot::sample.config))
  expect_error(perform.PCA(DEprot.object = "not an object"))
})


test_that("counts that are not available cannot be used", {
  bare <- tb.dpo.norm
  bare@imputed.counts <- NULL

  expect_error(perform.PCA(DEprot.object = bare, which.data = "imputed"))
  expect_error(perform.PCA(DEprot.object = tb.dpo.imp, which.data = "wrong data"))
})


test_that("a sample subset that does not exist is rejected", {
  expect_error(perform.PCA(DEprot.object = tb.dpo.imp, sample.subset = c("not.a.sample")))
})



##########################################################################################

test_that("the PCA object carries the scores, the importance and the metadata", {
  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed")

  expect_s4_class(pca, "DEprot.PCA")
  expect_equal(pca@data.used, "imputed")
  expect_equal(nrow(pca@PCs), ncol(tb.dpo.imp@imputed.counts))
  expect_true(all(c("PC1", "PC2") %in% colnames(pca@PCs)))
  expect_s3_class(pca@importance, "data.frame")
  expect_true(inherits(pca@cumulative.PC.plot, "ggplot"))
  ## complete counts are handled by stats::prcomp
  expect_s3_class(pca@prcomp, "prcomp")
})


test_that("counts carrying missing values are handled by the nipals engine", {
  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "normalized")

  expect_s4_class(pca, "DEprot.PCA")
  ## pcaMethods returns an S4 'pcaRes' object instead of a 'prcomp' list
  expect_true(isS4(pca@prcomp) | inherits(pca@prcomp, "prcomp"))
})


test_that("the variance explained is decreasing and bounded", {
  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed")

  expect_true(all(c("Proportion.of.Variance", "Cumulative.Proportion") %in% colnames(pca@importance)))
  expect_true(all(pca@importance$Proportion.of.Variance >= 0 &
                    pca@importance$Proportion.of.Variance <= 1))
  ## the components are returned in decreasing order of explained variance
  expect_true(all(diff(pca@importance$Proportion.of.Variance) <= 1e-8))
  ## the cumulative proportion can only grow
  expect_true(all(diff(pca@importance$Cumulative.Proportion) >= -1e-8))
})


test_that("works on a subset of samples", {
  samples <- tb.dpo.norm@metadata$column.id[1:8]

  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed", sample.subset = samples)

  expect_equal(nrow(pca@PCs), length(samples))
  expect_equal(nrow(pca@PCA.metadata), length(samples))
  ## 'sample.subset' is stored as a single comma-separated string
  expect_type(pca@sample.subset, "character")
  expect_equal(sort(strsplit(pca@sample.subset, ", ")[[1]]), sort(samples))
})


test_that("every count type available can be used", {
  for (w in c("raw", "normalized", "imputed")) {
    pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = w)
    expect_s4_class(pca, "DEprot.PCA")
    expect_equal(pca@data.used, w)
  }
})


test_that("centering and scaling can be switched off", {
  scaled <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed",
                        center.data = TRUE, scale.data = TRUE)
  bare <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed",
                      center.data = FALSE, scale.data = FALSE)

  expect_false(isTRUE(all.equal(scaled@PCs$PC1, bare@PCs$PC1)))
})


test_that("'n.PCs' caps the components computed on incomplete counts", {
  ## the argument is used only by the nipals engine, hence on counts carrying NAs
  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "normalized", n.PCs = 3)

  expect_true(sum(grepl("^PC[0-9]+$", colnames(pca@PCs))) <= 3)
  expect_true(nrow(pca@importance) <= 3)
})


test_that("a DEprot.analyses object is accepted", {
  pca <- perform.PCA(DEprot.object = tb.limma, which.data = "imputed")

  expect_s4_class(pca, "DEprot.PCA")
})


test_that("the show method prints a summary", {
  pca <- perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed")

  expect_no_error(suppressWarnings(suppressMessages(show(pca))))
})
