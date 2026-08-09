## ----------------------------------------------------------------------------------------
##  load.counts()
##  Deprecated predecessor of load.counts2(): it is still exported, hence still tested.
##  Every call raises the deprecation warning, which is silenced where it is not the point
##  of the test.
## ----------------------------------------------------------------------------------------

counts.mat <- make.counts(n.prot = 20, n.samples = 6, n.missing = 0)
meta.tb <- make.metadata(counts.mat)


######################################    ERRORS    ######################################

test_that("the function warns that it is deprecated", {
  expect_warning(load.counts(counts = counts.mat, metadata = meta.tb, log.base = 2),
                 regexp = "deprecated")
})


test_that("counts that are not a table are rejected", {
  expect_error(suppressWarnings(load.counts(counts = "not a table",
                                            metadata = meta.tb,
                                            log.base = 2)))
  expect_error(suppressWarnings(load.counts(counts = 1:10,
                                            metadata = meta.tb,
                                            log.base = 2)))
})


test_that("empty or missing protein IDs are rejected", {
  bad <- counts.mat
  rownames(bad)[1] <- ""

  expect_error(suppressWarnings(load.counts(counts = bad, metadata = meta.tb, log.base = 2)))
})


test_that("duplicated protein IDs are made unique with a warning", {
  dup <- counts.mat
  rownames(dup)[2] <- rownames(dup)[1]

  expect_warning(load.counts(counts = dup, metadata = meta.tb, log.base = 2),
                 regexp = "duplicated")

  dpo <- suppressWarnings(load.counts(counts = dup, metadata = meta.tb, log.base = 2))
  expect_false(any(duplicated(rownames(any.counts(dpo)))))
})


test_that("a metadata of the wrong type is rejected", {
  expect_error(suppressWarnings(load.counts(counts = counts.mat,
                                            metadata = 1:6,
                                            log.base = 2)))
})


test_that("a metadata without the 'column.id' column is rejected", {
  bad.meta <- meta.tb
  colnames(bad.meta)[1] <- "sample"

  expect_error(suppressWarnings(load.counts(counts = counts.mat,
                                            metadata = bad.meta,
                                            log.base = 2)))
})


test_that("sample names that do not match the metadata are rejected", {
  bad.meta <- meta.tb
  bad.meta$column.id[1] <- "not.a.sample"

  expect_error(suppressWarnings(load.counts(counts = counts.mat,
                                            metadata = bad.meta,
                                            log.base = 2)))
})



##########################################################################################

test_that("raw counts are loaded when no normalization and no imputation are declared", {
  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = meta.tb,
                                      log.base = 2,
                                      normalization.method = NULL))

  expect_s4_class(dpo, "DEprot")
  expect_equal(dim(dpo@raw.counts), dim(counts.mat))
  expect_false(dpo@normalized)
  expect_false(dpo@imputed)
  expect_equal(dpo@normalization.method, "none")
  expect_true(inherits(dpo@boxplot.raw, "ggplot"))
})


test_that("normalized counts are loaded in the dedicated slot", {
  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = meta.tb,
                                      log.base = 2,
                                      normalization.method = "quantile"))

  expect_true(dpo@normalized)
  expect_equal(dim(dpo@norm.counts), dim(counts.mat))
  expect_true(is.null(dpo@raw.counts))
  expect_equal(dpo@normalization.method, "quantile")
  expect_true(inherits(dpo@boxplot.norm, "ggplot"))
})


test_that("imputed counts are loaded in the dedicated slot", {
  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = meta.tb,
                                      log.base = 2,
                                      imputation = "missForest",
                                      normalization.method = NULL))

  expect_true(dpo@imputed)
  expect_equal(dpo@imputation.method, "missForest")
  expect_equal(dim(dpo@imputed.counts), dim(counts.mat))
  expect_true(inherits(dpo@boxplot.imputed, "ggplot"))
})


test_that("the log base is recorded and cannot be empty", {
  log2.obj <- suppressWarnings(load.counts(counts = counts.mat, metadata = meta.tb, log.base = 2))
  log10.obj <- suppressWarnings(load.counts(counts = counts.mat, metadata = meta.tb, log.base = 10))

  expect_equal(log2.obj@log.base, 2)
  expect_equal(log10.obj@log.base, 10)
  expect_true(log2.obj@log.transformed)

  ## the slot is declared as 'numeric' in the class: a NULL cannot be stored
  expect_error(suppressWarnings(load.counts(counts = counts.mat,
                                            metadata = meta.tb,
                                            log.base = NULL)))
})


test_that("a data.frame is accepted as counts table", {
  dpo <- suppressWarnings(load.counts(counts = as.data.frame(counts.mat),
                                      metadata = meta.tb,
                                      log.base = 2))

  expect_s4_class(dpo, "DEprot")
  expect_equal(nrow(any.counts(dpo)), nrow(counts.mat))
})


test_that("the metadata can be provided as a file path", {
  dir <- local.tmpdir()
  path <- file.path(dir, "metadata.tsv")
  utils::write.table(meta.tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = path,
                                      log.base = 2))

  expect_equal(nrow(dpo@metadata), nrow(meta.tb))
})


test_that("a column other than 'column.id' can carry the sample IDs", {
  renamed <- meta.tb
  colnames(renamed)[colnames(renamed) == "column.id"] <- "sample.name"

  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = renamed,
                                      log.base = 2,
                                      column.id = "sample.name"))

  ## the column is renamed internally, so that the object stays consistent
  expect_true("column.id" %in% colnames(dpo@metadata))
  expect_equal(sort(dpo@metadata$column.id), sort(colnames(counts.mat)))
})


test_that("a protein.info table is attached and aligned", {
  info <- data.frame(prot.id = rownames(counts.mat),
                     gene.name = toupper(rownames(counts.mat)),
                     stringsAsFactors = FALSE)

  dpo <- suppressWarnings(load.counts(counts = counts.mat,
                                      metadata = meta.tb,
                                      log.base = 2,
                                      protein.info = info,
                                      protein.info.id.column = "prot.id"))

  expect_false(is.null(dpo@protein.info))
  expect_equal(nrow(dpo@protein.info), nrow(counts.mat))
})


test_that("the object built is usable by the rest of the package", {
  dpo <- suppressWarnings(load.counts(counts = counts.mat, metadata = meta.tb, log.base = 2))

  expect_no_error(suppressWarnings(suppressMessages(show(dpo))))
  expect_equal(nrow(get.metadata(dpo)), nrow(meta.tb))
})
