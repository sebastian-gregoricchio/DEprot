## ----------------------------------------------------------------------------------------
##  import.external()
##  The fixtures are written on the fly: a handful of proteins is enough to walk through
##  every reader, the point being the parsing and not the quantification.
## ----------------------------------------------------------------------------------------

## 'load.counts2' fills a single counts slot, the one matching the detected data type:
## the counts of a DIA-NN or Spectronaut report end up in 'norm.counts', those of a
## MaxQuant 'Intensity' column in 'raw.counts'. The tests read whichever slot was filled.
imported.counts <-
  function(DEprot.object) {
    for (slot in c("imputed.counts", "norm.counts", "random.counts", "raw.counts")) {
      cnt <- methods::slot(DEprot.object, slot)
      if (!is.null(cnt)) {return(cnt)}
    }
    return(NULL)
  }



test_that("a DIA-NN long report is imported and the source is detected automatically", {
  dir <- local.tmpdir()
  path <- write.diann.report(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          summarization = "sum",
                                          verbose = FALSE))

  cnt <- imported.counts(dpo)

  expect_s4_class(dpo, "DEprot")
  expect_equal(ncol(cnt), 2)
  ## the contaminant protein must have been dropped
  expect_false(any(grepl("CON__", rownames(cnt))))
})


test_that("the DIA-NN protein-level quantity can be used instead of the roll-up", {
  dir <- local.tmpdir()
  path <- write.diann.report(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          summarization = "none",
                                          verbose = FALSE))

  cnt <- imported.counts(dpo)

  expect_s4_class(dpo, "DEprot")
  expect_true(nrow(cnt) > 0)
  expect_true(dpo@normalized)
})


test_that("the median roll-up is accepted as well", {
  dir <- local.tmpdir()
  path <- write.diann.report(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          summarization = "median",
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(ncol(imported.counts(dpo)), 2)
})


test_that("MaxLFQ requires 'iq' and is skipped when the package is missing", {
  skip_if_not_installed("iq")

  dir <- local.tmpdir()
  path <- write.diann.report(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          summarization = "maxlfq",
                                          install.missing = "never",
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(dpo@log.base, 2)
})


test_that("a DIA-NN matrix is imported and the sample names are cleaned", {
  dir <- local.tmpdir()
  path <- write.diann.matrix(dir)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  cnt <- imported.counts(dpo)

  expect_s4_class(dpo, "DEprot")
  expect_equal(nrow(cnt), 3)
  ## paths and extensions are removed by .clean.run.names()
  expect_false(any(grepl("\\.raw$|/", colnames(cnt))))
})


test_that("a MaxQuant proteinGroups table is imported without contaminants", {
  dir <- local.tmpdir()
  path <- write.maxquant(dir)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  ## the reverse hit and the potential contaminant are flagged in dedicated columns
  expect_equal(nrow(dpo@norm.counts), 2)
  expect_equal(sort(colnames(dpo@norm.counts)), c("sampleA", "sampleB"))
})


test_that("the MaxQuant quantity column can be chosen", {
  dir <- local.tmpdir()
  path <- write.maxquant(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          quantity = "Intensity",
                                          remove.contaminants = FALSE,
                                          verbose = FALSE))

  ## a plain 'Intensity' column is not normalized: the counts land in the raw slot
  expect_equal(nrow(dpo@raw.counts), 4)
  expect_false(dpo@normalized)
})


test_that("a Spectronaut pivot report is imported", {
  dir <- local.tmpdir()
  path <- write.spectronaut.pivot(dir)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(nrow(imported.counts(dpo)), 3)
})


test_that("a generic table is imported when the columns are declared", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          source = "generic",
                                          id.column = "id",
                                          sample.columns = c("ctrl.1", "ctrl.2", "treat.1", "treat.2"),
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(dim(dpo@raw.counts), c(3, 4))
})


test_that("the data type can be forced by the user", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          source = "generic",
                                          id.column = "id",
                                          sample.columns = c("ctrl.1", "ctrl.2"),
                                          data.type = "normalized",
                                          verbose = FALSE))

  expect_true(dpo@normalized)
  expect_true(is.null(dpo@raw.counts))
})


test_that("a metadata table is attached and the counts are subset to it", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  meta <- data.frame(column.id = c("ctrl.1", "treat.1"),
                     condition = c("ctrl", "treat"),
                     stringsAsFactors = FALSE)

  dpo <- suppressWarnings(import.external(file = path,
                                          metadata = meta,
                                          source = "generic",
                                          id.column = "id",
                                          sample.columns = c("ctrl.1", "ctrl.2", "treat.1", "treat.2"),
                                          subset.to.metadata = TRUE,
                                          verbose = FALSE))

  expect_equal(ncol(dpo@raw.counts), 2)
  expect_equal(nrow(dpo@metadata), 2)
})


test_that("samples absent from the metadata stop the import when they must be kept", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  meta <- data.frame(column.id = "ctrl.1", condition = "ctrl", stringsAsFactors = FALSE)

  expect_error(import.external(file = path,
                               metadata = meta,
                               source = "generic",
                               id.column = "id",
                               sample.columns = c("ctrl.1", "ctrl.2"),
                               subset.to.metadata = FALSE,
                               verbose = FALSE))
})


test_that("a metadata listing samples absent from the data is reported", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  meta <- data.frame(column.id = c("ctrl.1", "not.a.sample"),
                     condition = c("ctrl", "treat"),
                     stringsAsFactors = FALSE)

  expect_error(import.external(file = path,
                               metadata = meta,
                               source = "generic",
                               id.column = "id",
                               sample.columns = c("ctrl.1", "ctrl.2"),
                               verbose = FALSE))
})



## ----------------------------------------------------------------------------------------
##  Error branches
## ----------------------------------------------------------------------------------------

test_that("a missing file is reported", {
  expect_error(import.external(file = file.path(tempdir(), "does.not.exist.tsv"), verbose = FALSE))
})


test_that("an unrecognizable table cannot be imported automatically", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  expect_error(import.external(file = path, verbose = FALSE))
})


test_that("the generic reader requires both 'id.column' and 'sample.columns'", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  expect_error(import.external(file = path, source = "generic", id.column = "id", verbose = FALSE))
  expect_error(import.external(file = path,
                               source = "generic",
                               id.column = "not.a.column",
                               sample.columns = "ctrl.1",
                               verbose = FALSE))
})


test_that("a quantity column that does not exist is reported", {
  dir <- local.tmpdir()
  path <- write.maxquant(dir)

  expect_error(import.external(file = path, quantity = "not.a.quantity", verbose = FALSE))
})


test_that("a metadata table without the 'column.id' column is rejected", {
  dir <- local.tmpdir()
  path <- write.generic(dir)

  expect_error(import.external(file = path,
                               metadata = data.frame(wrong.id = "ctrl.1"),
                               source = "generic",
                               id.column = "id",
                               sample.columns = c("ctrl.1", "ctrl.2"),
                               verbose = FALSE))
})



## ----------------------------------------------------------------------------------------
##  Internal helpers
## ----------------------------------------------------------------------------------------

test_that(".detect.source recognizes every supported flavour", {
  detect <- DEprot:::.detect.source

  expect_equal(detect(data.frame(Protein.Group = "P1", Run = "r1")), "diann")
  expect_equal(detect(data.frame(Protein.Group = "P1", sampleA = 1)), "diann.matrix")
  expect_equal(detect(data.frame(R.FileName = "r1", PG.Quantity = 1)), "spectronaut")
  expect_equal(detect(data.frame(PG.ProteinGroups = "P1", PG.Quantity = 1)), "spectronaut.pivot")
  expect_equal(detect(data.frame(`Majority protein IDs` = "P1", check.names = FALSE)), "maxquant")
  expect_equal(detect(data.frame(`Protein ID` = "P1", `Entry Name` = "N", check.names = FALSE)), "fragpipe")
  expect_equal(detect(data.frame(Accession = "P1", Abundance.F1 = 1)), "proteome.discoverer")
  expect_error(detect(data.frame(a = 1, b = 2)))
})


test_that(".clean.run.names strips paths, extensions and Spectronaut prefixes", {
  clean <- DEprot:::.clean.run.names

  expect_equal(clean("D:/data/20240101_sampleA.raw"), "20240101_sampleA")
  expect_equal(clean("[1] sampleB.raw"), "sampleB")
  expect_equal(clean("plain.name"), "plain.name")
})


test_that(".pick.column falls back on the candidates and reports the failures", {
  pick <- DEprot:::.pick.column
  cn <- c("Protein.Group", "Run", "Precursor.Normalised")

  ## explicit choice
  expect_equal(pick("Run", c("Protein.Group"), cn), "Run")
  ## automatic choice: the first available candidate wins
  expect_equal(pick("auto", c("Protein.Ids", "Protein.Group"), cn), "Protein.Group")
  ## nothing available
  expect_true(is.na(pick("auto", c("not.here"), cn)))
  ## an explicit column that does not exist is an error when the argument is named
  expect_error(pick("not.here", c("Protein.Group"), cn, "protein.id"))
})


test_that(".long.to.matrix and .wide.to.matrix build consistent matrices", {
  long <- data.frame(id = rep(c("P1", "P2"), each = 2),
                     sample = rep(c("s1", "s2"), times = 2),
                     value = c(1, 2, 3, 4),
                     stringsAsFactors = FALSE)

  m <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                quantity.col = "value", FUN = function(x) {mean(x, na.rm = TRUE)})

  expect_equal(dim(m), c(2, 2))
  expect_equal(sort(rownames(m)), c("P1", "P2"))

  wide <- data.frame(id = c("P1", "P2"), s1 = c(1, 3), s2 = c(2, 4), stringsAsFactors = FALSE)
  m2 <- DEprot:::.wide.to.matrix(wide, "id", c(s1 = "s1", s2 = "s2"))

  expect_equal(m[sort(rownames(m)), sort(colnames(m))],
               m2[sort(rownames(m2)), sort(colnames(m2))])
})
