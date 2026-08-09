## ----------------------------------------------------------------------------------------
##  import.external(): FragPipe, Proteome Discoverer and long Spectronaut reports
##  These three readers are the ones left uncovered by 'test-import.external.R'.
## ----------------------------------------------------------------------------------------

######################################    FragPipe    ####################################

test_that("a FragPipe combined_protein report is imported", {
  dir <- local.tmpdir()
  path <- write.fragpipe(dir)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  ## the default quantity is 'MaxLFQ Intensity', hence the data are normalized
  expect_true(dpo@normalized)
  expect_equal(sort(colnames(dpo@norm.counts)), c("sampleA", "sampleB"))
  ## the contaminant entry is dropped on the protein ID
  expect_equal(nrow(dpo@norm.counts), 3)
})


test_that("the shorter FragPipe suffixes do not swallow the longer ones", {
  dir <- local.tmpdir()
  path <- write.fragpipe(dir)

  ## 'Intensity' must not collect the 'MaxLFQ Intensity' columns
  dpo <- suppressWarnings(import.external(file = path,
                                          quantity = "Intensity",
                                          verbose = FALSE))

  cnt <- any.counts(dpo)

  expect_equal(sort(colnames(cnt)), c("sampleA", "sampleB"))
  ## the plain intensities are one order of magnitude below the MaxLFQ ones in the fixture
  expect_true(max(cnt, na.rm = TRUE) < 1000)
  expect_false(dpo@normalized)
})


test_that("the spectral counts can be imported as well", {
  dir <- local.tmpdir()
  path <- write.fragpipe(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          quantity = "Spectral Count",
                                          remove.contaminants = FALSE,
                                          verbose = FALSE))

  expect_equal(nrow(any.counts(dpo)), 4)
})


test_that("an unknown FragPipe suffix lists the ones available", {
  dir <- local.tmpdir()
  path <- write.fragpipe(dir)

  expect_error(import.external(file = path, quantity = "not.a.suffix", verbose = FALSE))
})


test_that("the FragPipe contaminant pattern can be customized", {
  dir <- local.tmpdir()
  path <- write.fragpipe(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          contaminant.pattern = "^P00003",
                                          verbose = FALSE))

  expect_false("P00003" %in% rownames(any.counts(dpo)))
  expect_true("contam_P99999" %in% rownames(any.counts(dpo)))
})



##############################    Proteome Discoverer    #################################

test_that("a Proteome Discoverer export is imported", {
  dir <- local.tmpdir()
  path <- write.proteome.discoverer(dir)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  cnt <- any.counts(dpo)

  expect_s4_class(dpo, "DEprot")
  ## 'Abundances (Normalized)' is the default quantity: the counts are normalized
  expect_true(dpo@normalized)
  expect_equal(ncol(cnt), 2)
  ## both the ID pattern and the 'Contaminant' column drop the last entry
  expect_equal(nrow(cnt), 3)
})


test_that("an unknown quantity falls back on the abundance columns", {
  dir <- local.tmpdir()
  path <- write.proteome.discoverer(dir)

  ## no column starts with the prefix requested: the reader falls back on 'Abundance'
  dpo <- suppressWarnings(import.external(file = path,
                                          quantity = "not.a.quantity",
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
  expect_equal(ncol(any.counts(dpo)), 2)
})


test_that("the abundance columns are found through the fallback pattern", {
  path <- file.path(local.tmpdir(), "PD.raw.txt")

  tb <- data.frame(Accession = c("P00001", "P00002"),
                   `Abundance: F1: Sample, WT` = c(100, 200),
                   `Abundance: F2: Sample, KO` = c(150, 250),
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
  utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  dpo <- suppressWarnings(import.external(file = path, verbose = FALSE))

  ## a plain 'Abundance' is not normalized
  expect_false(dpo@normalized)
  expect_equal(dim(dpo@raw.counts), c(2, 2))
})


test_that("a Proteome Discoverer export without abundances is rejected", {
  dir <- local.tmpdir()
  path <- file.path(dir, "PD.empty.txt")

  ## no abundance column at all: the source must be declared, since the automatic detection
  ## relies precisely on their presence
  tb <- data.frame(Accession = c("P00001", "P00002"),
                   Description = c("a", "b"),
                   Coverage = c(10, 20),
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
  utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_error(import.external(file = path,
                               source = "proteome.discoverer",
                               verbose = FALSE))
  ## without the source, the file cannot even be recognized
  expect_error(import.external(file = path, verbose = FALSE))
})



###############################    Spectronaut (long)    #################################

test_that("a long Spectronaut report is imported at the protein level", {
  dir <- local.tmpdir()
  path <- write.spectronaut.long(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          summarization = "none",
                                          verbose = FALSE))

  cnt <- any.counts(dpo)

  expect_s4_class(dpo, "DEprot")
  expect_equal(nrow(cnt), 3)
  expect_equal(ncol(cnt), 2)
  ## the '[1] ' prefixes and the '.raw' extensions are cleaned away
  expect_equal(sort(colnames(cnt)), c("sampleA", "sampleB"))
})


test_that("the Spectronaut fragment quantities can be rolled up", {
  dir <- local.tmpdir()
  path <- write.spectronaut.long(dir)

  summed <- suppressWarnings(import.external(file = path, summarization = "sum", verbose = FALSE))
  median <- suppressWarnings(import.external(file = path, summarization = "median", verbose = FALSE))

  expect_equal(dim(any.counts(summed)), dim(any.counts(median)))
  ## the sum of two precursors is above their median
  expect_true(all(any.counts(summed) >= any.counts(median), na.rm = TRUE))
  expect_false(summed@normalized)
})


test_that("the q-value filters are applied to the long report", {
  path <- file.path(local.tmpdir(), "spectronaut.qvalues.tsv")

  ## three proteins, one of which fails the precursor q-value: two of them survive, which
  ## keeps the counts matrix two-dimensional whatever the filter applied
  tb <- expand.grid(R.FileName = c("sampleA", "sampleB"),
                    PG.ProteinGroups = c("P00001", "P00002", "P00003"),
                    stringsAsFactors = FALSE)
  tb$EG.PrecursorId <- "PEPTIDEK.2"
  tb$EG.Qvalue <- ifelse(tb$PG.ProteinGroups == "P00003", 0.5, 0.001)
  tb$PG.Qvalue <- 0.001
  tb$PG.Quantity <- 1000
  tb$FG.Quantity <- 500
  utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  filtered <- suppressWarnings(import.external(file = path,
                                               summarization = "none",
                                               q.value = 0.01,
                                               verbose = FALSE))

  kept <- suppressWarnings(import.external(file = path,
                                           summarization = "none",
                                           q.value = NULL,
                                           pg.q.value = NULL,
                                           verbose = FALSE))

  expect_equal(nrow(any.counts(filtered)), 2)
  expect_equal(nrow(any.counts(kept)), 3)
})


test_that("a Spectronaut report without a protein-level quantity is reported", {
  dir <- local.tmpdir()
  path <- file.path(dir, "spectronaut.nopg.tsv")

  tb <- data.frame(R.FileName = c("sampleA", "sampleB"),
                   PG.ProteinGroups = c("P00001", "P00001"),
                   FG.Quantity = c(100, 200),
                   stringsAsFactors = FALSE)
  utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_error(import.external(file = path, summarization = "none", verbose = FALSE))
})


test_that("the source can be forced instead of being detected", {
  dir <- local.tmpdir()
  path <- write.spectronaut.long(dir)

  dpo <- suppressWarnings(import.external(file = path,
                                          source = "spectronaut",
                                          summarization = "none",
                                          verbose = FALSE))

  expect_s4_class(dpo, "DEprot")
})


test_that("a report that does not match the source declared is reported", {
  dir <- local.tmpdir()
  path <- write.maxquant(dir)

  expect_error(import.external(file = path, source = "spectronaut", verbose = FALSE))
  expect_error(import.external(file = path, source = "diann", verbose = FALSE))
})
