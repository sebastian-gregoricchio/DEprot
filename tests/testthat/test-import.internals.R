## ----------------------------------------------------------------------------------------
##  internal_functions.import.R
##  The helpers shared by every reader: file reading, sample-name cleaning, reshaping,
##  MaxLFQ roll-up and the finalization step that builds the DEprot object.
## ----------------------------------------------------------------------------------------

####################################    .read.table.any    ###############################

test_that("a delimited file is read whatever its separator", {
  read.any <- DEprot:::.read.table.any
  dir <- local.tmpdir()

  tsv <- file.path(dir, "table.tsv")
  csv <- file.path(dir, "table.csv")
  tb <- data.frame(id = c("P1", "P2"), s1 = c(1, 2), s2 = c(3, 4), stringsAsFactors = FALSE)

  utils::write.table(tb, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.csv(tb, csv, row.names = FALSE, quote = FALSE)

  expect_equal(dim(read.any(tsv)), c(2, 3))
  expect_equal(dim(read.any(csv)), c(2, 3))
  ## the column names must not be mangled: they carry the sample IDs
  expect_equal(colnames(read.any(tsv)), c("id", "s1", "s2"))
})


test_that("a file that does not exist is reported", {
  expect_error(DEprot:::.read.table.any(file.path(tempdir(), "nope.tsv")))
})


test_that("a parquet report is read when a reader is installed", {
  skip_if_not_installed("nanoparquet")

  dir <- local.tmpdir()
  path <- file.path(dir, "table.parquet")
  tb <- data.frame(id = c("P1", "P2"), s1 = c(1, 2), stringsAsFactors = FALSE)
  nanoparquet::write_parquet(tb, path)

  out <- DEprot:::.read.table.any(path)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
})



###################################    .clean.run.names    ###############################

test_that("the run names are cleaned in every flavour", {
  clean <- DEprot:::.clean.run.names

  expect_equal(clean(c("D:/data/run1.raw", "/home/data/run2.RAW")), c("run1", "run2"))
  expect_equal(clean("[12] sample.d"), "sample")
  expect_equal(clean("sample.mzML"), "sample")
  ## nothing to clean
  expect_equal(clean("sample_01"), "sample_01")
  ## the function must be vectorized and keep the input length
  expect_length(clean(c("a.raw", "b.raw", "c.raw")), 3)
})



##################################    .pick.column    ####################################

test_that(".pick.column handles the three situations", {
  pick <- DEprot:::.pick.column
  cn <- c("Protein.Group", "Run", "PG.Quantity")

  expect_equal(pick("Run", "Protein.Group", cn), "Run")
  expect_equal(pick("auto", c("missing", "Run"), cn), "Run")
  expect_true(is.na(pick("auto", "missing", cn)))
  expect_error(pick("missing", "Run", cn, "sample.id"))
})



#############################    .long.to.matrix / .wide.to.matrix    ####################

test_that("the long reshaping aggregates the duplicated entries", {
  long <- data.frame(id = c("P1", "P1", "P2", "P2"),
                     sample = c("s1", "s1", "s1", "s2"),
                     value = c(10, 20, 5, 7),
                     stringsAsFactors = FALSE)

  summed <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                     quantity.col = "value",
                                     FUN = function(x) {sum(x, na.rm = TRUE)})

  expect_equal(summed["P1", "s1"], 30)
  ## a combination without any value stays NA
  expect_true(is.na(summed["P1", "s2"]))
})


test_that("the zeros are converted to NA unless requested otherwise", {
  long <- data.frame(id = c("P1", "P1", "P2"),
                     sample = c("s1", "s2", "s1"),
                     value = c(0, 8, 5),
                     stringsAsFactors = FALSE)

  with.na <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                      quantity.col = "value", zero.to.na = TRUE)
  kept <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                   quantity.col = "value", zero.to.na = FALSE)

  ## a zero becomes NA before the reshaping, hence the entry disappears from the long table
  ## and the cell is left empty; with zero.to.na = FALSE the value is kept as it is
  expect_true(is.na(with.na["P1", "s1"]))
  expect_equal(unname(kept["P1", "s1"]), 0)
  expect_equal(unname(with.na["P1", "s2"]), 8)
})


test_that("a column without any finite value is reported", {
  long <- data.frame(id = c("P1", "P2"),
                     sample = c("s1", "s2"),
                     value = c(NA, NA),
                     stringsAsFactors = FALSE)

  expect_error(DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                        quantity.col = "value"))
})


test_that("the wide reshaping renames the samples", {
  wide <- data.frame(id = c("P1", "P2"),
                     `LFQ intensity A` = c(1, 2),
                     `LFQ intensity B` = c(3, 4),
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

  cols <- c(A = "LFQ intensity A", B = "LFQ intensity B")
  m <- DEprot:::.wide.to.matrix(wide, "id", cols)

  expect_equal(colnames(m), c("A", "B"))
  expect_equal(rownames(m), c("P1", "P2"))
  expect_true(is.numeric(m))
})



####################################    .maxlfq.rollup    ################################

test_that("the MaxLFQ roll-up returns a protein-by-sample matrix", {
  skip_if_not_installed("iq")

  long <- expand.grid(protein = c("P1", "P2"),
                      sample = c("s1", "s2"),
                      precursor = c("pep1", "pep2"),
                      stringsAsFactors = FALSE)
  long$quantity <- seq_len(nrow(long)) * 100

  m <- DEprot:::.maxlfq.rollup(tb = long,
                               primary.id = "protein",
                               secondary.id = "precursor",
                               sample.col = "sample",
                               quantity.col = "quantity",
                               install.missing = "never",
                               verbose = FALSE)

  expect_true(is.matrix(m))
  expect_equal(sort(rownames(m)), c("P1", "P2"))
  expect_equal(sort(colnames(m)), c("s1", "s2"))
})



###################################    .finalize.import    ###############################

test_that("duplicated sample names are made unique", {
  imported <- list(counts = matrix(1:4, nrow = 2,
                                   dimnames = list(c("P1", "P2"), c("D:/run.raw", "E:/run.raw"))),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = NA,
                   quantity = "test",
                   feature = "proteins")

  ## the cleaning would collapse the two runs into the same name: it must be skipped
  out <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, verbose = FALSE)))

  expect_s4_class(out, "DEprot")
  expect_equal(ncol(any.counts(out)), 2)
  expect_false(any(duplicated(colnames(any.counts(out)))))
})


test_that("the sample names can be left untouched", {
  imported <- list(counts = matrix(1:4, nrow = 2,
                                   dimnames = list(c("P1", "P2"), c("D:/a.raw", "D:/b.raw"))),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = NA,
                   quantity = "test",
                   feature = "proteins")

  cleaned <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, clean.sample.names = TRUE, verbose = FALSE)))
  raw.names <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, clean.sample.names = FALSE, verbose = FALSE)))

  expect_equal(sort(colnames(any.counts(cleaned))), c("a", "b"))
  expect_true(all(grepl("\\.raw$", colnames(any.counts(raw.names)))))
})


test_that("an empty counts matrix is reported", {
  imported <- list(counts = matrix(numeric(0), nrow = 0, ncol = 0),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = NA,
                   quantity = "test",
                   feature = "proteins")

  expect_error(DEprot:::.finalize.import(imported = imported, verbose = FALSE))
  expect_error(DEprot:::.finalize.import(imported = list(counts = NULL), verbose = FALSE))
})


test_that("a minimal metadata is generated when none is provided", {
  imported <- list(counts = matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b"))),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = NA,
                   quantity = "test",
                   feature = "proteins")

  expect_warning(DEprot:::.finalize.import(imported = imported, verbose = FALSE))

  out <- suppressWarnings(DEprot:::.finalize.import(imported = imported, verbose = FALSE))
  expect_equal(nrow(out@metadata), 2)
  expect_true("column.id" %in% colnames(out@metadata))
})


test_that("the annotation carried by the imported object becomes the metadata", {
  imported <- list(counts = matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b"))),
                   log.base = 1,
                   data.type = "normalized",
                   normalization.method = "test",
                   quantity = "test",
                   feature = "proteins",
                   annotation = data.frame(column.id = c("a", "b"),
                                           condition = c("ctrl", "treat"),
                                           stringsAsFactors = FALSE))

  out <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, verbose = FALSE)))

  expect_true("condition" %in% colnames(out@metadata))
  expect_equal(nrow(out@metadata), 2)
})


test_that("the metadata can be given as a file path", {
  dir <- local.tmpdir()
  path <- file.path(dir, "meta.tsv")
  utils::write.table(data.frame(column.id = c("a", "b"), condition = c("x", "y")),
                     path, sep = "\t", row.names = FALSE, quote = FALSE)

  imported <- list(counts = matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b"))),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = NA,
                   quantity = "test",
                   feature = "proteins")

  out <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, metadata = path, verbose = FALSE)))

  expect_equal(nrow(out@metadata), 2)
  expect_true("condition" %in% colnames(out@metadata))
})


test_that("the data type declared overrides the one detected", {
  imported <- list(counts = matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b"))),
                   log.base = 1,
                   data.type = "raw",
                   normalization.method = "test",
                   quantity = "test",
                   feature = "proteins")

  meta <- data.frame(column.id = c("a", "b"), stringsAsFactors = FALSE)

  as.raw <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, metadata = meta, verbose = FALSE)))
  as.norm <- suppressWarnings(suppressMessages(
    DEprot:::.finalize.import(imported = imported, metadata = meta,
                              data.type = "normalized", verbose = FALSE)))

  expect_false(as.raw@normalized)
  expect_true(as.norm@normalized)
})
