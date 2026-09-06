## ----------------------------------------------------------------------------------------
##  internal_functions.import.R: the branches left uncovered
##  '.require.package' is the largest one, and it is never reached by a normal call: the
##  packages it guards are all installed when the suite runs. It is therefore called directly
##  on a package name that cannot exist, so that the whole decision tree is walked without
##  installing anything.
## ----------------------------------------------------------------------------------------

require.package <-
  function(...) {
    DEprot:::.require.package(...)
  }

## a name no repository will ever serve
absent.package <- "aPackageThatDoesNotExist.DEprot"



## ----------------------------------------------------------------------------------------
##  .require.package
## ----------------------------------------------------------------------------------------

test_that("an available package is confirmed without touching the repositories", {
  expect_true(require.package(package = "stats", install.missing = "never"))
  expect_invisible(require.package(package = "stats", install.missing = "never"))
  ## the repo is irrelevant when the package is already there
  expect_true(require.package(package = "stats", repo = "bioc", install.missing = "never"))
})


test_that("a missing package stops with the command needed to install it", {
  ## 'never' must never reach the installer, whatever the repository
  expect_error(require.package(package = absent.package,
                               repo = "CRAN",
                               install.missing = "never"),
               "install.packages")

  expect_error(require.package(package = absent.package,
                               repo = "bioc",
                               install.missing = "never"),
               "BiocManager::install")

  expect_error(require.package(package = absent.package,
                               repo = "bioconductor",
                               install.missing = "never"),
               "BiocManager::install")

  ## anything else is read as a GitHub slug
  expect_error(require.package(package = absent.package,
                               repo = "user/repo",
                               install.missing = "never"),
               "remotes::install_github")
})


test_that("the reason for the requirement is reported when provided", {
  expect_error(require.package(package = absent.package,
                               install.missing = "never",
                               reason = "read the parquet reports"),
               "read the parquet reports")

  ## without a reason the message stays generic but still names the package
  expect_error(require.package(package = absent.package, install.missing = "never"),
               absent.package, fixed = TRUE)
})


test_that("a non-interactive session never installs, whatever the policy", {
  skip_if(interactive(), "this expectation describes the non-interactive branch")

  ## 'ask' would need a prompt and 'always' an installation: both are refused outside an
  ## interactive session, which is what makes the function safe under R CMD check
  expect_error(require.package(package = absent.package, install.missing = "ask"))
  expect_error(require.package(package = absent.package, install.missing = "always"))
})


test_that("an unknown installation policy is rejected", {
  expect_error(require.package(package = absent.package, install.missing = "maybe"))
  expect_error(require.package(package = absent.package, install.missing = TRUE))
})



## ----------------------------------------------------------------------------------------
##  .read.table.any
## ----------------------------------------------------------------------------------------

test_that("the parquet branch reports the reader it needs", {
  dir <- local.tmpdir()
  path <- file.path(dir, "report.parquet")
  writeLines("not a real parquet file", path)

  skip_if(requireNamespace("nanoparquet", quietly = TRUE) | requireNamespace("arrow", quietly = TRUE),
          "a parquet reader is installed: the guard is not reached")

  ## without a reader the function must name the package to install instead of failing on
  ## the binary content
  expect_error(DEprot:::.read.table.any(path, install.missing = "never"), "nanoparquet")
})


test_that("the extension is matched whatever its case", {
  dir <- local.tmpdir()
  tb <- data.frame(id = c("P1", "P2"), s1 = c(1, 2), stringsAsFactors = FALSE)

  upper <- file.path(dir, "TABLE.TSV")
  utils::write.table(tb, upper, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_s3_class(DEprot:::.read.table.any(upper), "data.frame")
})


test_that("the column names are preserved exactly as they are in the file", {
  dir <- local.tmpdir()
  path <- file.path(dir, "table.tsv")

  tb <- data.frame(`Protein Group` = "P1",
                   `LFQ intensity A` = 1,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
  utils::write.table(tb, path, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- DEprot:::.read.table.any(path)

  ## 'check.names' must stay off: the sample names live in the headers
  expect_equal(colnames(out), c("Protein Group", "LFQ intensity A"))
})



## ----------------------------------------------------------------------------------------
##  .clean.run.names
## ----------------------------------------------------------------------------------------

test_that("every extension of the acquisition files is stripped", {
  clean <- DEprot:::.clean.run.names

  for (ext in c("raw", "mzML", "mzXML", "d", "wiff", "dia", "htrms", "timsTOF")) {
    expect_equal(clean(paste0("sample.", ext)), "sample")
    ## the match is case-insensitive
    expect_equal(clean(paste0("sample.", toupper(ext))), "sample")
  }
})


test_that("the windows separators are handled", {
  clean <- DEprot:::.clean.run.names

  expect_equal(clean("D:\\\\data\\\\run\\\\sampleA.raw"), "sampleA")
  expect_equal(clean("D:/data/sampleB.raw"), "sampleB")
})


test_that("the surrounding blanks are removed and the vector length is kept", {
  clean <- DEprot:::.clean.run.names

  expect_equal(clean("  sample  "), "sample")
  expect_length(clean(c("a.raw", "b.raw", "c.raw")), 3)
  expect_equal(clean(character(0)), character(0))
})


test_that("a trailing blank currently prevents the extension from being stripped", {
  clean <- DEprot:::.clean.run.names

  ## NOTE: trimws() is applied AFTER the extension is removed, hence the '$' anchor of the
  ## extension pattern does not match when the name ends with a blank. Trimming first would
  ## fix it; the expectation below documents the behaviour as it is today.
  expect_equal(clean("sample.raw  "), "sample.raw")
  expect_equal(clean("sample.raw"), "sample")
})


test_that("a name carrying no extension is left as it is", {
  clean <- DEprot:::.clean.run.names

  ## an internal dot is not an extension and must survive
  expect_equal(clean("sample.01.replicate"), "sample.01.replicate")
})



## ----------------------------------------------------------------------------------------
##  .long.to.matrix and .wide.to.matrix
## ----------------------------------------------------------------------------------------

test_that("the long reshaping drops the rows without an identifier", {
  long <- data.frame(id = c("P1", NA, "", "P2"),
                     sample = c("s1", "s1", "s1", "s1"),
                     value = c(10, 20, 25, 30),
                     stringsAsFactors = FALSE)

  ## the filtering is silent here: only the wide reader reports it, since a missing ID in a
  ## wide report means a whole protein is lost rather than a single measurement
  out <- DEprot:::.long.to.matrix(df = long, id.col = "id",
                                  sample.col = "sample", quantity.col = "value")

  expect_equal(sort(rownames(out)), c("P1", "P2"))
})


test_that("the wide reshaping reports the rows without an identifier", {
  wide <- data.frame(id = c("P1", NA, "P2"),
                     A = c(1, 2, 3),
                     B = c(4, 5, 6),
                     stringsAsFactors = FALSE)

  expect_warning(DEprot:::.wide.to.matrix(wide, "id", c(sampleA = "A", sampleB = "B")))

  out <- suppressWarnings(DEprot:::.wide.to.matrix(wide, "id", c(sampleA = "A", sampleB = "B")))
  expect_equal(sort(rownames(out)), c("P1", "P2"))
})


test_that("an empty wide table is reported", {
  empty <- data.frame(id = character(0), A = numeric(0), stringsAsFactors = FALSE)

  expect_error(DEprot:::.wide.to.matrix(empty, "id", c(sampleA = "A")))
})


test_that("the duplicated identifiers of a wide report are made unique", {
  wide <- data.frame(id = c("P1", "P1", "P2"),
                     A = c(1, 2, 3),
                     stringsAsFactors = FALSE)

  out <- DEprot:::.wide.to.matrix(wide, "id", c(sampleA = "A"))

  expect_equal(rownames(out), c("P1", "P1.1", "P2"))
})


test_that("a table left empty by the filters is reported", {
  ## every value is zero: with 'zero.to.na' nothing survives the conversion
  long <- data.frame(id = c("P1", "P2"),
                     sample = c("s1", "s2"),
                     value = c(0, 0),
                     stringsAsFactors = FALSE)

  expect_error(DEprot:::.long.to.matrix(df = long, id.col = "id",
                                        sample.col = "sample", quantity.col = "value"))
})


test_that("the aggregating function can be changed", {
  long <- data.frame(id = rep("P1", 3),
                     sample = rep("s1", 3),
                     value = c(10, 20, 60),
                     stringsAsFactors = FALSE)

  summed <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                     quantity.col = "value")
  median <- DEprot:::.long.to.matrix(df = long, id.col = "id", sample.col = "sample",
                                     quantity.col = "value",
                                     FUN = function(x) {stats::median(x, na.rm = TRUE)})

  expect_equal(unname(summed["P1", "s1"]), 90)
  expect_equal(unname(median["P1", "s1"]), 20)
})


test_that("the matrix is numeric even when the quantities are read as text", {
  long <- data.frame(id = c("P1", "P2"),
                     sample = c("s1", "s1"),
                     value = c("10", "20"),
                     stringsAsFactors = FALSE)

  out <- DEprot:::.long.to.matrix(df = long, id.col = "id",
                                  sample.col = "sample", quantity.col = "value")

  expect_true(is.numeric(out))
  expect_equal(unname(out["P2", "s1"]), 20)
})


test_that("the wide reshaping keeps the row names and converts the values", {
  wide <- data.frame(id = c("P1", "P2"),
                     A = c("1", "3"),
                     B = c("2", "4"),
                     stringsAsFactors = FALSE)

  out <- DEprot:::.wide.to.matrix(wide, "id", c(sampleA = "A", sampleB = "B"))

  expect_true(is.numeric(out))
  expect_equal(rownames(out), c("P1", "P2"))
  expect_equal(colnames(out), c("sampleA", "sampleB"))
})



## ----------------------------------------------------------------------------------------
##  .finalize.import
## ----------------------------------------------------------------------------------------

imported.list <-
  function(counts = matrix(1:4, nrow = 2, dimnames = list(c("P1", "P2"), c("a", "b"))),
           ...) {
    utils::modifyList(list(counts = counts,
                           log.base = 1,
                           data.type = "raw",
                           normalization.method = NA,
                           quantity = "test",
                           feature = "proteins"),
                      list(...))
  }


finalize <-
  function(imported = imported.list(), ...) {
    DEprot:::.finalize.import(imported = imported, verbose = FALSE, ...)
  }


test_that("a metadata listing samples absent from the counts is reported", {
  meta <- data.frame(column.id = c("a", "b", "c"), stringsAsFactors = FALSE)

  expect_error(finalize(metadata = meta))
})


test_that("samples absent from the metadata are either dropped or refused", {
  meta <- data.frame(column.id = "a", stringsAsFactors = FALSE)

  ## without 'subset.to.metadata' the mismatch is an error
  expect_error(finalize(metadata = meta, subset.to.metadata = FALSE))

  ## with it, the extra samples are silently removed
  out <- suppressMessages(finalize(metadata = meta, subset.to.metadata = TRUE))
  expect_equal(ncol(any.counts(out)), 1)
  expect_equal(nrow(out@metadata), 1)
})


test_that("the metadata column matching the samples can be renamed", {
  meta <- data.frame(sample.name = c("a", "b"), condition = c("x", "y"), stringsAsFactors = FALSE)

  out <- suppressMessages(finalize(metadata = meta, column.id = "sample.name"))

  expect_equal(nrow(out@metadata), 2)
  expect_error(finalize(metadata = meta, column.id = "not.a.column"))
})


test_that("the log base proposed by the reader can be overridden", {
  linear <- suppressWarnings(suppressMessages(finalize()))
  log2 <- suppressWarnings(suppressMessages(finalize(log.base = 2)))

  ## a base of 1 declares linear intensities: 'load.counts2' log2-transforms them and the
  ## object therefore reports base 2, on the scale the counts are actually stored in
  expect_equal(linear@log.base, 2)
  expect_true(linear@log.transformed)

  ## with the counts already in log2 nothing is transformed
  expect_equal(log2@log.base, 2)
  expect_equal(as.numeric(any.counts(log2)), as.numeric(matrix(1:4, nrow = 2)))
})


test_that("a normalized or imputed import fills the matching slot", {
  normalized <- suppressWarnings(suppressMessages(
    finalize(imported.list(data.type = "normalized", normalization.method = "test"))))

  imputed <- suppressWarnings(suppressMessages(
    finalize(imported.list(data.type = "imputed", imputation.method = "test"))))

  expect_true(normalized@normalized)
  expect_equal(normalized@normalization.method, "test")

  expect_true(imputed@imputed)
  expect_false(is.null(imputed@imputed.counts))
})


test_that("an imputed import without a recorded method is still flagged", {
  out <- suppressWarnings(suppressMessages(finalize(imported.list(data.type = "imputed"))))

  expect_true(out@imputed)
  expect_true(is.na(out@imputation.method) | is.character(out@imputation.method))
})


test_that("the messages describe what was imported", {
  expect_message(DEprot:::.finalize.import(imported = imported.list(),
                                           metadata = data.frame(column.id = c("a", "b")),
                                           verbose = TRUE),
                 "proteins")
})
