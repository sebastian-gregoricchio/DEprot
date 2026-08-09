## ----------------------------------------------------------------------------------------
##  export.analyses(), export.external(), export.report()
##  Everything is written to a temporary folder removed at the end of each test.
## ----------------------------------------------------------------------------------------

test_that("export.analyses writes the expected folder structure", {
  dir <- local.tmpdir()
  out <- file.path(dir, "export")

  invisible(export.analyses(DEprot.analyses.object = tb.limma,
                            output.folder = out,
                            verbose = FALSE))

  expect_true(dir.exists(out))
  expect_true(dir.exists(file.path(out, "counts")))
  expect_true(dir.exists(file.path(out, "differential_analyses")))
  expect_true(length(list.files(out, recursive = TRUE)) > 0)
})


test_that("export.analyses writes the metadata and the counts tables", {
  dir <- local.tmpdir()
  out <- file.path(dir, "export")

  invisible(export.analyses(DEprot.analyses.object = tb.limma,
                            output.folder = out,
                            verbose = FALSE))

  files <- list.files(out, recursive = TRUE)

  expect_true(any(grepl("metadata", files)))
  expect_true(any(grepl("counts", files)))
})


test_that("export.analyses can be restricted to a subset of contrasts", {
  dir <- local.tmpdir()
  out <- file.path(dir, "export")

  invisible(export.analyses(DEprot.analyses.object = tb.limma,
                            output.folder = out,
                            contrast.subset = 1,
                            verbose = FALSE))

  contrast.dirs <- list.dirs(file.path(out, "differential_analyses"), recursive = FALSE)
  expect_equal(length(contrast.dirs), 1)
})


test_that("export.analyses rejects an object of the wrong class", {
  dir <- local.tmpdir()

  expect_error(export.analyses(DEprot.analyses.object = tb.dpo.norm,
                               output.folder = file.path(dir, "export"),
                               verbose = FALSE))
})



## ----------------------------------------------------------------------------------------
##  export.external()
##  'install.missing' is a string: "never" keeps the tests offline and makes them fail
##  loudly rather than trying to install a package during a check.
## ----------------------------------------------------------------------------------------

test_that("a SummarizedExperiment is built", {
  skip_if_not_installed("SummarizedExperiment")

  se <- export.external(DEprot.object = tb.limma,
                        format = "SummarizedExperiment",
                        install.missing = "never",
                        verbose = FALSE)

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), nrow(tb.limma@metadata))
})


test_that("the differential results are added to the row annotation", {
  skip_if_not_installed("SummarizedExperiment")

  se <- export.external(DEprot.object = tb.limma,
                        format = "SummarizedExperiment",
                        add.results = TRUE,
                        install.missing = "never",
                        verbose = FALSE)

  expect_true(ncol(SummarizedExperiment::rowData(se)) > 0)

  se.bare <- export.external(DEprot.object = tb.limma,
                             format = "SummarizedExperiment",
                             add.results = FALSE,
                             add.protein.info = FALSE,
                             install.missing = "never",
                             verbose = FALSE)

  expect_true(ncol(SummarizedExperiment::rowData(se.bare)) <=
                ncol(SummarizedExperiment::rowData(se)))
})


test_that("a QFeatures object is built", {
  skip_if_not_installed("QFeatures")

  qf <- export.external(DEprot.object = tb.limma,
                        format = "QFeatures",
                        install.missing = "never",
                        verbose = FALSE)

  expect_s4_class(qf, "QFeatures")
})


test_that("an MSnSet object is built", {
  skip_if_not_installed("MSnbase")

  ms <- export.external(DEprot.object = tb.limma,
                        format = "MSnSet",
                        install.missing = "never",
                        verbose = FALSE)

  expect_s4_class(ms, "MSnSet")
})


test_that("an EList object is built", {
  skip_if_not_installed("limma")

  el <- export.external(DEprot.object = tb.limma,
                        format = "EList",
                        install.missing = "never",
                        verbose = FALSE)

  expect_s4_class(el, "EList")
})


test_that("the as.SummarizedExperiment alias dispatches correctly", {
  skip_if_not_installed("SummarizedExperiment")

  se <- as.SummarizedExperiment(DEprot.object = tb.limma, install.missing = "never", verbose = FALSE)

  expect_s4_class(se, "SummarizedExperiment")
})


test_that("a single assay can be exported", {
  skip_if_not_installed("SummarizedExperiment")

  se <- export.external(DEprot.object = tb.limma,
                        assays = "normalized",
                        install.missing = "never",
                        verbose = FALSE)

  expect_equal(length(SummarizedExperiment::assays(se)), 1)
})


test_that("an unknown format, assay or install policy is rejected", {
  expect_error(export.external(DEprot.object = tb.limma, format = "not.a.format", verbose = FALSE))
  expect_error(export.external(DEprot.object = tb.limma, assays = "not.an.assay", verbose = FALSE))
  expect_error(export.external(DEprot.object = tb.limma, install.missing = FALSE, verbose = FALSE))
  expect_error(export.external(DEprot.object = "not.an.object", verbose = FALSE))
})



## ----------------------------------------------------------------------------------------
##  export.report()
##  The rendering needs pandoc, which is not available on every checking machine.
## ----------------------------------------------------------------------------------------

test_that("an HTML report is rendered from a DEprot.analyses object", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")

  dir <- local.tmpdir()
  out <- file.path(dir, "report.html")

  invisible(suppressWarnings(export.report(DEprot.object = tb.limma,
                                           output.file = out,
                                           include.contrast.qc = FALSE,
                                           top.n.proteins = 5,
                                           quiet = TRUE)))

  expect_true(file.exists(out))
  expect_true(file.size(out) > 0)
})


test_that("the report can be built from a plain DEprot object and the Rmd can be kept", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")

  dir <- local.tmpdir()
  out <- file.path(dir, "report.qc.html")

  invisible(suppressWarnings(export.report(DEprot.object = tb.dpo.imp,
                                           output.file = out,
                                           report.title = "QC only",
                                           keep.Rmd = TRUE,
                                           quiet = TRUE)))

  expect_true(file.exists(out))
  expect_true(length(list.files(dir, pattern = "\\.Rmd$")) > 0)
})


test_that("export.report rejects an object of the wrong class", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")

  expect_error(export.report(DEprot.object = data.frame(a = 1),
                             output.file = file.path(tempdir(), "nope.html")))
})
