## ----------------------------------------------------------------------------------------
##  generate.mm()
##  The paragraphs are read back from the object they describe, hence the assertions compare
##  the text against the slots rather than against a hard-coded string: the wording can be
##  reworded without breaking the suite, the numbers cannot.
## ----------------------------------------------------------------------------------------

## Drafts reused across the tests
mm.dpa <- generate.mm(tb.limma, verbose = FALSE)
mm.dpo <- generate.mm(tb.dpo.imp, verbose = FALSE)
mm.raw <- generate.mm(tb.dpo.raw, verbose = FALSE)


#' Citation numbers actually used in the text, ranges expanded
#'
#' Only the brackets holding digits, commas and dashes are citations: every other bracket of
#' the draft carries at least one letter (sample sizes, seeds, versions, fitting methods).

cited.numbers <-
  function(text) {
    brackets <- unlist(regmatches(text, gregexpr("\\([0-9,-]+\\)", text)))
    items <- unlist(strsplit(gsub("[()]", "", brackets), ","))

    numbers <- lapply(items,
                      function(i) {
                        if (grepl("-", i)) {
                          range <- as.numeric(unlist(strsplit(i, "-")))
                          return(range[1]:range[2])
                        }
                        return(as.numeric(i))
                      })

    return(sort(unique(unlist(numbers))))
  }



## ----------------------------------------------------------------------------------------
##  Output structure
## ----------------------------------------------------------------------------------------

test_that("generate.mm returns the four elements of the output list", {
  expect_type(mm.dpa, "list")
  expect_named(mm.dpa, c("text", "references", "parameters", "full.text"))
  expect_type(mm.dpa$text, "character")
  expect_length(mm.dpa$text, 1)
  expect_s3_class(mm.dpa$references, "data.frame")
  expect_s3_class(mm.dpa$parameters, "data.frame")
})


test_that("the output is returned invisibly and printed only when verbose", {
  expect_invisible(generate.mm(tb.limma, verbose = FALSE))
  expect_output(generate.mm(tb.limma, verbose = TRUE), "MATERIAL AND METHODS")
  expect_silent(invisible(generate.mm(tb.limma, verbose = FALSE)))
})


test_that("generate.mm rejects objects of the wrong class", {
  expect_error(generate.mm(DEprot::test.toolbox$geneset))
  expect_error(generate.mm(DEprot::test.toolbox$unimputed.lfq))
  expect_error(generate.mm(tb.limma, enrichment.object = tb.dpo.imp, verbose = FALSE))
})


test_that("the parameters table collects the settings of each step", {
  expect_named(mm.dpa$parameters, c("step", "param", "value"))
  expect_type(mm.dpa$parameters$value, "character")
  expect_true(all(nchar(mm.dpa$parameters$value) > 0))

  expect_true(all(c("dataset", "differential analyses") %in% mm.dpa$parameters$step))
  expect_true(all(c("normalization", "imputation") %in% mm.dpo$parameters$step))

  ## the values are read from the object, not recomputed
  test <- mm.dpa$parameters$value[mm.dpa$parameters$step == "differential analyses" &
                                    mm.dpa$parameters$param == "test"]
  expect_equal(test, tb.limma@differential.analyses.params$stat.test)
})



## ----------------------------------------------------------------------------------------
##  Citations
## ----------------------------------------------------------------------------------------

test_that("no citation placeholder survives in the final text", {
  expect_false(grepl("[[", mm.dpa$text, fixed = TRUE))
  expect_false(grepl("]]", mm.dpa$text, fixed = TRUE))
  expect_false(any(grepl("missing reference", mm.dpa$references$citation, fixed = TRUE)))
  expect_false(any(grepl("missing reference", mm.dpo$references$citation, fixed = TRUE)))
  expect_false(any(grepl("missing reference", mm.raw$references$citation, fixed = TRUE)))
})


test_that("the reference list is numbered sequentially and matches the text", {
  expect_equal(mm.dpa$references$n, seq_len(nrow(mm.dpa$references)))
  expect_setequal(cited.numbers(mm.dpa$text), mm.dpa$references$n)
  expect_setequal(cited.numbers(mm.dpo$text), mm.dpo$references$n)
})


test_that("the numbering follows the order of first appearance", {
  ## R is the first tool mentioned, DEprot and its Zenodo record come right after
  expect_match(mm.dpa$references$citation[1], "R Core Team")
  expect_match(mm.dpa$references$citation[2], "nargab")
  expect_match(mm.dpa$references$citation[3], "zenodo", ignore.case = TRUE)
})


test_that("the DEprot paper and the Zenodo record are always cited", {
  for (draft in list(mm.dpa, mm.dpo, mm.raw)) {
    expect_true(any(grepl("10.1093/nargab/lqag015", draft$references$citation, fixed = TRUE)))
    expect_true(any(grepl("10.5281/zenodo.18233890", draft$references$citation, fixed = TRUE)))
  }
})


test_that("consecutive citations are merged in a single bracket", {
  expect_match(mm.dpa$text, "\\([0-9]+,[0-9]+\\)")
})


test_that("the tools of each step are cited", {
  ## normalization, randomization and imputation
  expect_true(any(grepl("Brombacher", mm.dpo$references$citation)))
  expect_true(any(grepl("Stekhoven", mm.dpo$references$citation)))

  ## differential analyses and multiple testing correction
  expect_true(any(grepl("Ritchie", mm.dpa$references$citation)))
  expect_true(any(grepl("Benjamini", mm.dpa$references$citation)))

  ## a plain object carries no differential analyses, hence no test to cite
  expect_false(any(grepl("Ritchie", mm.raw$references$citation)))
})


test_that("the version of the installed packages can be appended to the references", {
  with.version <- generate.mm(tb.limma, include.package.versions = TRUE, verbose = FALSE)
  no.version <- generate.mm(tb.limma, include.package.versions = FALSE, verbose = FALSE)

  expect_true(any(grepl("[R-package v", with.version$references$citation, fixed = TRUE)))
  expect_false(any(grepl("[R-package v", no.version$references$citation, fixed = TRUE)))

  ## the version reported is the one of the session
  expect_true(any(grepl(paste0("[R-package v", utils::packageVersion("DEprot"), "]"),
                        with.version$references$citation, fixed = TRUE)))
})


test_that("the quantification software is mentioned and cited when a reference is given", {
  with.ref <- generate.mm(tb.limma,
                          quantification.software = "MaxQuant (v2.4.2)",
                          quantification.reference = "Cox J., Mann M. MaxQuant. Nat Biotechnol 26, 2008.",
                          verbose = FALSE)

  no.ref <- generate.mm(tb.limma,
                        quantification.software = "DIA-NN (v1.8)",
                        verbose = FALSE)

  ## the citation provided by the user opens the list
  expect_match(with.ref$references$citation[1], "MaxQuant")
  expect_match(with.ref$text, "MaxQuant (v2.4.2) (1)", fixed = TRUE)

  ## without a citation the software is named but nothing is added to the list
  expect_match(no.ref$text, "DIA-NN (v1.8)", fixed = TRUE)
  expect_false(any(grepl("DIA-NN", no.ref$references$citation, fixed = TRUE)))
  expect_equal(nrow(no.ref$references), nrow(mm.dpa$references))
})


test_that("extra references are appended at the end and left uncited", {
  extra <- c("Extra reference A, 2001.", "Extra reference B, 2002.")
  draft <- generate.mm(tb.limma, extra.references = extra, verbose = FALSE)

  expect_equal(nrow(draft$references), nrow(mm.dpa$references) + length(extra))
  expect_equal(utils::tail(draft$references$citation, 2), extra)
  expect_false(any(utils::tail(draft$references$n, 2) %in% cited.numbers(draft$text)))
})



## ----------------------------------------------------------------------------------------
##  Content read from the object
## ----------------------------------------------------------------------------------------

test_that("the size of the dataset is reported", {
  counts <- tb.dpo.imp@imputed.counts

  expect_match(mm.dpo$text, paste(nrow(counts), "proteins"), fixed = TRUE)
  expect_match(mm.dpo$text, paste(ncol(counts), "samples"), fixed = TRUE)
  expect_match(mm.dpo$text, "log2-transformed", fixed = TRUE)
})


test_that("the normalization parameters are reported", {
  normalization <- tb.dpo.imp@normalization.method

  expect_match(mm.dpo$text, "MBQN", fixed = TRUE)
  expect_match(mm.dpo$text, normalization$value[normalization$param == "function"], fixed = TRUE)
  expect_match(mm.dpo$text,
               normalization$value[normalization$param == "NRI/RI ratio threshold"],
               fixed = TRUE)
})


test_that("the randomization parameters are reported", {
  randomization <- tb.dpo.imp@randomization.method

  expect_match(mm.dpo$text, randomization$group.column, fixed = TRUE)
  expect_match(mm.dpo$text, paste0(randomization$percentage.missing, "%"), fixed = TRUE)
  expect_match(mm.dpo$text, paste0("bottom ", randomization$tail.percentage, "%"), fixed = TRUE)
})


test_that("the imputation parameters are reported and the seeds are not rounded", {
  imputation <- tb.dpo.imp@imputation.method

  expect_match(mm.dpo$text, "missForest", fixed = TRUE)
  expect_match(mm.dpo$text, paste(imputation$max.iterations, "iterations"), fixed = TRUE)
  expect_match(mm.dpo$text, "out-of-bag", fixed = TRUE)

  ## a rounded seed would make the analyses impossible to reproduce
  expect_match(mm.dpo$text, format(imputation$seed, scientific = FALSE), fixed = TRUE)
  expect_match(mm.dpo$text,
               format(tb.dpo.imp@randomization.method$seed, scientific = FALSE),
               fixed = TRUE)
})


test_that("the differential analyses parameters and the contrasts are reported", {
  parameters <- tb.limma@differential.analyses.params

  expect_match(mm.dpa$text, "limma", fixed = TRUE)
  expect_match(mm.dpa$text, parameters$counts.used, fixed = TRUE)
  expect_match(mm.dpa$text, as.character(parameters$padj.th), fixed = TRUE)
  expect_match(mm.dpa$text, paste0("linear fold change of ", parameters$linear.FC.th), fixed = TRUE)

  ## the paired design is explicitly declared
  expect_match(mm.dpa$text, parameters$replicate.column, fixed = TRUE)

  ## every contrast of the object is described, with the size of the two groups
  for (contrast in tb.limma@contrasts) {
    expect_match(mm.dpa$text,
                 paste0(contrast$var.1, " versus ", contrast$var.2, " (n = ",
                        length(contrast$group.1), " and n = ", length(contrast$group.2)),
                 fixed = TRUE)
  }
})


test_that("the number of differential proteins matches the n.diff table", {
  contrast <- tb.limma@contrasts[[1]]
  n.diff <- tb.limma@analyses.result.list[[1]]$n.diff
  n.up <- n.diff$n[as.character(n.diff$diff.status) == contrast$var.1]

  expect_match(mm.dpa$text, "These thresholds returned", fixed = TRUE)
  expect_match(mm.dpa$text, paste0(n.up, " enriched in ", contrast$var.1), fixed = TRUE)

  no.summary <- generate.mm(tb.limma, include.results.summary = FALSE, verbose = FALSE)
  expect_false(grepl("These thresholds returned", no.summary$text, fixed = TRUE))
})


test_that("the contrasts to describe can be subset", {
  one <- generate.mm(tb.limma, contrasts = 1, verbose = FALSE)

  expect_match(one$text,
               paste0(tb.limma@contrasts[[1]]$var.1, " versus ", tb.limma@contrasts[[1]]$var.2),
               fixed = TRUE)
  expect_false(grepl(paste0(tb.limma@contrasts[[2]]$var.1, " versus ",
                            tb.limma@contrasts[[2]]$var.2),
                     one$text, fixed = TRUE))

  expect_error(generate.mm(tb.limma, contrasts = 100, verbose = FALSE))
  expect_error(generate.mm(tb.limma, contrasts = "1", verbose = FALSE))
})


test_that("a DEprot object gets no differential analyses paragraph", {
  expect_false(grepl("Differential protein abundance", mm.dpo$text, fixed = TRUE))
  expect_false(grepl("Differential protein abundance", mm.raw$text, fixed = TRUE))
})


test_that("the steps that were not applied are declared as such", {
  expect_match(mm.raw$text, "No normalization was applied", fixed = TRUE)
  expect_match(mm.raw$text, "No imputation was performed", fixed = TRUE)
  expect_false(grepl("MBQN", mm.raw$text, fixed = TRUE))
  expect_false(grepl("missForest", mm.raw$text, fixed = TRUE))
})


test_that("a synthetic object without recorded settings does not leak NA in the text", {
  dpo <- make.dpo()
  draft <- generate.mm(dpo, verbose = FALSE)

  expect_false(grepl("NA", draft$text, fixed = TRUE))
  expect_match(draft$text, paste(nrow(any.counts(dpo)), "proteins"), fixed = TRUE)
})



## ----------------------------------------------------------------------------------------
##  Geneset enrichment
## ----------------------------------------------------------------------------------------

test_that("no enrichment paragraph is written without an enrichment object", {
  expect_false(grepl("clusterProfiler", mm.dpa$text, fixed = TRUE))
  expect_false(any(grepl("clusterProfiler", mm.dpa$references$citation, fixed = TRUE)))
})


test_that("GSEA and ORA are described differently", {
  gsea <- generate.mm(tb.limma,
                      enrichment.object = DEprot::test.toolbox$gsea.results,
                      geneset.database = "CORUM v5.0",
                      verbose = FALSE)

  ora <- generate.mm(tb.limma,
                     enrichment.object = DEprot::test.toolbox$ora.results,
                     geneset.database = "CORUM v5.0",
                     verbose = FALSE)

  expect_match(gsea$text, "Gene set enrichment analysis", fixed = TRUE)
  expect_true(any(grepl("Subramanian", gsea$references$citation)))

  expect_match(ora$text, "Over-representation analysis", fixed = TRUE)
  expect_false(grepl("Gene set enrichment analysis", ora$text, fixed = TRUE))

  ## the collection and clusterProfiler are reported in both cases
  for (draft in list(gsea, ora)) {
    expect_match(draft$text, "CORUM v5.0", fixed = TRUE)
    expect_true(any(grepl("clusterProfiler", draft$references$citation)))
    expect_true("enrichment" %in% draft$parameters$step)
  }
})


test_that("the geneset collection can be cited", {
  citation <- "Tsitsiridis G. et al. CORUM. Nucleic Acids Research 51, D539-D545, 2023."

  draft <- generate.mm(tb.limma,
                       enrichment.object = DEprot::test.toolbox$ora.results,
                       geneset.database = "CORUM v5.0",
                       geneset.reference = citation,
                       verbose = FALSE)

  expect_true(any(grepl("Tsitsiridis", draft$references$citation)))
  expect_setequal(cited.numbers(draft$text), draft$references$n)
})



## ----------------------------------------------------------------------------------------
##  Formatting and export
## ----------------------------------------------------------------------------------------

test_that("the headers can be switched off", {
  expect_match(mm.dpa$full.text, "MATERIAL AND METHODS", fixed = TRUE)
  expect_match(mm.dpa$full.text, "REFERENCES", fixed = TRUE)

  no.headers <- generate.mm(tb.limma, add.headers = FALSE, verbose = FALSE)
  expect_false(grepl("MATERIAL AND METHODS", no.headers$full.text, fixed = TRUE))
  expect_false(grepl("REFERENCES", no.headers$full.text, fixed = TRUE))
})


test_that("the text is wrapped at the requested width", {
  wrapped <- generate.mm(tb.limma, wrap.width = 60, verbose = FALSE)
  unwrapped <- generate.mm(tb.limma, wrap.width = NULL, add.headers = FALSE, verbose = FALSE)

  expect_true(max(nchar(unlist(strsplit(wrapped$full.text, "\n")))) <= 60)
  expect_true(max(nchar(unlist(strsplit(unwrapped$full.text, "\n")))) > 60)

  ## the paragraphs stay separated by an empty line whatever the wrapping
  expect_true(grepl("\n\n", wrapped$full.text, fixed = TRUE))
})


test_that("the draft can be written to a file", {
  dir <- local.tmpdir()
  out <- file.path(dir, "material.and.methods.txt")

  draft <- generate.mm(tb.limma, output.file = out, verbose = FALSE)

  expect_true(file.exists(out))
  expect_equal(paste(readLines(out), collapse = "\n"), draft$full.text)
  expect_true(any(grepl("MATERIAL AND METHODS", readLines(out), fixed = TRUE)))
})
