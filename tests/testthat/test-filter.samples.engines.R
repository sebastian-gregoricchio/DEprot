## ----------------------------------------------------------------------------------------
##  filter.samples(): the engine dispatch and the parameter recovery
##  'test-filter.samples.R' covers the subsetting and the re-derivation of the counts. What
##  is left uncovered is the differential re-run: the function has to work out WHICH engine
##  produced the original analyses and to recover its arguments from the parameter slot,
##  which is a branch per engine and per legacy layout.
## ----------------------------------------------------------------------------------------

keep.most <- utils::head(tb.limma@metadata$column.id, -1)


#' Analyses object whose differential parameters have been overridden
#'
#' The engine is worked out from 'differential.analyses.params', hence editing the slot is
#' enough to reach the recovery code of any engine without re-running the analyses.

analyses.with <-
  function(...) {
    object <- tb.limma
    object@differential.analyses.params <- utils::modifyList(object@differential.analyses.params, list(...))
    return(object)
  }


filter.most <-
  function(object = tb.limma, ...) {
    suppressWarnings(suppressMessages(
      filter.samples(DEprot.object = object,
                     samples = keep.most,
                     mode = "keep",
                     verbose = FALSE,
                     ...)))
  }



## ----------------------------------------------------------------------------------------
##  Engine dispatch
## ----------------------------------------------------------------------------------------

test_that("the engine stored in the object is the one used by default", {
  out <- filter.most()

  expect_s4_class(out, "DEprot.analyses")
  expect_equal(tolower(out@differential.analyses.params$stat.test), "limma")
})


test_that("the engine can be forced whatever the one stored", {
  skip_on_cran()

  for (engine in c("limma", "t.test", "ttest", "wilcoxon")) {
    out <- filter.most(diff.method = engine)
    expect_s4_class(out, "DEprot.analyses")
  }
})


test_that("'wilcoxon' selects the t-test engine and passes the test down", {
  skip_on_cran()

  out <- filter.most(diff.method = "wilcoxon")

  ## the engine is 'diff.analyses', the test it runs is the Wilcoxon one
  expect_equal(tolower(out@differential.analyses.params$stat.test), "wilcoxon")
})


test_that("the test of the t-test engine can be chosen explicitly", {
  skip_on_cran()

  out <- filter.most(diff.method = "t.test", stat.test = "wilcoxon")

  expect_s4_class(out, "DEprot.analyses")
  expect_equal(tolower(out@differential.analyses.params$stat.test), "wilcoxon")
})


test_that("a legacy object without a recorded engine falls back on the t-test", {
  legacy <- analyses.with(stat.test = NULL)

  expect_warning(suppressMessages(
    filter.samples(DEprot.object = legacy, samples = keep.most, mode = "keep", verbose = FALSE)))

  out <- filter.most(object = legacy)
  expect_s4_class(out, "DEprot.analyses")
})


test_that("an unknown engine is rejected before anything is recomputed", {
  expect_error(filter.samples(DEprot.object = tb.limma,
                              samples = keep.most,
                              mode = "keep",
                              diff.method = "not.an.engine",
                              verbose = FALSE))
})



## ----------------------------------------------------------------------------------------
##  Parameter recovery
## ----------------------------------------------------------------------------------------

test_that("the thresholds of the original analyses are carried over", {
  out <- filter.most()
  original <- tb.limma@differential.analyses.params
  recomputed <- out@differential.analyses.params

  expect_equal(recomputed$linear.FC.th, original$linear.FC.th)
  expect_equal(recomputed$padj.th, original$padj.th)
  expect_equal(recomputed$padj.method, original$padj.method)
  expect_equal(recomputed$counts.used, original$counts.used)
})


test_that("the paired design is rebuilt from the stored replicate column", {
  out <- filter.most()

  expect_equal(out@differential.analyses.params$replicate.column,
               tb.limma@differential.analyses.params$replicate.column)
  expect_true(isTRUE(out@differential.analyses.params$rep.model) |
                isTRUE(out@differential.analyses.params$paired.test))
})


test_that("a paired analysis without a recoverable replicate column is refused", {
  ## the replicate column is the one piece the re-run cannot invent
  orphan <- analyses.with(replicate.column = NULL, rep.model = TRUE, paired.test = TRUE)
  orphan@metadata$replicate <- NULL

  expect_error(suppressMessages(
    filter.samples(DEprot.object = orphan, samples = keep.most, mode = "keep", verbose = FALSE)))
})


test_that("the replicate column can be provided when the object lost it", {
  skip_on_cran()

  orphan <- analyses.with(replicate.column = NULL, rep.model = TRUE)

  out <- filter.most(object = orphan, replicate.column = "replicate")

  expect_s4_class(out, "DEprot.analyses")
  expect_equal(out@differential.analyses.params$replicate.column, "replicate")
})


test_that("an unpaired analysis stays unpaired after the filtering", {
  skip_on_cran()

  unpaired <- analyses.with(rep.model = FALSE, paired.test = FALSE)

  out <- filter.most(object = unpaired)

  expect_false(isTRUE(out@differential.analyses.params$rep.model))
})


test_that("the counts used by the original analyses are used again", {
  skip_on_cran()

  on.normalized <- analyses.with(counts.used = "normalized")

  out <- filter.most(object = on.normalized)

  expect_equal(out@differential.analyses.params$counts.used, "normalized")
})



## ----------------------------------------------------------------------------------------
##  Contrast handling
## ----------------------------------------------------------------------------------------

test_that("every contrast is recomputed when the samples allow it", {
  out <- filter.most()

  expect_equal(length(out@analyses.result.list), length(tb.limma@analyses.result.list))
  expect_equal(names(out@contrasts), names(tb.limma@contrasts))
})


test_that("a contrast whose groups are too small is dropped with a warning", {
  meta <- tb.limma@metadata
  first.group <- as.character(meta$condition)[1]

  ## a single replicate is left in the first group: its contrasts cannot be tested
  keep <- c(meta$column.id[as.character(meta$condition) == first.group][1],
            meta$column.id[as.character(meta$condition) != first.group])

  expect_warning(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   min.samples.per.group = 2,
                   verbose = FALSE)))
})


test_that("an object left without any usable contrast returns a plain DEprot object", {
  meta <- tb.limma@metadata

  ## a single sample per condition: no contrast survives the threshold
  keep <- unlist(lapply(split(meta$column.id, as.character(meta$condition)),
                        function(x) {x[1]}))

  out <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   min.samples.per.group = 2,
                   verbose = FALSE)))

  expect_s4_class(out, "DEprot")
  expect_false("DEprot.analyses" %in% class(out))
})


test_that("'min.samples.per.group' drives how many contrasts survive", {
  permissive <- filter.most(min.samples.per.group = 2)
  strict <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep.most,
                   mode = "keep",
                   min.samples.per.group = 4,
                   verbose = FALSE)))

  expect_true(length(permissive@analyses.result.list) >=
                length(if (methods::is(strict, "DEprot.analyses")) {strict@analyses.result.list} else {list()}))
})



## ----------------------------------------------------------------------------------------
##  Messages and object bookkeeping
## ----------------------------------------------------------------------------------------

test_that("the steps of the rebuild are reported when verbose", {
  expect_message(suppressWarnings(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep.most,
                   mode = "keep",
                   verbose = TRUE)),
    "filter.samples")
})


test_that("the object returned is consistent with the samples retained", {
  out <- filter.most()

  expect_equal(sort(out@metadata$column.id), sort(keep.most))
  expect_equal(sort(colnames(any.counts(out))), sort(keep.most))

  ## every contrast lists only retained samples
  for (contrast in out@contrasts) {
    expect_true(all(c(contrast$group.1, contrast$group.2) %in% keep.most))
  }
})


test_that("the results tables cover the proteins of the filtered object", {
  out <- filter.most()

  results <- get.results(DEprot.analyses.object = out, contrast = 1)

  expect_true(nrow(results) > 0)
  expect_true(all(results$prot.id %in% rownames(any.counts(out))))
})
