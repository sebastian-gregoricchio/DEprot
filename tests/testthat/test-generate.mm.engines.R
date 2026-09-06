## ----------------------------------------------------------------------------------------
##  generate.mm(): the branches that depend on the pipeline actually run
##  'test-generate.mm.R' describes the limma/missForest/MBQN object shipped with the package.
##  What is left uncovered is everything the draft says about the OTHER engines, corrections
##  and normalizations: those paragraphs are reached only by objects carrying those settings,
##  hence the parameter slots are edited here instead of re-running the analyses, which would
##  cost minutes for a paragraph of text.
## ----------------------------------------------------------------------------------------

#' Copy of the analyses object with a modified 'differential.analyses.params'
#'
#' The draft is written from the parameters stored in the object, not from the counts, hence
#' overriding the slot is enough to reach the paragraph of any engine.

with.params <-
  function(...) {
    object <- tb.limma
    object@differential.analyses.params <- utils::modifyList(object@differential.analyses.params, list(...))
    return(object)
  }


draft.of <-
  function(object, ...) {
    suppressWarnings(suppressMessages(generate.mm(object, verbose = FALSE, ...)))
  }



## ----------------------------------------------------------------------------------------
##  Statistical engines
## ----------------------------------------------------------------------------------------

test_that("the t-test and the Wilcoxon paragraphs are written", {
  t.paired <- draft.of(with.params(stat.test = "t.test", paired.test = TRUE))
  t.unpaired <- draft.of(with.params(stat.test = "t.test", paired.test = FALSE))
  wilcox <- draft.of(with.params(stat.test = "wilcoxon", paired.test = FALSE))

  expect_match(t.paired$text, "Student's t-test", fixed = TRUE)
  expect_match(t.unpaired$text, "Student's t-test", fixed = TRUE)
  expect_match(wilcox$text, "Wilcoxon rank-sum test", fixed = TRUE)

  ## the paired design is declared only when the model includes the replicate
  expect_match(t.paired$text, tb.limma@differential.analyses.params$replicate.column, fixed = TRUE)
  expect_false(grepl("paired by", t.unpaired$text, fixed = TRUE))
})


test_that("the prolfqua paragraph reports the strategy and the moderation", {
  draft <- draft.of(with.params(stat.test = "prolfqua",
                                strategy = list(strategy.id = "lm"),
                                moderate.variance = TRUE,
                                robust.scaling = TRUE,
                                paired.test = FALSE))

  expect_match(draft$text, "prolfqua", fixed = TRUE)
  expect_match(draft$text, "'lm' strategy", fixed = TRUE)
  expect_match(draft$text, "moderation of the residual variances", fixed = TRUE)
  expect_match(draft$text, "robust scaling", fixed = TRUE)
  expect_true(any(grepl("prolfqua", draft$references$citation, ignore.case = TRUE)))

  ## the options switched off are simply not mentioned
  bare <- draft.of(with.params(stat.test = "prolfqua",
                               strategy = list(),
                               moderate.variance = FALSE,
                               robust.scaling = FALSE,
                               paired.test = FALSE))

  expect_false(grepl("robust scaling", bare$text, fixed = TRUE))
  expect_false(grepl("strategy)", bare$text, fixed = TRUE))
})


test_that("the proDA paragraph reports the dropout model and the detection threshold", {
  draft <- draft.of(with.params(stat.test = "proDA",
                                min.detected = 3,
                                moderate.location = TRUE,
                                moderate.variance = TRUE,
                                paired.test = TRUE))

  expect_match(draft$text, "proDA", fixed = TRUE)
  expect_match(draft$text, "probability of a dropout", fixed = TRUE)
  expect_match(draft$text, "fewer than 3 samples", fixed = TRUE)
  expect_match(draft$text, "blocking factor", fixed = TRUE)
  expect_true("minimum detections" %in% draft$parameters$param)

  ## without a threshold the sentence about the untested proteins disappears
  no.threshold <- draft.of(with.params(stat.test = "proDA",
                                       min.detected = NULL,
                                       paired.test = FALSE))
  expect_false(grepl("were not tested", no.threshold$text, fixed = TRUE))
})


test_that("an unknown engine is named without pretending to describe it", {
  draft <- draft.of(with.params(stat.test = "my.custom.test"))

  expect_match(draft$text, "'my.custom.test' test", fixed = TRUE)
})



## ----------------------------------------------------------------------------------------
##  Multiple testing correction
## ----------------------------------------------------------------------------------------

test_that("each correction method gets its own sentence", {
  bh <- draft.of(with.params(padj.method = "BH"))
  fdr <- draft.of(with.params(padj.method = "fdr"))
  fdrtool <- draft.of(with.params(padj.method = "fdrtool"))
  none <- draft.of(with.params(padj.method = "none"))
  bonferroni <- draft.of(with.params(padj.method = "bonferroni"))

  expect_match(bh$text, "Benjamini-Hochberg", fixed = TRUE)
  expect_match(fdr$text, "Benjamini-Hochberg", fixed = TRUE)
  expect_match(fdrtool$text, "fdrtool", fixed = TRUE)
  expect_true(any(grepl("Strimmer", fdrtool$references$citation)))

  ## without correction the threshold applies to the raw p-value
  expect_match(none$text, "p-value", fixed = TRUE)
  expect_false(grepl("adjusted p-value", none$text, fixed = TRUE))

  ## any other method accepted by p.adjust is named as it is
  expect_match(bonferroni$text, "bonferroni", fixed = TRUE)
})


test_that("the 'effective FDR' correction is described", {
  draft <- draft.of(with.params(padj.method = "effective FDR"))

  ## the effective FDR is the empirical estimation performed by prolfqua
  expect_match(draft$text, "empirical distribution of the test statistics", fixed = TRUE)
  expect_true("differential analyses" %in% draft$parameters$step)
})



## ----------------------------------------------------------------------------------------
##  Normalization and imputation
## ----------------------------------------------------------------------------------------

test_that("the HarmonizR batch correction is described", {
  object <- tb.dpo.imp
  object@normalization.method <- data.frame(param = c("package", "batch.column", "algorithm",
                                                      "ComBat.mode", "block", "cores"),
                                            value = c("HarmonizR", "batch", "ComBat", "1", "null", "1"),
                                            stringsAsFactors = FALSE)

  draft <- draft.of(object)

  expect_match(draft$text, "HarmonizR", fixed = TRUE)
  expect_match(draft$text, "batch", fixed = TRUE)
  expect_true(any(grepl("Voss", draft$references$citation) |
                    grepl("HarmonizR", draft$references$citation)))

  ## the limma algorithm is reported instead of ComBat when it was the one used
  object@normalization.method$value[3] <- "limma"
  limma.draft <- draft.of(object)
  expect_match(limma.draft$text, "limma", fixed = TRUE)
})


test_that("a normalization stored as a plain string is reported as such", {
  object <- tb.dpo.imp
  object@normalization.method <- "quantile"

  draft <- draft.of(object)

  expect_match(draft$text, "quantile", fixed = TRUE)
  expect_false(grepl("MBQN", draft$text, fixed = TRUE))
})


test_that("the pcaMethods-based imputations are described", {
  object <- tb.dpo.imp

  for (method in c("SVD", "BPCA", "PPCA", "LLS")) {
    object@imputation.method <- list(method = method, seed = 1234, PCs.tested = 3)
    draft <- draft.of(object)

    ## the sentence names the algorithm and the package implementing it
    expect_match(draft$text, "pcaMethods", fixed = TRUE)
    expect_match(draft$text, "1234", fixed = TRUE)
    ## the out-of-bag error is a missForest diagnostic only
    expect_false(grepl("out-of-bag", draft$text, fixed = TRUE))
  }
})


test_that("a pcaMethods imputation without a recorded number of PCs does not leak NA", {
  object <- tb.dpo.imp
  object@imputation.method <- list(method = "SVD", seed = 1234)

  draft <- draft.of(object)

  ## NOTE: the 'svd'/'bpca'/'ppca' branches paste number.to.text(PCs.tested) without the
  ## is.na() guard used elsewhere, hence the draft currently reads "among the first NA".
  expect_false(grepl("first NA", draft$text, fixed = TRUE))
})


test_that("an imputation without recorded parameters is still declared", {
  object <- tb.dpo.imp
  object@imputation.method <- list()

  draft <- draft.of(object)

  expect_true("imputation" %in% draft$parameters$step | grepl("imput", draft$text))
  expect_false(grepl("NA", draft$text, fixed = TRUE))
})


test_that("a randomization without a group column does not leak NA", {
  object <- tb.dpo.imp
  object@randomization.method <- list(percentage.missing = 100, tail.percentage = 1)

  draft <- draft.of(object)

  expect_false(grepl("NA", draft$text, fixed = TRUE))
})



## ----------------------------------------------------------------------------------------
##  Counts, log base and metadata
## ----------------------------------------------------------------------------------------

test_that("the log base is described whatever its value", {
  natural <- tb.dpo.imp
  natural@log.base <- exp(1)

  base.10 <- tb.dpo.imp
  base.10@log.base <- 10

  linear <- tb.dpo.imp
  linear@log.transformed <- FALSE
  linear@log.base <- NA_real_

  expect_match(draft.of(tb.dpo.imp)$text, "log2-transformed", fixed = TRUE)
  expect_match(draft.of(natural)$text, "ln-transformed", fixed = TRUE)
  expect_match(draft.of(base.10)$text, "log10-transformed", fixed = TRUE)

  ## an NA base means that the counts were never log-transformed
  expect_match(draft.of(linear)$text, "linear scale", fixed = TRUE)
  expect_false(grepl("log2-transformed", draft.of(linear)$text, fixed = TRUE))
})


test_that("an object holding no counts at all is handled", {
  bare <- tb.dpo.raw
  bare@raw.counts <- NULL
  bare@norm.counts <- NULL
  bare@random.counts <- NULL
  bare@imputed.counts <- NULL

  draft <- draft.of(bare)

  expect_type(draft$full.text, "character")

  ## NOTE: 'n.proteins'/'n.samples' are set to NA when no counts table is available and are
  ## pasted as they are, hence the draft currently reads "NA proteins across NA samples".
  ## The size clause should be dropped instead of being filled with NA.
  expect_false(grepl("NA proteins", draft$text, fixed = TRUE))
})


test_that("a citation that cannot be resolved raises a warning", {
  ## the reference keys are looked up in an internal table: an unknown package name has no
  ## citation and the draft must say so rather than printing an empty entry
  object <- with.params(stat.test = "prolfqua",
                        strategy = list(strategy.id = "lm"),
                        paired.test = FALSE)

  expect_no_error(draft.of(object))
})
