## ----------------------------------------------------------------------------------------
##  filter.samples()
##  The function does not simply subset the counts: it re-derives every step recorded in the
##  object (normalization, randomization, imputation) on the retained samples and, for a
##  DEprot.analyses object, re-runs the differential analyses. Most of its body is therefore
##  reachable only through objects that actually carry those steps.
## ----------------------------------------------------------------------------------------

all.samples <- tb.dpo.imp@metadata$column.id
half.samples <- all.samples[1:(length(all.samples) / 2)]


######################################    ERRORS    ######################################

test_that("an object of the wrong class is rejected", {
  expect_error(filter.samples(DEprot.object = DEprot::sample.config,
                              samples = "a",
                              verbose = FALSE))
  expect_error(filter.samples(DEprot.object = "not an object",
                              samples = "a",
                              verbose = FALSE))
})


test_that("'mode' must be either 'remove' or 'keep'", {
  expect_error(filter.samples(DEprot.object = tb.dpo.imp,
                              samples = all.samples[1],
                              mode = "not.a.mode",
                              verbose = FALSE))
})


test_that("'samples' is mandatory and must not be empty", {
  expect_error(filter.samples(DEprot.object = tb.dpo.imp, verbose = FALSE))
  expect_error(filter.samples(DEprot.object = tb.dpo.imp, samples = NULL, verbose = FALSE))
  expect_error(filter.samples(DEprot.object = tb.dpo.imp, samples = character(0), verbose = FALSE))
})


test_that("'diff.method' is validated", {
  expect_error(filter.samples(DEprot.object = tb.dpo.imp,
                              samples = all.samples[1],
                              diff.method = "not.an.engine",
                              verbose = FALSE))
})


test_that("samples that are not in the metadata are reported", {
  expect_error(filter.samples(DEprot.object = tb.dpo.imp,
                              samples = c(all.samples[1], "not.a.sample"),
                              verbose = FALSE))
})


test_that("removing every sample is refused", {
  expect_error(filter.samples(DEprot.object = tb.dpo.imp,
                              samples = all.samples,
                              mode = "remove",
                              verbose = FALSE))
  expect_error(filter.samples(DEprot.object = tb.dpo.imp,
                              samples = character(0),
                              mode = "keep",
                              verbose = FALSE))
})


test_that("keeping a single sample raises a warning", {
  expect_warning(suppressMessages(filter.samples(DEprot.object = tb.dpo.raw,
                                                 samples = all.samples[1],
                                                 mode = "keep",
                                                 verbose = FALSE)))
})


test_that("an object without any counts table is rejected", {
  bare <- tb.dpo.raw
  bare@raw.counts <- NULL
  bare@norm.counts <- NULL
  bare@random.counts <- NULL
  bare@imputed.counts <- NULL

  expect_error(suppressWarnings(filter.samples(DEprot.object = bare,
                                               samples = all.samples[1],
                                               verbose = FALSE)))
})



##########################################################################################

test_that("'remove' and 'keep' are complementary", {
  removed <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.raw, samples = half.samples, mode = "remove", verbose = FALSE)))

  kept <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.raw, samples = half.samples, mode = "keep", verbose = FALSE)))

  expect_s4_class(removed, "DEprot")
  expect_s4_class(kept, "DEprot")
  expect_equal(nrow(removed@metadata), length(all.samples) - length(half.samples))
  expect_equal(nrow(kept@metadata), length(half.samples))
  expect_equal(sort(kept@metadata$column.id), sort(half.samples))
  expect_length(intersect(removed@metadata$column.id, kept@metadata$column.id), 0)
})


test_that("the counts tables follow the metadata", {
  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.raw, samples = half.samples, mode = "keep", verbose = FALSE)))

  expect_equal(sort(colnames(any.counts(filtered))), sort(half.samples))
})


test_that("the normalization is re-derived from the raw counts", {
  ## the object carries raw counts and a recoverable normalization method, hence the
  ## normalization can be recomputed on the retained samples instead of being subset
  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.norm, samples = half.samples, mode = "keep", verbose = FALSE)))

  expect_s4_class(filtered, "DEprot")
  expect_true(filtered@normalized)
  expect_equal(ncol(filtered@norm.counts), length(half.samples))
})


test_that("the randomization is re-derived when the object carries one", {
  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.random, samples = half.samples, mode = "keep", verbose = FALSE)))

  expect_s4_class(filtered, "DEprot")
  expect_true(filtered@randomized)
  expect_equal(ncol(filtered@random.counts), length(half.samples))
})


test_that("the imputation is re-derived when the object carries one", {
  skip_on_cran()

  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.dpo.imp, samples = half.samples, mode = "keep", verbose = FALSE)))

  expect_s4_class(filtered, "DEprot")
  expect_true(filtered@imputed)
  expect_equal(ncol(filtered@imputed.counts), length(half.samples))
})


test_that("an object flagged as randomized without usable parameters is rejected", {
  broken <- tb.dpo.random
  broken@randomization.method <- "not a list of parameters"

  expect_error(suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = broken, samples = half.samples, mode = "keep", verbose = FALSE))))
})


test_that("an object flagged as imputed without usable parameters is rejected", {
  broken <- tb.dpo.imp
  broken@imputation.method <- "not a list of parameters"

  expect_error(suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = broken, samples = half.samples, mode = "keep", verbose = FALSE))))
})


test_that("a batch-corrected object re-runs the harmonization", {
  skip_on_cran()

  dpo <- tb.dpo.raw
  ## the batches must be blocks of samples: an alternating design is too sparse for HarmonizR
  dpo@metadata$batch <- c(rep("A", nrow(dpo@metadata) / 2), rep("B", nrow(dpo@metadata) / 2))

  harm <- suppressWarnings(suppressMessages(
    harmonize.batches(DEprot.object = dpo, batch.column = "batch", cores = 1, verbose = FALSE)))

  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = harm,
                   samples = harm@metadata$column.id[1:8],
                   mode = "keep",
                   batch.column = "batch",
                   verbose = FALSE)))

  expect_s4_class(filtered, "DEprot")
  expect_equal(ncol(filtered@norm.counts), 8)
})


test_that("a batch column that does not exist stops the harmonization", {
  dpo <- tb.dpo.raw
  dpo@metadata$batch <- c(rep("A", nrow(dpo@metadata) / 2), rep("B", nrow(dpo@metadata) / 2))

  harm <- suppressWarnings(suppressMessages(
    harmonize.batches(DEprot.object = dpo, batch.column = "batch", cores = 1, verbose = FALSE)))

  expect_error(suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = harm,
                   samples = harm@metadata$column.id[1:8],
                   mode = "keep",
                   batch.column = "not.a.column",
                   verbose = FALSE))))
})



##################################    DEprot.analyses    #################################

test_that("the differential analyses are recomputed on the retained samples", {
  skip_on_cran()

  ## a single sample is dropped: removing more of them leaves proteins with no variance,
  ## which the PCA computed by the differential functions cannot scale
  keep <- utils::head(tb.limma@metadata$column.id, -1)

  filtered <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   diff.method = "limma",
                   verbose = FALSE)))

  expect_s4_class(filtered, "DEprot.analyses")
  expect_equal(nrow(filtered@metadata), length(keep))
  expect_true(length(filtered@analyses.result.list) >= 1)
})


test_that("the engine can be forced to the t-test/Wilcoxon implementation", {
  skip_on_cran()

  keep <- utils::head(tb.limma@metadata$column.id, -1)

  t.based <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   diff.method = "t.test",
                   verbose = FALSE)))

  wilcox <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   diff.method = "wilcoxon",
                   verbose = FALSE)))

  expect_s4_class(t.based, "DEprot.analyses")
  expect_s4_class(wilcox, "DEprot.analyses")
})


test_that("contrasts left with too few samples are skipped", {
  skip_on_cran()

  ## keeping a single replicate of the first condition makes its contrasts unusable
  meta <- tb.limma@metadata
  first.group <- meta$condition[1]
  keep <- c(meta$column.id[meta$condition == first.group][1],
            meta$column.id[meta$condition != first.group])

  expect_warning(suppressMessages(
    filter.samples(DEprot.object = tb.limma,
                   samples = keep,
                   mode = "keep",
                   min.samples.per.group = 3,
                   verbose = FALSE)))
})


test_that("an analyses object without contrasts returns a plain DEprot object", {
  no.contrast <- tb.limma
  no.contrast@contrasts <- list()
  no.contrast@analyses.result.list <- list()

  out <- suppressWarnings(suppressMessages(
    filter.samples(DEprot.object = no.contrast,
                   samples = utils::head(tb.limma@metadata$column.id, -1),
                   mode = "keep",
                   verbose = FALSE)))

  expect_s4_class(out, "DEprot")
})


test_that("the messages can be silenced and printed", {
  expect_message(filter.samples(DEprot.object = tb.dpo.raw,
                                samples = half.samples,
                                mode = "keep",
                                verbose = TRUE))

  expect_no_message(suppressWarnings(filter.samples(DEprot.object = tb.dpo.raw,
                                                    samples = half.samples,
                                                    mode = "keep",
                                                    verbose = FALSE)))
})
