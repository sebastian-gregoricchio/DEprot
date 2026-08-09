## Real example objects
dpo.norm   <- DEprot::test.toolbox$dpo.norm
dpo.random <- DEprot::test.toolbox$dpo.random
dpo.imp    <- DEprot::test.toolbox$dpo.imp
dpo.limma  <- DEprot::test.toolbox$diff.exp.limma


## Diagnostic on the normalized object (no randomization run before)
miss <- missingness.diagnostic(DEprot.object = dpo.norm,
                               group.column = "combined.id",
                               verbose = FALSE)


## Helper: synthetic DEprot object with a fully controlled missingness pattern.
## Groups of 4 replicates: A = s1-s4, B = s5-s8.
##   P1 -> no missing value                    (complete)
##   P2 -> missing in all the replicates of A  (MNAR)
##   P3 -> one missing value in A              (MCAR)
##   P4 -> missing everywhere                  (all.missing)
##   P5 -> one missing value in A and one in B (MCAR)
##   P6 -> missing in all the replicates of B  (MNAR)
make_synth_dpo <- function(seed = 11) {
  set.seed(seed)
  m <- matrix(stats::rnorm(6 * 8, mean = 20, sd = 2), nrow = 6,
              dimnames = list(paste0("P", 1:6), paste0("s", 1:8)))
  m[2, 1:4]    <- NA
  m[3, 3]      <- NA
  m[4, ]       <- NA
  m[5, c(1, 6)] <- NA
  m[6, 5:8]    <- NA

  meta <- data.frame(column.id = colnames(m),
                     group = rep(c("A", "B"), each = 4),
                     stringsAsFactors = FALSE)

  obj <- dpo.norm
  obj@raw.counts <- m
  obj@norm.counts <- m
  obj@random.counts <- NULL
  obj@imputed.counts <- NULL
  obj@randomization.method <- NULL
  obj@metadata <- meta
  obj@log.base <- 2
  return(obj)
}

synth <- make_synth_dpo()

## Helper: run an expression on a null device (plots are not written to disk)
on_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}


# ---------------------------------------------------------------------
# missingness.diagnostic: object and slots
# ---------------------------------------------------------------------

test_that("missingness.diagnostic returns a DEprot.missingness object", {
  expect_s4_class(miss, "DEprot.missingness")
})

test_that("missingness.diagnostic populates all the expected slots", {
  expect_true(is.character(miss@data.used))
  expect_true(is.character(miss@counts.available))
  expect_true(is.data.frame(miss@metadata))
  expect_equal(miss@group.column, "combined.id")
  expect_true(is.matrix(miss@missing.matrix))
  expect_type(miss@missing.matrix, "logical")
  expect_true(is.matrix(miss@imputation.map))
  expect_type(miss@imputation.map, "character")
  expect_true(is.data.frame(miss@protein.stats))
  expect_true(is.data.frame(miss@sample.stats))
  expect_true(is.data.frame(miss@group.summary))
  expect_true(is.data.frame(miss@pattern.summary))
  expect_true(is.list(miss@global.stats))
  expect_true(is.list(miss@plots))
  expect_true(is.list(miss@parameters))
})

test_that("the tables have the expected dimensions", {
  n.prot <- nrow(miss@missing.matrix)
  n.samp <- ncol(miss@missing.matrix)
  expect_equal(nrow(miss@protein.stats), n.prot)
  expect_equal(nrow(miss@sample.stats), n.samp)
  expect_equal(nrow(miss@pattern.summary), 4L)
  expect_equal(nrow(miss@group.summary),
               length(unique(as.character(miss@metadata$combined.id))))
})


# ---------------------------------------------------------------------
# Selection of the counts table
# ---------------------------------------------------------------------

test_that("the lowest level of counts available is used by default", {
  expect_true(miss@data.used %in% c("raw", "normalized"))
  expect_equal(miss@data.used, miss@counts.available[1])
})

test_that("randomized counts are not used when lower levels are available", {
  m <- missingness.diagnostic(dpo.random, verbose = FALSE)
  expect_true("randomized" %in% m@counts.available)
  expect_false(m@data.used == "randomized")
})

test_that("an explicit which.data is respected", {
  m <- suppressWarnings(missingness.diagnostic(dpo.random,
                                               which.data = "randomized",
                                               verbose = FALSE))
  expect_equal(m@data.used, "randomized")
})

test_that("using randomized or imputed counts raises a warning", {
  expect_warning(missingness.diagnostic(dpo.random,
                                        which.data = "randomized",
                                        verbose = FALSE))
})

test_that("a counts type that is not available raises an error", {
  expect_error(missingness.diagnostic(dpo.norm,
                                      which.data = "imputed",
                                      group.column = "combined.id",
                                      verbose = FALSE))
})

test_that("an unrecognized which.data raises an error", {
  expect_error(missingness.diagnostic(dpo.norm,
                                      which.data = "something.else",
                                      group.column = "combined.id",
                                      verbose = FALSE))
})

test_that("an object containing only imputed counts is rejected", {
  only.imp <- dpo.imp
  only.imp@raw.counts <- NULL
  only.imp@norm.counts <- NULL
  only.imp@random.counts <- NULL
  expect_error(missingness.diagnostic(only.imp, verbose = FALSE))
})


# ---------------------------------------------------------------------
# Retrieval of the randomization parameters
# ---------------------------------------------------------------------

test_that("the randomization parameters are retrieved from the object", {
  m <- missingness.diagnostic(dpo.random, verbose = FALSE)
  rand <- dpo.random@randomization.method
  expect_equal(m@parameters$group.column, rand$group.column)
  expect_equal(m@parameters$percentage.missing, rand$percentage.missing)
  expect_equal(m@parameters$tail.percentage, rand$tail.percentage)
  expect_true(grepl("randomization", m@parameters$parameters.source))
})

test_that("explicit parameters override the ones stored in the object", {
  m <- missingness.diagnostic(dpo.random,
                              group.column = "combined.id",
                              percentage.missing = 50,
                              tail.percentage = 10,
                              verbose = FALSE)
  expect_equal(m@group.column, "combined.id")
  expect_equal(m@parameters$percentage.missing, 50)
  expect_equal(m@parameters$tail.percentage, 10)
})

test_that("the default parameters are used when no randomization was run", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  expect_equal(m@parameters$percentage.missing, 100)
  expect_equal(m@parameters$tail.percentage, 3)
})

test_that("invalid parameters or group columns are rejected", {
  expect_error(missingness.diagnostic(synth, group.column = "not.a.column",
                                      verbose = FALSE))
  expect_error(missingness.diagnostic(synth, group.column = "group",
                                      percentage.missing = 150, verbose = FALSE))
  expect_error(missingness.diagnostic(synth, group.column = "group",
                                      tail.percentage = 0, verbose = FALSE))
  expect_error(missingness.diagnostic(synth, group.column = "group",
                                      tail.percentage = 100, verbose = FALSE))
})


# ---------------------------------------------------------------------
# Classification of the missing values (deterministic, synthetic)
# ---------------------------------------------------------------------

test_that("the proteins are assigned to the expected missingness class", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  cls <- stats::setNames(as.character(m@protein.stats$missing.class),
                         m@protein.stats$prot.id)
  expect_equal(unname(cls["P1"]), "complete")
  expect_equal(unname(cls["P2"]), "MNAR")
  expect_equal(unname(cls["P3"]), "MCAR")
  expect_equal(unname(cls["P4"]), "all.missing")
  expect_equal(unname(cls["P5"]), "MCAR")
  expect_equal(unname(cls["P6"]), "MNAR")
})

test_that("the group in which a protein is MNAR is reported", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  gr <- stats::setNames(m@protein.stats$groups.MNAR, m@protein.stats$prot.id)
  expect_equal(unname(gr["P2"]), "A")
  expect_equal(unname(gr["P6"]), "B")
  expect_equal(unname(gr["P1"]), "")
})

test_that("a lower percentage.missing moves proteins from MCAR to MNAR", {
  m <- missingness.diagnostic(synth, group.column = "group",
                              percentage.missing = 25, verbose = FALSE)
  cls <- stats::setNames(as.character(m@protein.stats$missing.class),
                         m@protein.stats$prot.id)
  expect_equal(unname(cls["P3"]), "MNAR")
  expect_equal(unname(cls["P5"]), "MNAR")
  expect_equal(unname(cls["P1"]), "complete")
})

test_that("the imputation map counts the cells of each type", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  expect_true(all(m@imputation.map %in% c("detected", "MNAR", "MCAR")))
  expect_equal(sum(m@imputation.map == "MNAR"), 16L)   # P2 (4) + P4 (8) + P6 (4)
  expect_equal(sum(m@imputation.map == "MCAR"), 3L)    # P3 (1) + P5 (2)
  expect_equal(sum(m@imputation.map != "detected"), sum(m@missing.matrix))
})

test_that("the missing matrix matches the NAs of the counts used", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  expect_equal(unname(m@missing.matrix), unname(is.na(synth@raw.counts)))
})

test_that("proteins never detected have no intensity estimate", {
  m <- missingness.diagnostic(synth, group.column = "group", verbose = FALSE)
  expect_true(is.na(m@protein.stats$mean.intensity[m@protein.stats$prot.id == "P4"]))
  expect_false(is.na(m@protein.stats$mean.intensity[m@protein.stats$prot.id == "P1"]))
})


# ---------------------------------------------------------------------
# Internal consistency of the summaries
# ---------------------------------------------------------------------

test_that("the per-class counts sum to the number of proteins", {
  expect_equal(sum(miss@pattern.summary$n.proteins), nrow(miss@protein.stats))
  expect_equal(sum(miss@pattern.summary$percentage), 100)
})

test_that("the per-group counts sum to the number of proteins", {
  tot <- miss@group.summary$proteins.complete +
    miss@group.summary$proteins.MCAR +
    miss@group.summary$proteins.MNAR
  expect_true(all(tot == nrow(miss@protein.stats)))
})

test_that("the global statistics agree with the missing matrix", {
  expect_equal(miss@global.stats$n.missing, sum(miss@missing.matrix))
  expect_equal(miss@global.stats$n.proteins, nrow(miss@missing.matrix))
  expect_equal(miss@global.stats$n.samples, ncol(miss@missing.matrix))
  expect_equal(miss@global.stats$perc.missing,
               100 * sum(miss@missing.matrix) / length(miss@missing.matrix))
  expect_equal(miss@global.stats$n.cells.MNAR + miss@global.stats$n.cells.MCAR,
               sum(miss@missing.matrix))
})

test_that("the per-sample counts agree with the missing matrix", {
  expect_equal(miss@sample.stats$n.missing, unname(colSums(miss@missing.matrix)))
  expect_equal(miss@sample.stats$n.detected, unname(colSums(!miss@missing.matrix)))
  expect_true(all(miss@sample.stats$perc.missing >= 0 &
                    miss@sample.stats$perc.missing <= 100))
})

test_that("the per-protein counts agree with the missing matrix", {
  expect_equal(miss@protein.stats$n.missing, unname(rowSums(miss@missing.matrix)))
  expect_equal(miss@protein.stats$n.detected + miss@protein.stats$n.missing,
               rep(ncol(miss@missing.matrix), nrow(miss@missing.matrix)))
})


# ---------------------------------------------------------------------
# Dropout model
# ---------------------------------------------------------------------

test_that("the dropout model is a binomial glm", {
  expect_s3_class(miss@dropout.model, "glm")
  expect_equal(stats::family(miss@dropout.model)$family, "binomial")
})

test_that("the dropout statistics are coherent with the fitted model", {
  coefs <- stats::coef(miss@dropout.model)
  expect_equal(miss@global.stats$dropout.slope, unname(coefs[2]))
  if (!is.na(miss@global.stats$LOD50)) {
    expect_equal(miss@global.stats$LOD50, unname(-coefs[1] / coefs[2]))
    rng <- range(miss@protein.stats$mean.intensity, na.rm = TRUE)
    expect_gte(miss@global.stats$LOD50, rng[1])
    expect_lte(miss@global.stats$LOD50, rng[2])
  }
})

test_that("the tail threshold is the quantile of the counts used", {
  expect_equal(miss@global.stats$tail.threshold, miss@parameters$tail.threshold)
  expect_true(is.numeric(miss@global.stats$tail.threshold))
  expect_length(miss@global.stats$tail.threshold, 1L)
})


# ---------------------------------------------------------------------
# Jaccard similarity and clustering
# ---------------------------------------------------------------------

test_that("the Jaccard matrix is square, symmetric and bounded", {
  jm <- miss@jaccard.matrix
  expect_true(is.matrix(jm))
  expect_equal(nrow(jm), ncol(jm))
  expect_equal(nrow(jm), ncol(miss@missing.matrix))
  expect_equal(jm, t(jm))
  expect_true(all(diag(jm) == 1))
  expect_true(all(jm >= 0 & jm <= 1))
})

test_that("the samples are clustered on the Jaccard dissimilarity", {
  expect_s3_class(miss@jaccard.cluster, "hclust")
  expect_length(miss@jaccard.cluster$order, ncol(miss@missing.matrix))
  expect_true(all(miss@jaccard.cluster$labels %in% colnames(miss@missing.matrix)))
})


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

test_that("the expected plots are generated", {
  expect_named(miss@plots,
               c("detection.density", "dropout.curve", "missingness.heatmap",
                 "missing.per.sample", "detection.frequency", "pattern.barplot",
                 "sample.similarity", "upset"))
  expect_s3_class(miss@plots$detection.density, "ggplot")
  expect_s3_class(miss@plots$dropout.curve, "ggplot")
  expect_s3_class(miss@plots$missing.per.sample, "ggplot")
  expect_s3_class(miss@plots$detection.frequency, "ggplot")
  expect_s3_class(miss@plots$pattern.barplot, "ggplot")
  expect_s3_class(miss@plots$sample.similarity, "ggplot")
})

test_that("the missingness heatmap is present when missing values exist", {
  if (sum(miss@missing.matrix) > 0) {
    expect_s3_class(miss@plots$missingness.heatmap, "ggplot")
  } else {
    expect_null(miss@plots$missingness.heatmap)
  }
})

test_that("the upset plot is either a plot or NULL", {
  expect_true(is.null(miss@plots$upset) ||
                inherits(miss@plots$upset, c("ggplot", "patchwork", "upset")))
})

test_that("the number of proteins in the heatmap can be capped", {
  m <- missingness.diagnostic(dpo.norm,
                              group.column = "combined.id",
                              max.proteins.heatmap = 10,
                              verbose = FALSE)
  n.rows <- length(unique(as.character(m@plots$missingness.heatmap$data$prot.id)))
  expect_lte(n.rows, 10L)
})


# ---------------------------------------------------------------------
# Sample subset
# ---------------------------------------------------------------------

test_that("the analysis can be restricted to a subset of samples", {
  keep <- colnames(dpo.norm@norm.counts)[1:4]
  m <- missingness.diagnostic(dpo.norm,
                              group.column = "combined.id",
                              sample.subset = keep,
                              verbose = FALSE)
  expect_equal(ncol(m@missing.matrix), 4L)
  expect_equal(nrow(m@sample.stats), 4L)
  expect_true(all(as.character(m@sample.stats$column.id) %in% keep))
})

test_that("an unknown sample raises an error", {
  expect_error(missingness.diagnostic(dpo.norm,
                                      group.column = "combined.id",
                                      sample.subset = c("not.a.sample"),
                                      verbose = FALSE))
})


# ---------------------------------------------------------------------
# Contrast-level diagnostics
# ---------------------------------------------------------------------

test_that("no contrast-level result is produced from a plain DEprot object", {
  expect_null(miss@contrast.stats)
})

test_that("the contrasts of a DEprot.analyses object are analyzed", {
  skip_if(identical(sort(names(which(!vapply(
    c("raw.counts", "norm.counts", "random.counts"),
    function(s) !is.null(methods::slot(dpo.limma, s)), logical(1))))),
    c("norm.counts", "random.counts", "raw.counts")),
    "no unimputed counts available in the test object")

  m <- missingness.diagnostic(dpo.limma, verbose = FALSE)
  expect_true(is.list(m@contrast.stats))
  expect_length(m@contrast.stats, length(dpo.limma@contrasts))
  expect_named(m@contrast.stats, names(dpo.limma@contrasts))

  ct <- m@contrast.stats[[1]]
  expect_true(is.data.frame(ct$protein.stats))
  expect_true(is.data.frame(ct$summary))
  expect_type(ct$protein.stats$testable, "logical")
  expect_true(all(c("group.1", "group.2") %in% names(ct$samples)))
  expect_s3_class(ct$plots$pattern.barplot, "ggplot")
  expect_s3_class(ct$plots$detection.density, "ggplot")

  ## the classes are directional and the proteins missing everywhere are not testable
  lev <- levels(ct$protein.stats$missing.class)
  expect_true(any(grepl("^MNAR in ", lev)))
  expect_equal(sum(ct$summary$n.proteins), nrow(ct$protein.stats))
  expect_false(any(ct$protein.stats$testable[ct$protein.stats$missing.class == "all.missing"]))
})

test_that("the contrast-level analyses can be skipped", {
  m <- missingness.diagnostic(dpo.limma, contrasts = "none", verbose = FALSE)
  expect_null(m@contrast.stats)
})

test_that("a subset of contrasts can be selected", {
  skip_if(length(dpo.limma@contrasts) < 2, "less than two contrasts available")
  m <- missingness.diagnostic(dpo.limma, contrasts = 2, verbose = FALSE)
  expect_length(m@contrast.stats, 1L)
  expect_named(m@contrast.stats, names(dpo.limma@contrasts)[2])
})

test_that("an unavailable contrast raises an error", {
  expect_error(missingness.diagnostic(dpo.limma, contrasts = 99, verbose = FALSE))
  expect_error(missingness.diagnostic(dpo.limma, contrasts = "not.a.contrast",
                                      verbose = FALSE))
})


# ---------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------

test_that("the summary method returns a table", {
  s <- getMethod("summary", "DEprot.missingness")(miss)
  expect_true(is.data.frame(s))
  expect_equal(nrow(s), nrow(miss@group.summary))
})

test_that("the plot method returns the combined plots", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  p <- getMethod("plot", "DEprot.missingness")(miss)
  expect_s3_class(p, "patchwork")
})

test_that("the plot method returns a single plot when one is requested", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  p <- getMethod("plot", "DEprot.missingness")(miss, plot.type = "dropout")
  expect_s3_class(p, "ggplot")
})

