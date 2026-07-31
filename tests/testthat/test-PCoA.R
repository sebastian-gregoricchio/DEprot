## Real example object
dpo <- DEprot::test.toolbox$dpo.imp

## Correlation object from the real data (spearman)
corr <- plot.correlation.heatmap(DEprot.object = dpo,
                                 correlation.method = "spearman")

n.samples <- nrow(corr@corr.matrix)

## PCoA on the real data, WITHOUT the DEprot object (no protein projections)
pcoa <- perform.PCoA(DEprot.correlation.object = corr, verbose = FALSE)

## PCoA on the real data, WITH the DEprot object (protein projections available)
pcoa.load <- perform.PCoA(DEprot.correlation.object = corr,
                          DEprot.object = dpo,
                          verbose = FALSE)


## Helper: build a small synthetic DEprot.correlation with GUARANTEED
## negative eigenvalues (strong anticorrelation makes 1 - r non-euclidean).
make_synth_corr <- function(n.features = 80, n.samples = 7, seed = 7) {
  set.seed(seed)
  m <- matrix(rnorm(n.features * n.samples), nrow = n.features,
              dimnames = list(paste0("PROT", seq_len(n.features)),
                              paste0("s", seq_len(n.samples))))
  # force a couple of near-perfect anticorrelations
  m[, 2] <- -m[, 1] + rnorm(n.features, sd = 0.15)
  if (n.samples >= 5) m[, 5] <- -m[, 4] + rnorm(n.features, sd = 0.15)

  cm <- stats::cor(m, method = "pearson")
  meta <- data.frame(column.id = colnames(m),
                     condition = rep(c("a", "b"), length.out = n.samples),
                     replicate = paste0("rep", seq_len(n.samples)),
                     stringsAsFactors = FALSE)

  methods::new(Class = "DEprot.correlation",
               heatmap = list(labels = list(title = "**Pearson correlation**")),
               corr.metadata = meta,
               sample.subset = colnames(m),
               data.used = "imputed",
               corr.matrix = cm,
               distance = stats::as.dist(1 - cm),
               cluster = NULL,
               method = "pearson")
}

## Helper: number of arrows drawn in a biplot (robust to layer ordering)
n_arrows <- function(p) {
  hits <- vapply(p$layers, function(l) {
    d <- l$data
    is.data.frame(d) && all(c("variable", "loading.x.scaled") %in% names(d))
  }, logical(1))
  if (!any(hits)) return(0L)
  nrow(p$layers[[which(hits)[1]]]$data)
}


# ---------------------------------------------------------------------
# perform.PCoA: object and slots
# ---------------------------------------------------------------------

test_that("perform.PCoA returns a DEprot.PCoA object", {
  expect_s4_class(pcoa, "DEprot.PCoA")
})

test_that("perform.PCoA populates all the expected slots", {
  expect_true(is.data.frame(pcoa@PCos))
  expect_s3_class(pcoa@distance, "dist")
  expect_s3_class(pcoa@distance.input, "dist")
  expect_type(pcoa@eigenvalues, "double")
  expect_true(is.data.frame(pcoa@importance))
  expect_true(is.list(pcoa@euclidean.diagnostics))
  expect_true(is.list(pcoa@parameters))
  expect_false(is.null(pcoa@cmdscale))
  # stored plots
  expect_s3_class(pcoa@scatter.plot, "ggplot")
  expect_s3_class(pcoa@cumulative.PCo.plot, "ggplot")
  expect_s3_class(pcoa@shepard.plot, "ggplot")
})

test_that("the scores table has one row per sample and the coordinate columns", {
  expect_equal(nrow(pcoa@PCos), n.samples)
  expect_true("PCo1" %in% colnames(pcoa@PCos))
  expect_true("PCo2" %in% colnames(pcoa@PCos))
  expect_true(all(c("column.id", "condition", "replicate") %in% colnames(pcoa@PCos)))
})

test_that("the combined scatter is present when at least 3 axes exist", {
  if (pcoa@parameters$n.axes >= 3) {
    expect_s3_class(pcoa@scatter.plot.123, "patchwork")
  } else {
    expect_null(pcoa@scatter.plot.123)
  }
})


# ---------------------------------------------------------------------
# perform.PCoA: correlation-method recovery (deterministic, synthetic)
# ---------------------------------------------------------------------

test_that("the correlation method is read from the 'method' slot", {
  sc <- make_synth_corr()
  sc@method <- "spearman"
  sc@heatmap$labels$title <- "**Pearson correlation**"   # title disagrees on purpose
  p <- perform.PCoA(sc, verbose = FALSE)
  expect_equal(p@correlation.method, "spearman")          # slot wins over title
})

test_that("an explicit correlation.method argument overrides the slot", {
  sc <- make_synth_corr()
  sc@method <- "pearson"
  p <- perform.PCoA(sc, correlation.method = "kendall", verbose = FALSE)
  expect_equal(p@correlation.method, "kendall")
})

test_that("the method falls back to the title, then to 'unknown'", {
  sc <- make_synth_corr()
  sc@method <- NA
  sc@heatmap$labels$title <- "**Spearman correlation**"
  expect_equal(perform.PCoA(sc, verbose = FALSE)@correlation.method, "spearman")

  sc@heatmap$labels$title <- "**correlation**"            # uninformative
  expect_equal(perform.PCoA(sc, verbose = FALSE)@correlation.method, "unknown")
})

test_that("an unsupported correlation.method argument errors", {
  expect_error(perform.PCoA(corr, correlation.method = "kendal"),
               "not supported")
})


# ---------------------------------------------------------------------
# perform.PCoA: number of axes and axis validation
# ---------------------------------------------------------------------

test_that("by default all computable coordinates are returned", {
  # with non-euclidean dissimilarities (1 - correlation) the ordination keeps only the
  # coordinates that correspond to a positive eigenvalue, so their number can be lower
  # than n_samples - 1; it can never exceed it.
  expect_true(pcoa@parameters$n.axes <= n.samples - 1)
  expect_true(pcoa@parameters$n.axes >= 2)

  # the number of coordinate columns matches the reported number of axes
  pco.cols <- grep("^PCo[0-9]+$", colnames(pcoa@PCos), value = TRUE)
  expect_equal(length(pco.cols), pcoa@parameters$n.axes)
})

test_that("a custom number of axes is respected", {
  p <- perform.PCoA(corr, n.axes = 3, verbose = FALSE)
  expect_equal(p@parameters$n.axes, 3)
  expect_true("PCo3" %in% colnames(p@PCos))
  expect_false("PCo4" %in% colnames(p@PCos))
})

test_that("requesting fewer than 2 axes errors", {
  expect_error(perform.PCoA(corr, n.axes = 1, verbose = FALSE),
               "At least 2 axes")
})

test_that("identical PCo.x and PCo.y error", {
  expect_error(perform.PCoA(corr, PCo.x = 2, PCo.y = 2, verbose = FALSE),
               "different principal coordinates")
})

test_that("scatter axes beyond the computed coordinates error", {
  expect_error(perform.PCoA(corr, n.axes = 3, PCo.x = 5, PCo.y = 1, verbose = FALSE),
               "exceed")
})


# ---------------------------------------------------------------------
# perform.PCoA: distance transformation and corrections
# ---------------------------------------------------------------------

test_that("distance.transformation is recorded in distance.method", {
  expect_match(perform.PCoA(corr, verbose = FALSE)@distance.method,
               "1 - correlation.matrix")
  expect_match(perform.PCoA(corr, distance.transformation = "sqrt", verbose = FALSE)@distance.method,
               "sqrt")
})

test_that("the sqrt transformation removes negative eigenvalues", {
  sc <- make_synth_corr()
  none <- perform.PCoA(sc, verbose = FALSE)
  sqrt <- perform.PCoA(sc, distance.transformation = "sqrt", verbose = FALSE)
  expect_gt(none@euclidean.diagnostics$n.negative.eigenvalues, 0)   # guaranteed non-euclidean
  expect_equal(sqrt@euclidean.diagnostics$n.negative.eigenvalues, 0)
  expect_true(sqrt@euclidean.diagnostics$euclidean.embeddable)
})

test_that("the Cailliez correction removes negative eigenvalues", {
  sc <- make_synth_corr()
  ca <- perform.PCoA(sc, correction = "cailliez", verbose = FALSE)
  expect_equal(ca@euclidean.diagnostics$n.negative.eigenvalues, 0)
  expect_match(ca@distance.method, "Cailliez")
})

test_that("the Lingoes correction removes negative eigenvalues", {
  sc <- make_synth_corr()
  li <- perform.PCoA(sc, correction = "lingoes", verbose = FALSE)
  expect_equal(li@euclidean.diagnostics$n.negative.eigenvalues, 0)
  expect_true(li@euclidean.diagnostics$euclidean.embeddable)
  expect_match(li@distance.method, "Lingoes")
})

test_that("unsupported transformation or correction error", {
  expect_error(perform.PCoA(corr, distance.transformation = "log", verbose = FALSE),
               "transformation")
  expect_error(perform.PCoA(corr, correction = "banana", verbose = FALSE),
               "correction")
})


# ---------------------------------------------------------------------
# perform.PCoA: importance table
# ---------------------------------------------------------------------

test_that("the importance table is coherent", {
  imp <- pcoa@importance
  expect_true(all(c("Eigenvalue", "Proportion.of.Variance", "Cumulative.Proportion",
                    "Percentage.of.Variance", "Broken.stick", "Positive.eigenvalue") %in% colnames(imp)))

  # proportions of the positive eigenvalues sum to 1
  pos <- imp$Positive.eigenvalue
  expect_equal(sum(imp$Proportion.of.Variance[pos]), 1, tolerance = 1e-8)

  # cumulative proportion is monotonically non-decreasing over the positive eigenvalues
  # (with non-euclidean dissimilarities the negative eigenvalues make it decrease afterwards)
  expect_true(all(diff(imp$Cumulative.Proportion[pos]) >= -1e-8))

  # eigenvalues are sorted in decreasing order
  expect_true(all(diff(imp$Eigenvalue) <= 1e-8))
})

test_that("the naive proportion differs from the positive-only proportion when negatives exist", {
  sc <- make_synth_corr()
  p <- perform.PCoA(sc, verbose = FALSE)
  naive <- p@eigenvalues[1] / sum(p@eigenvalues)
  correct <- p@importance$Proportion.of.Variance[1]
  expect_false(isTRUE(all.equal(naive, correct)))
})


# ---------------------------------------------------------------------
# perform.PCoA: protein projections (axis.loadings)
# ---------------------------------------------------------------------

test_that("axis.loadings is NULL when no DEprot object is provided", {
  expect_null(pcoa@axis.loadings)
})

test_that("axis.loadings is a proper table when the DEprot object is provided", {
  expect_true(is.data.frame(pcoa.load@axis.loadings))
  expect_true("protein" %in% colnames(pcoa.load@axis.loadings))
  expect_true(all(c("PCo1", "PCo2") %in% colnames(pcoa.load@axis.loadings)))
  # projections are correlations, hence bounded in [-1, 1]
  vals <- as.matrix(pcoa.load@axis.loadings[, c("PCo1", "PCo2")])
  expect_true(all(vals >= -1 - 1e-8 & vals <= 1 + 1e-8, na.rm = TRUE))
})


# ---------------------------------------------------------------------
# perform.PCoA: euclidean diagnostics
# ---------------------------------------------------------------------

test_that("euclidean diagnostics expose the expected fields", {
  d <- pcoa@euclidean.diagnostics
  expect_true(all(c("n.negative.eigenvalues", "euclidean.embeddable",
                    "GOF", "shepard.correlation") %in% names(d)))
  expect_type(d$euclidean.embeddable, "logical")
  expect_true(is.numeric(d$shepard.correlation))
})


# ---------------------------------------------------------------------
# perform.PCoA: input validation
# ---------------------------------------------------------------------

test_that("perform.PCoA rejects a non-DEprot.correlation input", {
  expect_error(perform.PCoA("not an object"), "DEprot.correlation")
  expect_error(perform.PCoA(dpo), "DEprot.correlation")
})


# ---------------------------------------------------------------------
# plot.PCoA.scatter
# ---------------------------------------------------------------------

test_that("plot.PCoA.scatter returns a ggplot", {
  expect_s3_class(plot.PCoA.scatter(pcoa, PCo.x = 1, PCo.y = 2), "ggplot")
})

test_that("plot.PCoA.scatter honours color/shape/label columns", {
  p <- plot.PCoA.scatter(pcoa, PCo.x = 1, PCo.y = 2,
                         color.column = "condition",
                         shape.column = "replicate",
                         label.column = "replicate")
  expect_s3_class(p, "ggplot")
})

test_that("plot.PCoA.scatter axis labels report the coordinate and the variance", {
  p <- plot.PCoA.scatter(pcoa, PCo.x = 1, PCo.y = 2)
  expect_match(p$labels$x, "PCo1")
  expect_match(p$labels$x, "%")
  expect_match(p$labels$y, "PCo2")
})

test_that("plot.PCoA.scatter rejects a wrong object class", {
  expect_error(plot.PCoA.scatter(pcoa.load@axis.loadings), "DEprot.PCoA")
})

test_that("plot.PCoA.scatter errors on non-existent aesthetic columns", {
  expect_error(plot.PCoA.scatter(pcoa, color.column = "nope"), "color")
  expect_error(plot.PCoA.scatter(pcoa, shape.column = "nope"), "shape")
  expect_error(plot.PCoA.scatter(pcoa, label.column = "nope"), "label")
})

test_that("plot.PCoA.scatter errors on unavailable axes", {
  expect_error(plot.PCoA.scatter(pcoa, PCo.x = 999, PCo.y = 1))
})


# ---------------------------------------------------------------------
# plot.PCoA.scatter.123
# ---------------------------------------------------------------------

test_that("plot.PCoA.scatter.123 returns a patchwork", {
  skip_if_not(pcoa@parameters$n.axes >= 3, "fewer than 3 axes available")
  expect_s3_class(plot.PCoA.scatter.123(pcoa, color.column = "condition"), "patchwork")
})

test_that("plot.PCoA.scatter.123 rejects a wrong object class", {
  expect_error(plot.PCoA.scatter.123("x"), "DEprot.PCoA")
})

test_that("plot.PCoA.scatter.123 errors when fewer than 3 axes are available", {
  p2 <- perform.PCoA(corr, n.axes = 2, verbose = FALSE)
  expect_error(plot.PCoA.scatter.123(p2), "3 axes")
})


# ---------------------------------------------------------------------
# plot.PCoA.cumulative
# ---------------------------------------------------------------------

test_that("plot.PCoA.cumulative returns a ggplot", {
  expect_s3_class(plot.PCoA.cumulative(pcoa), "ggplot")
  expect_s3_class(plot.PCoA.cumulative(pcoa, broken.stick = FALSE), "ggplot")
})

test_that("plot.PCoA.cumulative rejects a wrong object class", {
  expect_error(plot.PCoA.cumulative("x"), "DEprot.PCoA")
})


# ---------------------------------------------------------------------
# plot.PCoA.biplot
# ---------------------------------------------------------------------

test_that("plot.PCoA.biplot returns a ggplot when projections are available", {
  expect_s3_class(plot.PCoA.biplot(pcoa.load, PCo.x = 1, PCo.y = 2,
                                   color.column = "condition"), "ggplot")
})

test_that("plot.PCoA.biplot draws the requested number of arrows", {
  p <- plot.PCoA.biplot(pcoa.load, PCo.x = 1, PCo.y = 2, n.loadings = 8)
  expect_equal(n_arrows(p), 8)
})

test_that("plot.PCoA.biplot errors when the axis.loadings slot is empty", {
  expect_error(plot.PCoA.biplot(pcoa), "axis.loadings")   # pcoa built without DEprot object
})

test_that("plot.PCoA.biplot rejects a wrong object class", {
  expect_error(plot.PCoA.biplot("x"), "DEprot.PCoA")
})

test_that("plot.PCoA.biplot errors on unavailable axes", {
  expect_error(plot.PCoA.biplot(pcoa.load, PCo.x = 999, PCo.y = 1))
})

test_that("plot.PCoA.biplot min.abs.correlation filtering behaves", {
  # a threshold above 1 removes every protein
  expect_error(plot.PCoA.biplot(pcoa.load, min.abs.correlation = 1.1),
               "min.abs.correlation")
  # a mild threshold still returns a plot
  expect_s3_class(plot.PCoA.biplot(pcoa.load, min.abs.correlation = 0.1,
                                   n.loadings = 5), "ggplot")
})


# ---------------------------------------------------------------------
# DEprot.PCoA class and show-method
# ---------------------------------------------------------------------

test_that("the DEprot.PCoA class is valid and the show-method prints", {
  expect_true(isVirtualClass("DEprot.PCoA") == FALSE)
  expect_true(isTRUE(methods::validObject(pcoa, test = TRUE)))
  expect_output(print(pcoa), "DEprot.PCoA object")
  expect_output(print(pcoa), "Principal coordinates")
})


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

test_that(".ordination.slots harmonizes PCA and PCoA and rejects others", {
  s <- DEprot:::.ordination.slots(pcoa)
  expect_equal(s$axis.prefix, "PCo")
  expect_true(is.data.frame(s$coordinates))
  expect_error(DEprot:::.ordination.slots(list()), "DEprot.PCA")
})

test_that(".check.aes.columns validates the requested columns", {
  expect_error(DEprot:::.check.aes.columns(pcoa@PCos, color.column = "zzz",
                                           axis.prefix = "PCo"), "color")
  ok <- DEprot:::.check.aes.columns(pcoa@PCos, color.column = "condition",
                                    shape.column = "replicate", axis.prefix = "PCo")
  expect_true(ok$show.labels == FALSE)
})

test_that(".scale.loadings ranks and rescales the loadings", {
  ld <- data.frame(variable = paste0("P", 1:20),
                   loading.x = rnorm(20),
                   loading.y = rnorm(20))
  res <- DEprot:::.scale.loadings(ld, scores.x = rnorm(6), scores.y = rnorm(6),
                                  n.loadings = 5)
  expect_equal(nrow(res$top), 5)
  expect_true(is.numeric(res$scale.factor))
  # the top loadings are the ones with the largest distance from origin
  expect_equal(res$top$variable,
               ld$variable[order(sqrt(ld$loading.x^2 + ld$loading.y^2),
                                 decreasing = TRUE)][1:5])
})

test_that(".lingoes.correction returns a euclidean-embeddable dist", {
  sc <- make_synth_corr()
  d.corrected <- DEprot:::.lingoes.correction(sc@distance)
  expect_s3_class(d.corrected, "dist")
  # after correction the Gower matrix has no negative eigenvalues
  D <- as.matrix(d.corrected)
  n <- nrow(D)
  J <- diag(n) - matrix(1/n, n, n)
  G <- -0.5 * (J %*% (D^2) %*% J)
  ev <- eigen((G + t(G))/2, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) > -1e-8)
})
