## ----------------------------------------------------------------------------------------
##  Argument branches of the visualization functions
##  'test-visualizations.R' calls each of them once with the default arguments: what is left
##  uncovered are the optional branches, i.e. most of the body of these functions.
##  Building the plot is not enough to walk through every branch: 'ggplot_build' forces the
##  evaluation of the layers, hence the plots are built here rather than only created.
## ----------------------------------------------------------------------------------------

build.plot <-
  function(p) {
    if (inherits(p, "ggplot")) {return(suppressWarnings(suppressMessages(ggplot2::ggplot_build(p))))}
    return(p)
  }

proteins <- rownames(tb.dpo.imp@imputed.counts)[1:3]



##################################    expression.boxplot    ##############################

test_that("a single protein is plotted", {
  for (prot in proteins) {
    p <- expression.boxplot(DEprot.object = tb.dpo.imp, protein.id = prot)
    expect_true(inherits(p, "ggplot") | inherits(p, "patchwork"))
    expect_no_error(build.plot(p))
  }
})


test_that("several proteins are plotted in one panel each", {
  p <- expression.boxplot(DEprot.object = tb.dpo.imp, protein.id = proteins)

  expect_true(inherits(p, "ggplot"))
  expect_no_error(build.plot(p))

  ## one panel per protein, in the order in which they were requested
  panels <- levels(ggplot2::ggplot_build(p)$data[[1]]$PANEL)
  expect_equal(length(panels), length(proteins))
  expect_equal(levels(p$data$prot.id), proteins)

  ## the layout arguments are honoured
  expect_no_error(build.plot(expression.boxplot(DEprot.object = tb.dpo.imp,
                                                protein.id = proteins,
                                                ncol = 1,
                                                free.y = FALSE)))
})


test_that("the Z-score is computed protein by protein", {
  p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                          protein.id = proteins,
                          scale.expression = TRUE)

  ## each panel is centred on its own mean, hence every protein averages to zero
  means <- tapply(p$data$expression, p$data$prot.id, mean, na.rm = TRUE)
  expect_true(all(abs(means) < 1e-8))
})


test_that("the proteins requested must exist", {
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp,
                                  protein.id = c(proteins[1], "not.a.protein")))
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp, protein.id = character(0)))
})


test_that("the pairwise comparisons are computed within each panel", {
  skip_if_not_installed("ggpubr")

  p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                          protein.id = proteins[1:2],
                          group.by.metadata.column = "condition",
                          pairwise.comparisons = TRUE,
                          pairwise.p.label = "p.value")

  expect_no_error(build.plot(p))
})


test_that("the counts type and the sample subset are honoured", {
  for (w in c("raw", "normalized", "imputed")) {
    p <- expression.boxplot(DEprot.object = tb.dpo.imp, protein.id = proteins[1], which.data = w)
    expect_no_error(build.plot(p))
  }

  subset.p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                                 protein.id = proteins[1],
                                 sample.subset = tb.dpo.imp@metadata$column.id[1:6])
  expect_no_error(build.plot(subset.p))
})


test_that("the grouping, the levels and the shape column are honoured", {
  p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                          protein.id = proteins[1],
                          group.by.metadata.column = "condition",
                          shape.column = "replicate")

  expect_no_error(build.plot(p))

  levels.p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                                 protein.id = proteins[1],
                                 group.by.metadata.column = "condition",
                                 group.levels = rev(sort(unique(as.character(tb.dpo.imp@metadata$condition)))))

  expect_no_error(build.plot(levels.p))
})


test_that("the expression can be scaled and the labels rotated", {
  p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                          protein.id = proteins[1],
                          scale.expression = TRUE,
                          x.label.angle = 90)

  expect_no_error(build.plot(p))
})


test_that("the pairwise comparisons are computed", {
  skip_if_not_installed("ggpubr")

  p <- expression.boxplot(DEprot.object = tb.dpo.imp,
                          protein.id = proteins[1],
                          group.by.metadata.column = "condition",
                          pairwise.comparisons = TRUE,
                          pairwise.test.type = "wilcox.test",
                          pairwise.include.ns = FALSE,
                          pairwise.p.label = "p.format",
                          pairwise.p.decimals = 3)

  expect_no_error(build.plot(p))

  t.based <- expression.boxplot(DEprot.object = tb.dpo.imp,
                                protein.id = proteins[1],
                                group.by.metadata.column = "condition",
                                pairwise.comparisons = TRUE,
                                pairwise.test.type = "t.test")

  expect_no_error(build.plot(t.based))
})


test_that("wrong inputs are rejected", {
  expect_error(expression.boxplot(DEprot.object = DEprot::sample.config, protein.id = proteins[1]))
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp, protein.id = "not.a.protein"))
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp,
                                  protein.id = proteins[1],
                                  which.data = "not.a.type"))
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp,
                                  protein.id = proteins[1],
                                  group.by.metadata.column = "not.a.column"))
  expect_error(expression.boxplot(DEprot.object = tb.dpo.imp,
                                  protein.id = proteins[1],
                                  sample.subset = "not.a.sample"))
})



#####################################    heatmap.counts    ###############################

test_that("the heatmap is generated with the default arguments", {
  hm <- heatmap.counts(DEprot.object = tb.dpo.imp, which.data = "imputed")

  expect_s4_class(hm, "DEprot.counts.heatmap")
  expect_no_error(build.plot(hm@heatmap))
})


test_that("the protein and sample subsets are honoured", {
  hm <- heatmap.counts(DEprot.object = tb.dpo.imp,
                       which.data = "imputed",
                       protein.subset = rownames(tb.dpo.imp@imputed.counts)[1:20],
                       sample.subset = tb.dpo.imp@metadata$column.id[1:6],
                       show.protein.names = TRUE,
                       protein.names.pattern = "^prot")

  expect_s4_class(hm, "DEprot.counts.heatmap")
  expect_no_error(build.plot(hm@heatmap))
})


test_that("the clustering can be switched off and tuned", {
  no.clust <- heatmap.counts(DEprot.object = tb.dpo.imp,
                             which.data = "imputed",
                             clust.rows = FALSE,
                             clust.columns = FALSE)

  tuned <- heatmap.counts(DEprot.object = tb.dpo.imp,
                          which.data = "imputed",
                          distance.method = "manhattan",
                          clustering.method = "average")

  expect_no_error(build.plot(no.clust@heatmap))
  expect_no_error(build.plot(tuned@heatmap))
})


test_that("the scaling and the colors are honoured", {
  for (sc in list(NULL, "row", "column")) {
    hm <- suppressWarnings(heatmap.counts(DEprot.object = tb.dpo.imp,
                                          which.data = "imputed",
                                          scale = sc))
    expect_s4_class(hm, "DEprot.counts.heatmap")
  }

  colored <- heatmap.counts(DEprot.object = tb.dpo.imp,
                            which.data = "imputed",
                            scale = "row",
                            high.color = "darkred",
                            low.color = "navy",
                            mid.color = "grey95",
                            na.color = "black",
                            color.limits = c(-2, 2),
                            cell.border.color = "white",
                            cell.border.width = 0.2,
                            title = "custom title")

  expect_no_error(build.plot(colored@heatmap))
})


test_that("the heatmap can be restricted to the proteins of a contrast", {
  hm <- heatmap.counts(DEprot.object = tb.limma,
                       which.data = "imputed",
                       contrast = 1,
                       top.n = 10,
                       group.by.metadata.column = "condition")

  expect_s4_class(hm, "DEprot.counts.heatmap")
  expect_no_error(build.plot(hm@heatmap))

  uncorrected <- heatmap.counts(DEprot.object = tb.limma,
                                which.data = "imputed",
                                contrast = 1,
                                use.uncorrected.pvalue = TRUE)

  expect_s4_class(uncorrected, "DEprot.counts.heatmap")
})


test_that("wrong inputs are rejected", {
  expect_error(heatmap.counts(DEprot.object = DEprot::sample.config))
  expect_error(heatmap.counts(DEprot.object = tb.dpo.imp, which.data = "not.a.type"))
  expect_error(heatmap.counts(DEprot.object = tb.dpo.imp, sample.subset = "not.a.sample"))
  expect_error(heatmap.counts(DEprot.object = tb.limma, contrast = 100))
  ## a contrast can only be requested on an object carrying differential analyses
  expect_error(heatmap.counts(DEprot.object = tb.dpo.imp, contrast = 1))
})



################################    plot.correlation.heatmap    ##########################

test_that("the three correlation methods are accepted", {
  for (m in c("pearson", "spearman", "kendall")) {
    corr <- plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                     which.data = "imputed",
                                     correlation.method = m)

    expect_s4_class(corr, "DEprot.correlation")
    expect_equal(dim(corr@corr.matrix), rep(nrow(tb.dpo.imp@metadata), 2))
  }
})


test_that("the dendrogram can be moved and styled", {
  ## 'left'/'right' go through legendry::scale_y_dendro
  for (pos in c("left", "right")) {
    corr <- suppressWarnings(plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                                      which.data = "imputed",
                                                      dendrogram.position = pos,
                                                      dendrogram.color = "steelblue",
                                                      dendrogram.linewidth = 1))
    expect_s4_class(corr, "DEprot.correlation")
  }
})


test_that("the dendrogram can be placed on the x axis", {
  ## 'top'/'bottom' go through legendry::scale_x_dendro
  for (pos in c("top", "bottom")) {
    corr <- suppressWarnings(plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                                      which.data = "imputed",
                                                      dendrogram.position = pos,
                                                      dendrogram.color = "steelblue",
                                                      dendrogram.linewidth = 1))
    expect_s4_class(corr, "DEprot.correlation")
    expect_no_error(build.plot(corr@heatmap))
  }
})


test_that("the values displayed on the cells are configurable", {
  shown <- plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                    which.data = "imputed",
                                    display.values = TRUE,
                                    values.color = "black",
                                    values.decimals = 3,
                                    values.font.size = 3,
                                    values.transparency = 0.7)

  hidden <- plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                     which.data = "imputed",
                                     display.values = FALSE)

  expect_no_error(build.plot(shown@heatmap))
  expect_no_error(build.plot(hidden@heatmap))
})


test_that("the diagonal can be excluded and the scale tuned", {
  corr <- plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                   which.data = "imputed",
                                   exclude.diagonal = TRUE,
                                   correlation.scale.limits = c(NA, 1),
                                   plot.title = "**custom**",
                                   plot.subtitle = "a subtitle",
                                   clustering.method = "average")

  expect_s4_class(corr, "DEprot.correlation")
  expect_no_error(build.plot(corr@heatmap))
})


test_that("wrong inputs are rejected", {
  expect_error(plot.correlation.heatmap(DEprot.object = DEprot::sample.config))
  expect_error(plot.correlation.heatmap(DEprot.object = tb.dpo.imp, which.data = "not.a.type"))
  expect_error(plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                        correlation.method = "not.a.method"))
  expect_error(plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                        sample.subset = "not.a.sample"))
})



######################################    plot.upset    ##################################

test_that("the upset plot is generated and the contrasts can be subset", {
  up <- plot.upset(DEprot.analyses.object = tb.limma)
  subset.up <- plot.upset(DEprot.analyses.object = tb.limma, contrast.subset = 1)

  expect_s4_class(up, "DEprot.upset")
  expect_s4_class(subset.up, "DEprot.upset")
})


test_that("the sorting options are honoured", {
  for (int in c("cardinality", "degree")) {
    for (sets in list("descending", "ascending", FALSE)) {
      up <- suppressWarnings(plot.upset(DEprot.analyses.object = tb.limma,
                                        sort.intersections = int,
                                        sort.sets = sets))
      expect_s4_class(up, "DEprot.upset")
    }
  }
})


test_that("the aesthetics and the thresholds are honoured", {
  up <- plot.upset(DEprot.analyses.object = tb.limma,
                   title = "custom title",
                   subtitle = "custom subtitle",
                   intersection.bar.color = "steelblue",
                   setsize.bar.color = "indianred",
                   show.counts = FALSE,
                   min.size = 2,
                   height.ratio = 0.8,
                   width.ratio = 0.4)

  expect_s4_class(up, "DEprot.upset")

  uncorrected <- plot.upset(DEprot.analyses.object = tb.limma, use.uncorrected.pvalue = TRUE)
  expect_s4_class(uncorrected, "DEprot.upset")
})


test_that("wrong inputs are rejected", {
  expect_error(plot.upset(DEprot.analyses.object = tb.dpo.imp))
  expect_error(plot.upset(DEprot.analyses.object = DEprot::sample.config))
  expect_error(plot.upset(DEprot.analyses.object = tb.limma, contrast.subset = 100))
})



##############################    remove.undetected.proteins    ##########################

test_that("the proteins detected in too few samples are removed", {
  ## a real object is used: 'filter.proteins' rebuilds the boxplots stored in the object,
  ## which a synthetic object assembled by hand does not carry
  dpo <- tb.dpo.norm

  ## one protein quantified in a single sample only
  dpo@norm.counts[1, 2:ncol(dpo@norm.counts)] <- NA

  filtered <- suppressWarnings(suppressMessages(
    remove.undetected.proteins(DEprot.object = dpo, min.n.samples = 3, which.data = "normalized")))

  expect_s4_class(filtered, "DEprot")
  expect_false(rownames(dpo@norm.counts)[1] %in% rownames(filtered@norm.counts))
  expect_true(nrow(filtered@norm.counts) < nrow(dpo@norm.counts))
})


test_that("a threshold of one keeps every protein quantified at least once", {
  dpo <- tb.dpo.norm

  filtered <- suppressWarnings(suppressMessages(
    remove.undetected.proteins(DEprot.object = dpo, min.n.samples = 1, which.data = "normalized")))

  quantified <- sum(rowSums(!is.na(dpo@norm.counts)) >= 1)
  expect_equal(nrow(filtered@norm.counts), quantified)
})


test_that("every count type can be filtered and wrong inputs are rejected", {
  dpo <- tb.dpo.random

  for (w in c("raw", "normalized", "randomized")) {
    expect_s4_class(suppressWarnings(suppressMessages(
      remove.undetected.proteins(DEprot.object = dpo, min.n.samples = 2, which.data = w))),
      "DEprot")
  }

  expect_error(remove.undetected.proteins(DEprot.object = DEprot::sample.config))
  expect_error(remove.undetected.proteins(DEprot.object = dpo, which.data = "not.a.type"))
})
