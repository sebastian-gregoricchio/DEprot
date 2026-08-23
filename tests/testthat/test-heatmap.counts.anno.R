## ----------------------------------------------------------------------------------------
##  heatmap.counts.anno()
##  'ComplexHeatmap' and 'circlize' are optional dependencies: every test is skipped when
##  they are not installed. 'install.missing' is a string: "never" keeps the tests offline
##  and makes them fail loudly rather than trying to install a package during a check.
## ----------------------------------------------------------------------------------------

skip.no.heatmap <-
  function() {
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("circlize")
  }


## Object without protein annotation, used for the row-annotation errors.
## The objects of the toolbox have been generated before the introduction of the
## 'protein.info' slot, the removal makes the test independent from their content.
tb.limma.noinfo <- suppressMessages(add.protein.info(DEprot.object = tb.limma,
                                                     protein.info = NULL))


## Protein annotation built on the results of the first contrast: one discrete column
## (the differential status) and one numeric column (the ranking of the protein)
res.limma <- get.results(DEprot.analyses.object = tb.limma, contrast = 1)

info.limma <- data.frame(diff.status = res.limma$diff.status,
                         ranking = seq_len(nrow(res.limma)),
                         row.names = res.limma$prot.id,
                         stringsAsFactors = FALSE)

tb.limma.info <- add.protein.info(DEprot.object = tb.limma,
                                  protein.info = info.limma,
                                  overwrite = TRUE)


## Values of the metadata columns used in the tests
conditions <- unique(tb.limma@metadata$condition)



## ----------------------------------------------------------------------------------------
##  Object returned
## ----------------------------------------------------------------------------------------

test_that("heatmap.counts.anno returns a DEprot.counts.heatmap object", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            install.missing = "never")

  expect_s4_class(ht, "DEprot.counts.heatmap")
  expect_s4_class(ht@heatmap, "Heatmap")
  expect_true(nrow(ht@heatmap@matrix) > 1)
  expect_equal(ncol(ht@heatmap@matrix), ncol(tb.limma@imputed.counts))
  expect_s3_class(ht@row.cluster, "hclust")
  expect_s3_class(ht@column.cluster, "hclust")
})


test_that("the clusters are NULL when the clustering is disabled", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            clust.rows = FALSE,
                            clust.columns = FALSE,
                            install.missing = "never")

  expect_null(ht@row.cluster)
  expect_null(ht@column.cluster)
})


test_that("an object of the wrong class is rejected", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = data.frame(a = 1), install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Column annotation (metadata)
## ----------------------------------------------------------------------------------------

test_that("the metadata columns are attached as column annotation", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            column.annotation = c("condition", "replicate"),
                            install.missing = "never")

  expect_s4_class(ht@heatmap@top_annotation, "HeatmapAnnotation")
  expect_equal(length(ht@heatmap@top_annotation@anno_list), 2)
  expect_null(ht@heatmap@bottom_annotation)
})


test_that("the column annotation can be moved to the bottom", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            column.annotation = "condition",
                            column.annotation.side = "bottom",
                            install.missing = "never")

  expect_null(ht@heatmap@top_annotation)
  expect_s4_class(ht@heatmap@bottom_annotation, "HeatmapAnnotation")
})


test_that("a column absent from the metadata is rejected", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   column.annotation = "not.a.column",
                                   install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Row annotation (protein.info)
## ----------------------------------------------------------------------------------------

test_that("the row annotation requires a protein.info table", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma.noinfo,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   row.annotation = "diff.status",
                                   install.missing = "never"))
})


test_that("the protein.info columns are attached as row annotation", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma.info,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            row.annotation = c("diff.status", "ranking"),
                            install.missing = "never")

  expect_s4_class(ht@heatmap@left_annotation, "HeatmapAnnotation")
  expect_equal(length(ht@heatmap@left_annotation@anno_list), 2)
  expect_null(ht@heatmap@right_annotation)
})


test_that("the row annotation can be moved to the right", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma.info,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            row.annotation = "diff.status",
                            row.annotation.side = "right",
                            install.missing = "never")

  expect_null(ht@heatmap@left_annotation)
  expect_s4_class(ht@heatmap@right_annotation, "HeatmapAnnotation")
})


test_that("a column absent from the protein.info is rejected", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma.info,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   row.annotation = "not.a.column",
                                   install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Colors of the annotations (pheatmap syntax)
## ----------------------------------------------------------------------------------------

test_that("a named vector of colors is accepted for a discrete annotation", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            column.annotation = "condition",
                            annotation.colors = list(condition = stats::setNames(rep("grey", length(conditions)),
                                                                                 conditions)),
                            install.missing = "never")

  expect_s4_class(ht@heatmap@top_annotation, "HeatmapAnnotation")
})


test_that("a color missing for one of the values raises an error", {
  skip.no.heatmap()

  incomplete.colors <- stats::setNames(rep("grey", length(conditions) - 1), conditions[-1])

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   column.annotation = "condition",
                                   annotation.colors = list(condition = incomplete.colors),
                                   install.missing = "never"))
})


test_that("an unnamed vector of colors is rejected for a discrete annotation", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   column.annotation = "condition",
                                   annotation.colors = list(condition = rep("grey", length(conditions))),
                                   install.missing = "never"))
})


test_that("'annotation.colors' must be a named list", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   column.annotation = "condition",
                                   annotation.colors = c(condition = "grey"),
                                   install.missing = "never"))
})


test_that("colors not matching any annotation raise a warning", {
  skip.no.heatmap()

  expect_warning(heatmap.counts.anno(DEprot.object = tb.limma,
                                     contrast = 1,
                                     top.n = 10,
                                     use.uncorrected.pvalue = TRUE,
                                     column.annotation = "condition",
                                     annotation.colors = list(not.an.annotation = c(a = "grey")),
                                     install.missing = "never"))
})


test_that("a numeric annotation accepts both a gradient and a color function", {
  skip.no.heatmap()

  expect_no_error(heatmap.counts.anno(DEprot.object = tb.limma.info,
                                      contrast = 1,
                                      top.n = 10,
                                      use.uncorrected.pvalue = TRUE,
                                      row.annotation = "ranking",
                                      annotation.colors = list(ranking = c("white", "black")),
                                      install.missing = "never"))

  expect_no_error(heatmap.counts.anno(DEprot.object = tb.limma.info,
                                      contrast = 1,
                                      top.n = 10,
                                      use.uncorrected.pvalue = TRUE,
                                      row.annotation = "ranking",
                                      annotation.colors = list(ranking = circlize::colorRamp2(breaks = c(0, 100),
                                                                                              colors = c("white", "black"))),
                                      install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Averaged counts: the annotation must be constant within each group
## ----------------------------------------------------------------------------------------

test_that("the constant metadata columns annotate the averaged counts", {
  skip.no.heatmap()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            group.by.metadata.column = "combined.id",
                            column.annotation = "condition",
                            install.missing = "never")

  expect_equal(ncol(ht@heatmap@matrix), length(unique(tb.limma@metadata$combined.id)))
  expect_s4_class(ht@heatmap@top_annotation, "HeatmapAnnotation")
  expect_equal(length(ht@heatmap@top_annotation@anno_list), 1)
})


test_that("the columns varying within a group are dropped with a warning", {
  skip.no.heatmap()

  expect_warning(heatmap.counts.anno(DEprot.object = tb.limma,
                                     contrast = 1,
                                     top.n = 10,
                                     use.uncorrected.pvalue = TRUE,
                                     group.by.metadata.column = "combined.id",
                                     column.annotation = c("condition", "replicate"),
                                     install.missing = "never"))
})


test_that("a split on a column varying within a group is an error", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   group.by.metadata.column = "combined.id",
                                   column.split = "replicate",
                                   install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Splitting
## ----------------------------------------------------------------------------------------

test_that("the heatmap can be split by a metadata and a protein.info column", {
  skip.no.heatmap()

  expect_no_error(heatmap.counts.anno(DEprot.object = tb.limma.info,
                                      contrast = 1,
                                      top.n = 10,
                                      use.uncorrected.pvalue = TRUE,
                                      column.split = "condition",
                                      row.split = "diff.status",
                                      install.missing = "never"))
})


test_that("the dendrograms can be cut in a given number of blocks", {
  skip.no.heatmap()

  expect_no_error(heatmap.counts.anno(DEprot.object = tb.limma,
                                      contrast = 1,
                                      top.n = 10,
                                      use.uncorrected.pvalue = TRUE,
                                      column.split = 2,
                                      row.split = 2,
                                      install.missing = "never"))
})


test_that("a numeric split without clustering is ignored with a warning", {
  skip.no.heatmap()

  expect_warning(heatmap.counts.anno(DEprot.object = tb.limma,
                                     contrast = 1,
                                     top.n = 10,
                                     use.uncorrected.pvalue = TRUE,
                                     clust.rows = FALSE,
                                     row.split = 2,
                                     install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Data handling shared with heatmap.counts
## ----------------------------------------------------------------------------------------

test_that("the scaling by row and by column give different matrices", {
  skip.no.heatmap()

  ht.row <- heatmap.counts.anno(DEprot.object = tb.limma,
                                contrast = 1,
                                top.n = 10,
                                use.uncorrected.pvalue = TRUE,
                                scale = "row",
                                install.missing = "never")

  ht.column <- heatmap.counts.anno(DEprot.object = tb.limma,
                                   contrast = 1,
                                   top.n = 10,
                                   use.uncorrected.pvalue = TRUE,
                                   scale = "column",
                                   install.missing = "never")

  expect_false(identical(ht.row@heatmap@matrix, ht.column@heatmap@matrix))
})


test_that("the samples can be subset and the pattern removed from the protein names", {
  skip.no.heatmap()

  samples <- tb.limma@metadata$column.id[1:4]

  ht <- heatmap.counts.anno(DEprot.object = tb.limma,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            sample.subset = samples,
                            protein.names.pattern = "protein[.]",
                            install.missing = "never")

  expect_equal(colnames(ht@heatmap@matrix), samples)
  expect_false(any(grepl("protein[.]", rownames(ht@heatmap@matrix))))
})


test_that("counts that are not available are rejected", {
  skip.no.heatmap()

  expect_error(heatmap.counts.anno(DEprot.object = tb.dpo.raw,
                                   which.data = "imputed",
                                   install.missing = "never"))
})



## ----------------------------------------------------------------------------------------
##  Drawing
## ----------------------------------------------------------------------------------------

test_that("the heatmap is drawn without errors", {
  skip.no.heatmap()

  dir <- local.tmpdir()

  ht <- heatmap.counts.anno(DEprot.object = tb.limma.info,
                            contrast = 1,
                            top.n = 10,
                            use.uncorrected.pvalue = TRUE,
                            scale = "row",
                            column.annotation = c("condition", "replicate"),
                            row.annotation = c("diff.status", "ranking"),
                            row.split = "diff.status",
                            show.protein.names = TRUE,
                            title = "test",
                            install.missing = "never")

  grDevices::pdf(file.path(dir, "heatmap.pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(ComplexHeatmap::draw(ht@heatmap))
  expect_no_error(show(ht))
})
