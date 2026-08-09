## Real enrichment objects from the test toolbox
ora <- DEprot::test.toolbox$ora.results
gsea <- DEprot::test.toolbox$gsea.results


## Helper: minimal clusterProfiler-like ORA table, with controlled IDs and p-values.
## GeneRatio (count/50) and BgRatio (2*count/200) are set so that the FoldEnrichment
## is always exactly 2 and the geneset size exactly 2*count.
make_ora_df <- function(ids, p.adjust = 0.001, count = 5) {
  data.frame(ID = ids,
             Description = ids,
             GeneRatio = paste0(count, "/50"),
             BgRatio = paste0(count * 2, "/200"),
             pvalue = p.adjust / 10,
             p.adjust = p.adjust,
             qvalue = p.adjust,
             geneID = paste(paste0("prot.", 1:count), collapse = "/"),
             Count = count,
             stringsAsFactors = FALSE)
}


## Helper: minimal clusterProfiler-like GSEA table (no Count, no GeneRatio)
make_gsea_df <- function(ids, p.adjust = 0.001, NES = 2, setSize = 20, core = 6) {
  data.frame(ID = ids,
             Description = ids,
             setSize = setSize,
             enrichmentScore = NES / 2,
             NES = NES,
             pvalue = p.adjust / 10,
             p.adjust = p.adjust,
             qvalue = p.adjust,
             rank = 100,
             leading_edge = "tags=10%",
             core_enrichment = paste(paste0("prot.", 1:core), collapse = "/"),
             stringsAsFactors = FALSE)
}


## Helper: does the plot contain a text layer? (robust to layer ordering)
has_text_layer <- function(p) {
  any(vapply(p$layers, function(l){inherits(l$geom, "GeomText")}, logical(1)))
}


## Two dummy discoveries sharing part of the genesets
ora.A <- make_ora_df(ids = c("set.A", "set.B", "set.C", "set.D"),
                     p.adjust = c(0.001, 0.002, 0.003, 0.004))

ora.B <- make_ora_df(ids = c("set.C", "set.D", "set.E", "set.F"),
                     p.adjust = c(0.001, 0.002, 0.003, 0.004),
                     count = 8)

gsea.A <- make_gsea_df(ids = c("set.A", "set.B"))



# ---------------------------------------------------------------------
# .get.enrichment.table: harmonization of the metrics
# ---------------------------------------------------------------------

test_that(".get.enrichment.table returns the harmonized columns", {
  tb <- .get.enrichment.table(enrichment = ora.A, name = "test")

  expect_true(is.data.frame(tb))
  expect_true(all(c("discovery", "enrichment.type", "ID", "alias", "Count", "set.size",
                    "GeneRatio.numeric", "FoldEnrichment", "NES", "pvalue", "p.adjust") %in% colnames(tb)))
  expect_equal(unique(tb$discovery), "test")
  expect_equal(unique(tb$enrichment.type), "ORA")
})


test_that(".get.enrichment.table computes the FoldEnrichment when it is missing", {
  tb <- .get.enrichment.table(enrichment = ora.A, name = "test")

  expect_equal(unique(tb$FoldEnrichment), 2)
  expect_equal(unique(tb$GeneRatio.numeric), 0.1)
  expect_equal(unique(tb$set.size), 10)
})


test_that(".get.enrichment.table uses the leading edge as count for GSEA results", {
  tb <- .get.enrichment.table(enrichment = gsea.A, name = "test")

  expect_equal(unique(tb$enrichment.type), "GSEA")
  expect_equal(unique(tb$Count), 6)
  expect_equal(unique(tb$set.size), 20)
  expect_true(all(is.na(tb$FoldEnrichment)))
})


test_that(".get.enrichment.table accepts DEprot and clusterProfiler objects", {
  expect_true(is.data.frame(.get.enrichment.table(enrichment = ora, name = "DEprot")))
  expect_true(is.data.frame(.get.enrichment.table(enrichment = ora@enrichment.discovery, name = "clusterProfiler")))
  expect_true(is.data.frame(.get.enrichment.table(enrichment = gsea, name = "DEprot.gsea")))
})


test_that(".get.enrichment.table skips empty or unsupported enrichments", {
  expect_warning(expect_null(.get.enrichment.table(enrichment = NULL, name = "empty")))
  expect_warning(expect_null(.get.enrichment.table(enrichment = "not.an.enrichment", name = "wrong")))
  expect_warning(expect_null(.get.enrichment.table(enrichment = ora.A[,c("ID", "Count")], name = "incomplete")))
})


test_that(".collect.enrichments keeps the order of the input list", {
  tb <- .collect.enrichments(enrichment.list = list(second = ora.B, first = ora.A))

  expect_equal(levels(tb$discovery), c("second", "first"))
})


test_that(".collect.enrichments requires a list", {
  expect_error(.collect.enrichments(enrichment.list = ora.A))
  expect_error(.collect.enrichments(enrichment.list = list()))
})


test_that(".collect.enrichments assigns generic names to unnamed lists", {
  expect_warning(tb <- .collect.enrichments(enrichment.list = list(ora.A, ora.B)))
  expect_equal(levels(tb$discovery), c("enrichment.1", "enrichment.2"))
})



# ---------------------------------------------------------------------
# combine.enrichments: output
# ---------------------------------------------------------------------

test_that("combine.enrichments returns the results and the dotplot", {
  combined <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B))

  expect_type(combined, "list")
  expect_true(is.data.frame(combined$results))
  expect_s3_class(combined$dotplot, "ggplot")
  expect_equal(levels(combined$results$discovery), c("A", "B"))
})


test_that("combine.enrichments works with mixed input classes", {
  expect_no_error(combine.enrichments(enrichment.list = list(DEprot = ora,
                                                             clusterProfiler = ora@enrichment.discovery,
                                                             table = ora.A),
                                      padj.cutoff = 1,
                                      size.by = "GeneRatio"))
})


test_that("combine.enrichments displays the union of the top genesets", {
  combined <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                  dotplot.n = 2)

  ## top-2 of A (set.A, set.B) and top-2 of B (set.C, set.D)
  expect_equal(sort(levels(combined$dotplot$data$alias)), c("set.A", "set.B", "set.C", "set.D"))
})


test_that("combine.enrichments keeps only the requested genesets", {
  combined <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                  terms = c("set.C", "set.E"))

  expect_equal(sort(levels(combined$dotplot$data$alias)), c("set.C", "set.E"))
  expect_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                   terms = "set.not.existing"))
})


test_that("combine.enrichments writes the counts inside the dots only when asked", {
  expect_true(has_text_layer(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B))$dotplot))
  expect_false(has_text_layer(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                                  show.numbers = FALSE)$dotplot))
})



# ---------------------------------------------------------------------
# combine.enrichments: parameters
# ---------------------------------------------------------------------

test_that("combine.enrichments accepts all the ordering methods", {
  for (method in c("discovery", "significance", "clustering", "alphabetical")) {
    expect_no_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                        order.by = method))
  }
})


test_that("combine.enrichments accepts both the size metrics", {
  expect_no_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                      size.by = "FoldEnrichment"))
  expect_no_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                      size.by = "GeneRatio"))
})


test_that("combine.enrichments stops when the parameters are not supported", {
  expect_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B), size.by = "pvalue"))
  expect_error(combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B), order.by = "random"))
})


test_that("combine.enrichments warns when the FoldEnrichment is not available", {
  expect_warning(combine.enrichments(enrichment.list = list(ORA = ora.A, GSEA = gsea.A),
                                     size.by = "FoldEnrichment"))
  expect_no_warning(combine.enrichments(enrichment.list = list(ORA = ora.A, GSEA = gsea.A),
                                        size.by = "GeneRatio"))
})


test_that("combine.enrichments does not generate the dotplot when nothing is significant", {
  expect_warning(combined <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                                 padj.cutoff = 0))
  expect_null(combined$dotplot)
  expect_true(is.data.frame(combined$results))
})


test_that("combine.enrichments shows the non-significant genesets only when asked", {
  ## 'set.C' and 'set.D' are shared: in B they are above the cutoff used here
  full <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                              padj.cutoff = 0.0025)
  signif.only <- combine.enrichments(enrichment.list = list(A = ora.A, B = ora.B),
                                     padj.cutoff = 0.0025,
                                     show.non.significant = FALSE)

  expect_true(any(full$dotplot$data$significant == FALSE))
  expect_true(all(signif.only$dotplot$data$significant == TRUE))
})


test_that("combine.enrichments expands a time-course enrichment into one discovery per cluster", {
  skip_if_not(!is.null(methods::getClassDef("DEprot.timecourse.enrichment")),
              "the time-course module is not available")

  tc.results <- rbind(cbind(make_ora_df(ids = c("set.A", "set.B")), cluster = 1),
                      cbind(make_ora_df(ids = c("set.B", "set.C")), cluster = 2))

  tc <- methods::new(Class = "DEprot.timecourse.enrichment",
                     results = tc.results,
                     enrichment.per.cluster = NULL,
                     dotplot = NULL,
                     universe = paste0("prot.", 1:100),
                     parameters = list())

  combined <- combine.enrichments(enrichment.list = list(TC = tc))

  expect_equal(levels(combined$results$discovery), c("TC.cluster.1", "TC.cluster.2"))
})



# ---------------------------------------------------------------------
# divergent.enrichment: output
# ---------------------------------------------------------------------

test_that("divergent.enrichment returns the results and the plot", {
  divergent <- divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B))

  expect_type(divergent, "list")
  expect_true(is.data.frame(divergent$results))
  expect_s3_class(divergent$divergent.plot, "ggplot")
})


test_that("divergent.enrichment inverts the sign of the second enrichment", {
  divergent <- divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B))

  expect_true(all(divergent$results$signed.value[divergent$results$discovery == "up"] >= 0))
  expect_true(all(divergent$results$signed.value[divergent$results$discovery == "down"] <= 0))
})


test_that("divergent.enrichment keeps the top genesets of each side", {
  divergent <- divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B),
                                    top.n = 2)

  expect_equal(nrow(divergent$divergent.plot$data), 4)
  expect_equal(as.numeric(table(divergent$divergent.plot$data$discovery)), c(2, 2))
})


test_that("divergent.enrichment places the shared genesets on a single row", {
  divergent <- divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B),
                                    top.n = 4)

  ## 'set.C' and 'set.D' are found on both sides: 8 bars, 6 rows
  expect_equal(nrow(divergent$divergent.plot$data), 8)
  expect_equal(nlevels(divergent$divergent.plot$data$alias), 6)
})


test_that("divergent.enrichment works with mixed input classes", {
  expect_no_error(divergent.enrichment(enrichment.list = list(DEprot = ora,
                                                              table = ora.A),
                                       padj.cutoff = 1,
                                       value = "GeneRatio"))
})



# ---------------------------------------------------------------------
# divergent.enrichment: parameters
# ---------------------------------------------------------------------

test_that("divergent.enrichment requires exactly two enrichments", {
  expect_error(divergent.enrichment(enrichment.list = list(A = ora.A)))
  expect_error(divergent.enrichment(enrichment.list = list(A = ora.A, B = ora.B, C = ora.A)))
})


test_that("divergent.enrichment accepts all the metrics", {
  for (metric in c("FoldEnrichment", "GeneRatio", "Count", "padj")) {
    expect_no_error(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B),
                                         value = metric))
  }

  expect_no_error(divergent.enrichment(enrichment.list = list(up = gsea.A, down = gsea.A),
                                       value = "NES"))
})


test_that("divergent.enrichment stops when the parameters are not supported", {
  expect_error(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B), value = "qvalue"))
  expect_error(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B), top.by = "random"))
  expect_error(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B), terms = "set.not.existing"))
})


test_that("divergent.enrichment stops when the metric is not available", {
  expect_error(divergent.enrichment(enrichment.list = list(up = gsea.A, down = gsea.A),
                                    value = "FoldEnrichment"))
})


test_that("divergent.enrichment adds the counts only when asked", {
  expect_true(has_text_layer(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B))$divergent.plot))
  expect_false(has_text_layer(divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B),
                                                   add.counts = FALSE)$divergent.plot))
})


test_that("divergent.enrichment does not generate the plot when nothing is significant", {
  expect_warning(divergent <- divergent.enrichment(enrichment.list = list(up = ora.A, down = ora.B),
                                                   padj.cutoff = 0))
  expect_null(divergent$divergent.plot)
  expect_true(is.data.frame(divergent$results))
})
