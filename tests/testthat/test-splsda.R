## ----------------------------------------------------------------------------------------
##  sPLS-DA module
##  The shared model is fitted with a fixed 'keepX' and a minimal resampling, so that the
##  suite stays fast; the tuning is exercised once, apart.
## ----------------------------------------------------------------------------------------

dpo <- DEprot::test.toolbox$dpo.imp

splsda <- perform.sPLSDA(DEprot.object = dpo,
                         group.column = "condition",
                         which.data = "imputed",
                         ncomp = 3,
                         keepX = c(10, 10, 5),
                         validate = TRUE,
                         folds = 3,
                         nrepeat = 2,
                         seed = 1234)



######################################    ERRORS    ######################################

test_that("errors are returned if the object is not of class DEprot or DEprot.analyses", {
  expect_error(perform.sPLSDA(DEprot.object = DEprot::sample.config, group.column = "condition"))
  expect_error(perform.sPLSDA(DEprot.object = "not an object", group.column = "condition"))
})


test_that("counts that are not available cannot be used", {
  bare <- dpo
  bare@imputed.counts <- NULL

  expect_error(perform.sPLSDA(DEprot.object = bare, group.column = "condition", which.data = "imputed"))
  expect_error(perform.sPLSDA(DEprot.object = dpo, group.column = "condition", which.data = "wrong data"))
})


test_that("a group column that does not exist is rejected", {
  expect_error(perform.sPLSDA(DEprot.object = dpo, group.column = "not.a.column"))
})


test_that("a sample subset that does not exist is rejected", {
  expect_error(perform.sPLSDA(DEprot.object = dpo,
                              group.column = "condition",
                              sample.subset = c("not.a.sample")))
})


test_that("a reference group that is not one of the classes is rejected", {
  expect_error(perform.sPLSDA(DEprot.object = dpo,
                              group.column = "condition",
                              keepX = 5,
                              reference.group = "not.a.class"))
})


test_that("counts carrying missing values are refused", {
  ## both the tuning and the validation resample the samples: a missing value would be
  ## handled differently at every fold
  expect_error(perform.sPLSDA(DEprot.object = dpo,
                              group.column = "condition",
                              which.data = "normalized",
                              keepX = 5))
})


test_that("a keepX of the wrong length is rejected", {
  expect_error(perform.sPLSDA(DEprot.object = dpo,
                              group.column = "condition",
                              ncomp = 3,
                              keepX = c(10, 10)))
})


test_that("a single class cannot be discriminated", {
  single <- dpo
  single@metadata$only.one <- "same.group"

  expect_error(perform.sPLSDA(DEprot.object = single, group.column = "only.one", keepX = 5))
})



##########################################################################################
##  Object built by perform.sPLSDA
##########################################################################################

test_that("the sPLSDA object carries the scores, the loadings and the metadata", {
  expect_s4_class(splsda, "DEprot.sPLSDA")
  expect_equal(splsda@data.used, "imputed")
  expect_equal(splsda@group.column, "condition")
  expect_equal(nrow(splsda@components), ncol(dpo@imputed.counts))
  expect_true(all(c("comp1", "comp2", "comp3") %in% colnames(splsda@components)))
  expect_s3_class(splsda@loadings, "data.frame")
  expect_s3_class(splsda@importance, "data.frame")
})


test_that("ncomp is raised to 3 when a lower value is requested", {
  expect_message(small <- perform.sPLSDA(DEprot.object = dpo,
                                         group.column = "condition",
                                         ncomp = 2,
                                         keepX = 5,
                                         validate = FALSE))

  expect_equal(small@parameters$ncomp, 3)
  expect_true("comp3" %in% colnames(small@components))
})


test_that("keepX is respected and the loadings of the discarded proteins are exactly zero", {
  selected <- splsda@selected.proteins

  expect_equal(length(selected), 3)
  expect_equal(vapply(selected, length, integer(1), USE.NAMES = FALSE), c(10, 10, 5))

  comp1 <- splsda@loadings[splsda@loadings$component == 1, ]
  expect_equal(sum(comp1$loading != 0), 10)
  expect_true(all(comp1$loading[!comp1$selected] == 0))
})


test_that("a keepX larger than the number of proteins is capped", {
  expect_warning(capped <- perform.sPLSDA(DEprot.object = dpo,
                                          group.column = "condition",
                                          keepX = nrow(dpo@imputed.counts) + 100,
                                          validate = FALSE))

  expect_true(all(capped@parameters$keepX <= nrow(capped@counts.used)))
})


test_that("the reference group sits on the positive side of every component", {
  reference <- splsda@parameters$reference.group
  scores <- splsda@components

  for (k in seq_len(splsda@parameters$ncomp)) {
    expect_gte(mean(scores[[paste0("comp", k)]][scores$condition == reference]), 0)
  }
})


test_that("the orientation of the mixOmics object is left untouched", {
  ## the sign correction is applied to the derived tables only: flipping the loadings of the
  ## fitted object without its deflation coefficients would break mixOmics::predict
  expect_equal(abs(as.numeric(splsda@splsda$variates$X[, 1])),
               abs(splsda@components$comp1),
               tolerance = 1e-8)
})


test_that("the companion non-sparse model covers every protein", {
  expect_equal(nrow(splsda@full.loadings),
               nrow(splsda@counts.used) * splsda@parameters$ncomp)

  comp1 <- splsda@full.loadings[splsda@full.loadings$component == 1, ]
  expect_equal(sum(comp1$loading == 0), 0)
})


test_that("the explained variance is bounded and cumulative", {
  expect_true(all(c("Proportion.of.Variance", "Cumulative.Proportion") %in% colnames(splsda@importance)))
  expect_true(all(splsda@importance$Proportion.of.Variance >= 0 &
                    splsda@importance$Proportion.of.Variance <= 1))
  expect_true(!is.unsorted(splsda@importance$Cumulative.Proportion))
})


test_that("the selected proteins are annotated with the class in which they are the highest", {
  selected <- splsda@loadings[splsda@loadings$selected %in% TRUE, ]

  expect_false(any(is.na(selected$contrib.group)))
  expect_true(all(selected$contrib.group %in% splsda@parameters$groups))
})


test_that("the validation fills the performance slot", {
  expect_false(is.null(splsda@performance))
  expect_s3_class(splsda@performance$error.rate, "data.frame")
  expect_true(all(c("Overall error rate", "Balanced error rate") %in% splsda@performance$error.rate$metric))
  expect_true(all(splsda@performance$error.rate$error.rate >= 0 &
                    splsda@performance$error.rate$error.rate <= 1))
})


test_that("the validation can be skipped", {
  quick <- perform.sPLSDA(DEprot.object = dpo,
                          group.column = "condition",
                          keepX = 5,
                          validate = FALSE)

  expect_null(quick@performance)
  expect_null(quick@performance.plot)
  expect_error(plot.sPLSDA.performance(DEprot.sPLSDA.object = quick))
})


test_that("no tuning is stored when keepX is provided", {
  expect_null(splsda@tuning)
  expect_null(splsda@tuning.plot)
  expect_error(plot.sPLSDA.tuning(DEprot.sPLSDA.object = splsda))
})


test_that("the number of folds is capped to the smallest class", {
  expect_warning(capped <- perform.sPLSDA(DEprot.object = dpo,
                                          group.column = "condition",
                                          keepX = 5,
                                          folds = 50,
                                          nrepeat = 1))

  expect_lte(capped@parameters$folds, min(table(capped@groups)))
})


test_that("the tuning explores the grid and returns one value per component", {
  skip_on_cran()

  tuned <- perform.sPLSDA(DEprot.object = dpo,
                          group.column = "condition",
                          ncomp = 3,
                          keepX = NULL,
                          test.keepX = c(5, 10, 20),
                          validate = FALSE,
                          folds = 3,
                          nrepeat = 1,
                          seed = 1234)

  expect_false(is.null(tuned@tuning))
  expect_equal(length(tuned@tuning$choice.keepX), 3)
  expect_true(all(tuned@tuning$choice.keepX %in% c(5, 10, 20)))
  expect_s3_class(tuned@tuning$error.rate, "data.frame")
  expect_true(inherits(tuned@tuning.plot, "ggplot"))
})



##########################################################################################
##  Getters
##########################################################################################

test_that("get.sPLSDA.results returns the selected proteins by default", {
  res <- get.sPLSDA.results(DEprot.sPLSDA.object = splsda, component = 1)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  expect_true(all(res$component == 1))
  ## the table is sorted by decreasing absolute loading
  expect_false(is.unsorted(rev(res$abs.loading)))
})


test_that("get.sPLSDA.results can return the complete loadings", {
  res <- get.sPLSDA.results(DEprot.sPLSDA.object = splsda, component = 1, selected.only = FALSE)

  expect_equal(nrow(res), nrow(splsda@counts.used))
})


test_that("a component that does not exist is rejected by the getters", {
  expect_error(get.sPLSDA.results(DEprot.sPLSDA.object = splsda, component = 99))
  expect_error(get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, component = 99))
  expect_error(get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda, component = 99))
})


test_that("get.sPLSDA.proteins combines the components as requested", {
  union.proteins <- get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, mode = "union")
  list.proteins <- get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, mode = "list")

  expect_type(union.proteins, "character")
  expect_type(list.proteins, "list")
  expect_equal(length(list.proteins), 3)
  expect_setequal(union.proteins, unique(unlist(list.proteins, use.names = FALSE)))
  expect_true(all(union.proteins %in% rownames(splsda@counts.used)))
})


test_that("get.sPLSDA.proteins filters by direction and by rank", {
  positive <- get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, component = 1, direction = "positive")
  negative <- get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, component = 1, direction = "negative")
  top.three <- get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, component = 1, top.n = 3)

  expect_equal(length(positive) + length(negative), 10)
  expect_length(top.three, 3)
  expect_error(get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, direction = "sideways"))
  expect_error(get.sPLSDA.proteins(DEprot.sPLSDA.object = splsda, mode = "whatever"))
})


test_that("the ranking of the non-sparse model has no ties at zero", {
  ranking <- get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda, component = 1)

  expect_type(ranking, "double")
  expect_equal(length(ranking), nrow(splsda@counts.used))
  expect_false(is.unsorted(rev(ranking)))
  expect_equal(sum(ranking == 0), 0)
})


test_that("the sparse ranking is available but warns about the ties", {
  expect_warning(sparse.ranking <- get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda,
                                                      component = 1,
                                                      metric = "sparse"))

  expect_true(sum(sparse.ranking == 0) > 0)
  expect_error(get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda, metric = "wrong.metric"))
})


test_that("the VIP ranking is unsigned", {
  vip.ranking <- get.sPLSDA.ranking(DEprot.sPLSDA.object = splsda, component = 1, metric = "vip")

  expect_true(all(vip.ranking >= 0))
})



##########################################################################################
##  Plots
##########################################################################################

test_that("the scatters are generated", {
  expect_true(inherits(plot.sPLSDA.scatter(DEprot.sPLSDA.object = splsda), "ggplot"))
  expect_true(inherits(plot.sPLSDA.scatter(DEprot.sPLSDA.object = splsda,
                                           comp.x = 1,
                                           comp.y = 3,
                                           shape.column = "replicate",
                                           ellipse = FALSE), "ggplot"))
  expect_true(inherits(splsda@scatter.plot.123, "patchwork"))
})


test_that("axes and aesthetic columns that do not exist are rejected", {
  expect_error(plot.sPLSDA.scatter(DEprot.sPLSDA.object = splsda, comp.x = 99))
  expect_error(plot.sPLSDA.scatter(DEprot.sPLSDA.object = splsda, color.column = "not.a.column"))
  expect_error(plot.sPLSDA.scatter(DEprot.sPLSDA.object = splsda, shape.column = "not.a.column"))
})


test_that("the other plots are generated", {
  expect_true(inherits(plot.sPLSDA.cumulative(DEprot.sPLSDA.object = splsda), "ggplot"))
  expect_true(inherits(plot.sPLSDA.biplot(DEprot.sPLSDA.object = splsda, n.loadings = 5), "ggplot"))
  expect_true(inherits(plot.sPLSDA.loadings(DEprot.sPLSDA.object = splsda, component = 1), "ggplot"))
  expect_true(inherits(plot.sPLSDA.performance(DEprot.sPLSDA.object = splsda), "ggplot"))
  expect_true(inherits(splsda@cumulative.plot, "ggplot"))
})


test_that("the stability plot is generated when the folds were recorded", {
  skip_if(is.null(splsda@performance$stability))

  expect_true(inherits(plot.sPLSDA.stability(DEprot.sPLSDA.object = splsda, component = 1), "ggplot"))
})


test_that("the AUC plot is generated when the AUC was computed", {
  skip_if(is.null(splsda@performance$auc))

  expect_true(inherits(plot.sPLSDA.auroc(DEprot.sPLSDA.object = splsda), "ggplot"))
})


test_that("the plotting functions reject the wrong class", {
  expect_error(plot.sPLSDA.scatter(DEprot.sPLSDA.object = dpo))
  expect_error(plot.sPLSDA.loadings(DEprot.sPLSDA.object = "not an object"))
})



##########################################################################################
##  Methods
##########################################################################################

test_that("the show and summary methods behave", {
  expect_output(show(splsda))

  smry <- summary(splsda)

  expect_s3_class(smry, "data.frame")
  expect_equal(nrow(smry), 3)
  expect_true(all(c("component", "keepX", "n.selected") %in% colnames(smry)))
  expect_equal(smry$n.selected, c(10, 10, 5))
})


test_that("the plot method dispatches on 'what'", {
  expect_true(inherits(plot(splsda, what = "cumulative"), "ggplot"))
  expect_true(inherits(plot(splsda, what = "performance"), "ggplot"))
  expect_error(plot(splsda, what = "not.a.plot"))
  expect_error(plot(splsda, what = "tuning"))
})


test_that("the slots are reachable with the dollar operator", {
  ## requires 'DEprot.sPLSDA' to be listed in '.deprot_classes'
  expect_equal(splsda$group.column, splsda@group.column)
  expect_equal(nrow(splsda$components), nrow(splsda@components))
})



##########################################################################################
##  Enrichment
##########################################################################################

test_that("the ORA runs on the selected proteins", {
  ora <- sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
                           TERM2GENE = DEprot::test.toolbox$geneset,
                           enrichment.type = "ORA",
                           component = 1,
                           pvalueCutoff = 1,
                           qvalueCutoff = 1)

  expect_s4_class(ora, "DEprot.enrichResult")
  expect_equal(ora@parameters$enrichment.type, "ORA")
  expect_equal(ora@parameters$source, "sPLS-DA")
  ## the universe is the list of proteins that entered the model
  expect_setequal(ora@parameters$universe, rownames(splsda@counts.used))
})


test_that("the GSEA runs on the complete ranking", {
  gsea <- sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
                            TERM2GENE = DEprot::test.toolbox$geneset,
                            enrichment.type = "GSEA",
                            component = 1,
                            pvalueCutoff = 1)

  expect_s4_class(gsea, "DEprot.enrichResult")
  expect_equal(gsea@parameters$enrichment.type, "GSEA")
  expect_equal(gsea@parameters$gsea.rank.method, "loading")
  ## the labels of the two sides are what the downstream plots read
  expect_equal(gsea@parameters$contrast$var.1, splsda@parameters$reference.group)
})


test_that("the enrichment rejects the wrong inputs", {
  expect_error(sPLSDA.enrichment(DEprot.sPLSDA.object = dpo,
                                 TERM2GENE = DEprot::test.toolbox$geneset,
                                 enrichment.type = "ORA"))

  expect_error(sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
                                 TERM2GENE = DEprot::test.toolbox$geneset,
                                 enrichment.type = "wrong.type"))

  expect_error(sPLSDA.enrichment(DEprot.sPLSDA.object = splsda,
                                 TERM2GENE = "not a data.frame",
                                 enrichment.type = "ORA"))
})
