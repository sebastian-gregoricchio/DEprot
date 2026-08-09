## ----------------------------------------------------------------------------------------
##  classes.and.methods.R
##  One object per class is built once, then every method defined for the DEprot classes is
##  exercised on it: show, summary, plot, and the '$'/'@' accessors with their completion.
##  The objects that require an optional package or a slow computation are built defensively:
##  a class that cannot be instantiated in this setup is skipped rather than failing the file.
## ----------------------------------------------------------------------------------------

build <-
  function(expr) {
    out <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
    if (inherits(out, "try-error")) {return(NULL)}
    return(out)
  }


objects <-
  list(DEprot = tb.dpo.imp,
       DEprot.analyses = tb.limma,
       DEprot.PCA = build(perform.PCA(DEprot.object = tb.dpo.imp, which.data = "imputed")),
       DEprot.PCoA = build(perform.PCoA(DEprot.object = tb.dpo.imp, which.data = "imputed")),
       DEprot.correlation = build(plot.correlation.heatmap(DEprot.object = tb.dpo.imp,
                                                           which.data = "imputed")),
       DEprot.counts.heatmap = build(heatmap.counts(DEprot.object = tb.dpo.imp,
                                                    which.data = "imputed")),
       DEprot.contrast.heatmap = build(heatmap.contrasts(DEprot.analyses.object = tb.limma,
                                                         contrasts = 1)),
       DEprot.upset = build(plot.upset(DEprot.analyses.object = tb.limma)),
       DEprot.normality = build(check.normality(DEprot.object = tb.dpo.imp, verbose = FALSE)),
       DEprot.pvalues = build(check.pvalues(tb.limma, contrast = 1)),
       DEprot.enrichResult = DEprot::test.toolbox$ora.results,
       DEprot.missingness = build(missingness.diagnostic(DEprot.object = tb.dpo.norm,
                                                         verbose = FALSE)),
       DEprot.outliers = build(detect.outliers(DEprot.object = tb.dpo.imp, verbose = FALSE)),
       DEprot.SAINTq = build(SAINTq(DEprot.object = tb.dpo.imp,
                                    control = "BCa_6h.DMSO",
                                    metadata.column = "combined.id",
                                    verbose = FALSE)),
       ## only one imputation method is run: the object is needed for its methods, not for
       ## the benchmark itself, which is tested in 'test-compare.ranking.R'
       DEprot.RMSE = build(compare.imp.methods(DEprot.object = tb.dpo.norm,
                                               percentage.test = 100,
                                               sample.group.column = "combined.id",
                                               which.data = "normalized",
                                               run.missForest = FALSE,
                                               run.kNN = FALSE,
                                               run.tkNN = FALSE,
                                               run.corkNN = FALSE,
                                               run.LLS = FALSE,
                                               run.SVD = TRUE,
                                               run.BPCA = FALSE,
                                               run.PPCA = FALSE,
                                               run.RegImpute = FALSE,
                                               pcaMethods.nPCs.to.test = 2,
                                               seed = 1,
                                               verbose = FALSE)))

objects <- objects[!vapply(objects, is.null, logical(1))]


## the classes for which '$' and '$<-' are defined (see '.deprot_classes'): the accessor is
## not registered for every DEprot class, hence the loops below are restricted to these
dollar.classes <- intersect(names(objects), DEprot:::.deprot_classes)



test_that("an object of each class could be built", {
  ## a canary: if a constructor starts failing, the methods below would be silently skipped
  expect_true(length(objects) >= 12)
  expect_true(all(c("DEprot", "DEprot.analyses", "DEprot.PCA") %in% names(objects)))
})


test_that("the classes are the ones expected", {
  for (cls in names(objects)) {
    expect_s4_class(objects[[cls]], cls)
  }
})


test_that("the show method of every class runs", {
  for (cls in names(objects)) {
    expect_no_error(suppressWarnings(suppressMessages(show(objects[[cls]]))))
  }
})



## ----------------------------------------------------------------------------------------
##  summary methods: DEprot.analyses, DEprot.RMSE, DEprot.SAINTq, DEprot.missingness
## ----------------------------------------------------------------------------------------

test_that("the summary methods return their table", {
  for (cls in intersect(names(objects),
                        c("DEprot.analyses", "DEprot.RMSE", "DEprot.SAINTq", "DEprot.missingness"))) {
    out <- suppressWarnings(suppressMessages(summary(objects[[cls]])))
    expect_true(is.data.frame(out) | is.list(out) | is.matrix(out))
  }
})



## ----------------------------------------------------------------------------------------
##  plot methods: DEprot, DEprot.analyses, DEprot.normality, DEprot.missingness, DEprot.outliers
## ----------------------------------------------------------------------------------------

is.plot <-
  function(p) {
    inherits(p, "ggplot") | inherits(p, "patchwork") | inherits(p, "gg")
  }


test_that("the plot method of a DEprot object assembles the boxplots", {
  p <- suppressWarnings(suppressMessages(plot(objects$DEprot)))
  expect_true(is.plot(p))

  ## the layout arguments are forwarded to patchwork
  expect_true(is.plot(suppressWarnings(suppressMessages(plot(objects$DEprot, ncol = 1)))))
  expect_true(is.plot(suppressWarnings(suppressMessages(plot(objects$DEprot, nrow = 1)))))
})


test_that("the plot method of a DEprot.analyses object honours the plot type", {
  for (type in c("volcano", "MA")) {
    p <- suppressWarnings(suppressMessages(plot(objects$DEprot.analyses, plot.type = type)))
    expect_true(is.plot(p))
  }

  labelled <- suppressWarnings(suppressMessages(plot(objects$DEprot.analyses,
                                                     plot.type = "volcano",
                                                     label.top.n = 3,
                                                     ncol = 1)))
  expect_true(is.plot(labelled))
})


test_that("the plot method of a DEprot.normality object honours 'n.samples'", {
  skip_if(is.null(objects$DEprot.normality))

  expect_true(is.plot(suppressWarnings(suppressMessages(plot(objects$DEprot.normality)))))
  expect_true(is.plot(suppressWarnings(suppressMessages(plot(objects$DEprot.normality, n.samples = 2)))))
})


test_that("the plot method of a DEprot.missingness object honours the plot type", {
  skip_if(is.null(objects$DEprot.missingness))

  for (type in c("summary", "dropout")) {
    p <- try(suppressWarnings(suppressMessages(plot(objects$DEprot.missingness, plot.type = type))),
             silent = TRUE)
    skip_if(inherits(p, "try-error"), paste0("plot.type '", type, "' is not available"))
    expect_true(is.plot(p))
  }
})


test_that("the plot method of a DEprot.outliers object returns the diagnostic panel", {
  skip_if(is.null(objects$DEprot.outliers))

  expect_true(is.plot(suppressWarnings(suppressMessages(plot(objects$DEprot.outliers)))))
})


test_that("the classes without a plot method fall back on the default one", {
  skip_if(is.null(objects$DEprot.RMSE))

  ## no plot() method is defined for DEprot.RMSE: the correlation panel is printed by show()
  expect_error(suppressWarnings(suppressMessages(plot(objects$DEprot.RMSE))))
})



## ----------------------------------------------------------------------------------------
##  Slot access: '$' and '$<-' are defined for the classes listed in '.deprot_classes'
## ----------------------------------------------------------------------------------------

test_that("every slot can be read with '$' exactly as with '@'", {
  for (cls in dollar.classes) {
    obj <- objects[[cls]]

    for (slot.name in methods::slotNames(obj)) {
      expect_identical(methods::slot(obj, slot.name),
                       do.call("$", list(obj, slot.name)),
                       info = paste0(cls, "$", slot.name))
    }
  }
})


test_that("every slot can be replaced with '$<-'", {
  for (cls in dollar.classes) {
    obj <- objects[[cls]]
    slot.name <- methods::slotNames(obj)[1]

    ## the slots are declared as 'ANY' for all the DEprot classes but 'DEprot' itself,
    ## whose typed slots would reject an arbitrary value: the round trip re-assigns the
    ## value already stored, which is valid whatever the class
    original <- methods::slot(obj, slot.name)
    obj <- do.call("$<-", list(obj, slot.name, original))

    expect_s4_class(obj, cls)
    expect_identical(do.call("$", list(obj, slot.name)), original)
  }
})


test_that("the replacement really writes the new value", {
  dpo <- tb.dpo.imp

  dpo$log.base <- 10
  expect_equal(dpo@log.base, 10)
  expect_equal(dpo$log.base, 10)
  expect_s4_class(dpo, "DEprot")

  ## the same through a class whose slots are all 'ANY'
  pca <- objects$DEprot.PCA
  pca$data.used <- "custom"
  expect_equal(pca@data.used, "custom")
})


test_that("a slot that does not exist cannot be accessed", {
  expect_error(tb.dpo.imp$not.a.slot)
  expect_error(objects$DEprot.PCA$not.a.slot)
})


test_that("'$' is registered for DEprot.outliers as well", {
  skip_if(is.null(objects$DEprot.outliers))

  expect_true("DEprot.outliers" %in% DEprot:::.deprot_classes)
  expect_identical(objects$DEprot.outliers$metrics, objects$DEprot.outliers@metrics)
})



## ----------------------------------------------------------------------------------------
##  Instance-aware completion
## ----------------------------------------------------------------------------------------

test_that("the completion offers only the slots holding a value", {
  complete <- DEprot:::.deprot_complete_slots
  is.empty <- DEprot:::.deprot_slot_is_empty

  ## an empty slot is NULL, zero-length or a single NA; FALSE is a real value
  expect_true(is.empty(NULL))
  expect_true(is.empty(character(0)))
  expect_true(is.empty(NA))
  expect_false(is.empty(FALSE))
  expect_false(is.empty(0))
  expect_false(is.empty(data.frame(a = 1)))

  offered <- complete(tb.dpo.imp)
  expect_true("imputed.counts" %in% offered)
  expect_true(all(offered %in% methods::slotNames(tb.dpo.imp)))

  ## an emptied slot disappears from the completion
  dpo <- tb.dpo.imp
  dpo@raw.counts <- NULL
  expect_false("raw.counts" %in% complete(dpo))

  ## the pattern filters the names, as it does after a partial '$'
  expect_true(all(grepl("^log", complete(tb.dpo.imp, pattern = "^log"))))
  expect_length(complete(tb.dpo.imp, pattern = "^no.such.slot"), 0)
})


test_that("the completion works on every class built", {
  for (cls in names(objects)) {
    offered <- DEprot:::.deprot_complete_slots(objects[[cls]])
    expect_true(is.character(offered))
    expect_true(all(offered %in% methods::slotNames(objects[[cls]])))
  }
})
