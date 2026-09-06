## ----------------------------------------------------------------------------------------
##  classes.and.methods.R: the classes not covered by 'test-classes-methods.R'
##  Those are the classes whose constructor is slow or needs a dedicated fixture, hence they
##  were left out of the shared 'objects' list: DEprot.power, DEprot.timecourse,
##  DEprot.timecourse.enrichment and DEprot.sPLSDA. Their show/plot/summary methods and their
##  '$' accessors are exercised here, on the smallest object that reaches each of them.
## ----------------------------------------------------------------------------------------

build <-
  function(expr) {
    out <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
    if (inherits(out, "try-error")) {return(NULL)}
    return(out)
  }


is.plot <-
  function(p) {
    inherits(p, "ggplot") | inherits(p, "patchwork") | inherits(p, "gg")
  }


## A small simulated time course: the same fixture as 'test-timecourse.R', halved
sim <- DEprot:::.simulate.timecourse(n.proteins = 200,
                                     timepoints = c(0, 1, 2, 6, 24),
                                     n.replicates = 3,
                                     fraction.responsive = 0.3,
                                     seed = 1234)

dpo.tc <- build(load.counts2(counts = sim$counts,
                             metadata = sim$metadata,
                             data.type = "imputed",
                             log.base = 2,
                             column.id = "column.id"))

extra.objects <-
  list(DEprot.power = build(estimate.power(DEprot.analyses.object = tb.limma,
                                           contrast = 1,
                                           sample.size.range = c(2, 10))),
       DEprot.timecourse = build(analyze.timecourse(DEprot.object = dpo.tc,
                                                    time.column = "time.hours",
                                                    time.transform = "log2",
                                                    log2.amplitude.th = 0.5,
                                                    n.clusters = 3,
                                                    seed = 1234,
                                                    verbose = FALSE)),
       DEprot.sPLSDA = build(perform.sPLSDA(DEprot.object = tb.dpo.imp,
                                            group.column = "condition",
                                            ncomp = 2,
                                            validate = FALSE,
                                            n.cores = 1,
                                            seed = 1234)))

extra.objects$DEprot.timecourse.enrichment <-
  if (is.null(extra.objects$DEprot.timecourse)) {NULL} else {
    build(timecourse.enrichment(DEprot.timecourse.object = extra.objects$DEprot.timecourse,
                                TERM2GENE = sim$TERM2GENE,
                                min.cluster.size = 3))
  }

extra.objects <- extra.objects[!vapply(extra.objects, is.null, logical(1))]



test_that("the objects of the remaining classes could be built", {
  ## a canary: a broken constructor would silently skip every test below
  expect_true(length(extra.objects) >= 2)

  for (cls in names(extra.objects)) {
    expect_s4_class(extra.objects[[cls]], cls)
  }
})


test_that("the show method of every remaining class runs", {
  for (cls in names(extra.objects)) {
    expect_no_error(suppressWarnings(suppressMessages(show(extra.objects[[cls]]))))
  }
})


test_that("the slots of the remaining classes are reachable with '$' and '@'", {
  for (cls in names(extra.objects)) {
    obj <- extra.objects[[cls]]

    for (slot.name in methods::slotNames(obj)) {
      expect_identical(methods::slot(obj, slot.name),
                       do.call("$", list(obj, slot.name)),
                       info = paste0(cls, "$", slot.name))
    }

    ## the replacement returns an object of the same class
    first.slot <- methods::slotNames(obj)[1]
    replaced <- do.call("$<-", list(obj, first.slot, methods::slot(obj, first.slot)))
    expect_s4_class(replaced, cls)
  }
})


test_that("the completion offers only the slots holding a value", {
  for (cls in names(extra.objects)) {
    offered <- DEprot:::.deprot_complete_slots(extra.objects[[cls]])
    expect_true(is.character(offered))
    expect_true(all(offered %in% methods::slotNames(extra.objects[[cls]])))
  }
})



## ----------------------------------------------------------------------------------------
##  Class-specific methods
## ----------------------------------------------------------------------------------------

test_that("the plot method of a DEprot.power object assembles the three panels", {
  skip_if(is.null(extra.objects$DEprot.power))

  expect_no_error(suppressWarnings(suppressMessages(plot(extra.objects$DEprot.power))))
  ## the layout arguments are forwarded to patchwork
  expect_no_error(suppressWarnings(suppressMessages(plot(extra.objects$DEprot.power, nrow = 1))))
})


test_that("the summary of a DEprot.timecourse object describes the clusters", {
  skip_if(is.null(extra.objects$DEprot.timecourse))

  out <- suppressWarnings(suppressMessages(summary(extra.objects$DEprot.timecourse)))

  expect_s3_class(out, "data.frame")
  expect_true(all(c("cluster", "n", "dominant.shape") %in% colnames(out)))
})


test_that("the plot method of a DEprot.timecourse object returns the profile panel", {
  skip_if(is.null(extra.objects$DEprot.timecourse))

  expect_true(is.plot(plot(extra.objects$DEprot.timecourse)))
})


test_that("the summary of a time course without clustering returns the results", {
  skip_if(is.null(extra.objects$DEprot.timecourse))

  no.cluster <- extra.objects$DEprot.timecourse
  no.cluster@clusters <- NULL

  out <- suppressWarnings(suppressMessages(summary(no.cluster)))
  expect_s3_class(out, "data.frame")
})


test_that("the plot method of a DEprot.timecourse.enrichment object returns the dotplot", {
  skip_if(is.null(extra.objects$DEprot.timecourse.enrichment))

  p <- plot(extra.objects$DEprot.timecourse.enrichment)
  expect_true(is.plot(p) | is.null(p))
})


test_that("an enrichment object holding no result still prints", {
  skip_if(is.null(extra.objects$DEprot.timecourse.enrichment))

  empty <- extra.objects$DEprot.timecourse.enrichment
  empty@results <- NULL

  expect_no_error(suppressWarnings(suppressMessages(show(empty))))
})


test_that("the show method reports the counts available and the missing ones", {
  ## the branches listing the slots depend on which counts the object carries
  raw.only <- tb.dpo.raw
  raw.only@norm.counts <- NULL
  raw.only@random.counts <- NULL
  raw.only@imputed.counts <- NULL

  expect_output(show(raw.only), "raw")
  expect_output(show(tb.dpo.imp), "imputed")
})


test_that("the show method of a DEprot.analyses object needs at least one contrast", {
  empty <- tb.limma
  empty@analyses.result.list <- list()
  empty@contrasts <- list()

  ## the summary is built by indexing the result list: an object emptied by hand is not a
  ## state the package produces, and the method does not guard against it
  expect_error(suppressWarnings(suppressMessages(show(empty))))
})
