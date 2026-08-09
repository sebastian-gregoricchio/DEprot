## Small simulated dataset, shared by all the tests of this file.
## 400 proteins are enough to exercise every code path and keep the suite fast.
sim <- DEprot:::.simulate.timecourse(n.proteins = 400,
                                     timepoints = c(0, 1, 2, 6, 24),
                                     n.replicates = 3,
                                     fraction.responsive = 0.3,
                                     seed = 1234)

dpo <- load.counts2(counts = sim$counts,
                    metadata = sim$metadata,
                    data.type = "imputed",
                    log.base = 2,
                    column.id = "column.id")

tc <- analyze.timecourse(DEprot.object = dpo,
                         time.column = "time.hours",
                         time.transform = "log2",
                         log2.amplitude.th = 0.5,
                         n.clusters = 5,
                         seed = 1234,
                         verbose = FALSE)


## Two-group version, used for the interaction model
sim.groups <- DEprot:::.simulate.timecourse(n.proteins = 300,
                                            timepoints = c(0, 1, 2, 6, 24),
                                            n.replicates = 3,
                                            groups = c("treated", "control"),
                                            fraction.responsive = 0.25,
                                            seed = 42)

dpo.groups <- load.counts2(counts = sim.groups$counts,
                           metadata = sim.groups$metadata,
                           data.type = "imputed",
                           log.base = 2,
                           column.id = "column.id")

tc.groups <- analyze.timecourse(DEprot.object = dpo.groups,
                                time.column = "time.hours",
                                group.column = "group",
                                time.transform = "log2",
                                log2.amplitude.th = 0.5,
                                n.clusters = 4,
                                seed = 1234,
                                verbose = FALSE)



##########################################################################################
####                                   SIMULATOR                                      ####
##########################################################################################

test_that("the simulator returns the expected structure", {
  expect_named(sim, c("counts", "metadata", "truth", "TERM2GENE"))
  expect_equal(nrow(sim$counts), 400)
  expect_equal(ncol(sim$counts), 15)
  expect_equal(colnames(sim$counts), sim$metadata$column.id)
})


test_that("the simulated sample IDs are unique, also with several groups", {
  expect_false(any(duplicated(sim$metadata$column.id)))
  expect_false(any(duplicated(sim.groups$metadata$column.id)))
  expect_equal(nrow(sim.groups$metadata), 30)
})


test_that("only the first group responds to the time", {
  ## a responsive protein must be flat in the control arm
  responsive <- sim.groups$truth$prot.id[sim.groups$truth$archetype != "flat"][1]
  ctrl <- sim.groups$metadata$column.id[sim.groups$metadata$group == "control"]

  expect_lt(diff(range(sim.groups$counts[responsive, ctrl])), 2)
})



##########################################################################################
####                                     ERRORS                                       ####
##########################################################################################

test_that("errors are returned if the object is not of class DEprot", {
  expect_error(analyze.timecourse(DEprot.object = DEprot::sample.config,
                                  time.column = "time.hours"))
})


test_that("errors are returned if the time column is absent or not numeric", {
  expect_error(analyze.timecourse(DEprot.object = dpo, time.column = "not.a.column"))
  expect_error(analyze.timecourse(DEprot.object = dpo, time.column = "replicate"))
})


test_that("errors are returned if the group column is absent", {
  expect_error(analyze.timecourse(DEprot.object = dpo,
                                  time.column = "time.hours",
                                  group.column = "not.a.column"))
})


test_that("errors are returned with less than three timepoints", {
  sim.short <- DEprot:::.simulate.timecourse(n.proteins = 100,
                                             timepoints = c(0, 6),
                                             n.replicates = 3,
                                             seed = 1)

  dpo.short <- load.counts2(counts = sim.short$counts,
                            metadata = sim.short$metadata,
                            data.type = "imputed",
                            log.base = 2,
                            column.id = "column.id")

  expect_error(analyze.timecourse(DEprot.object = dpo.short,
                                  time.column = "time.hours",
                                  verbose = FALSE))
})


test_that("errors are returned by the getters and plots with a wrong object class", {
  expect_error(rank.timecourse(DEprot.timecourse.object = dpo))
  expect_error(get.timecourse.results(DEprot.timecourse.object = dpo))
  expect_error(plot.timecourse.protein(DEprot.timecourse.object = dpo, protein.id = "protein.1"))
  expect_error(plot.timecourse.profiles(DEprot.timecourse.object = dpo))
  expect_error(heatmap.timecourse(DEprot.timecourse.object = dpo))
})


test_that("errors are returned for invalid argument values", {
  expect_error(rank.timecourse(DEprot.timecourse.object = tc, rank.by = "wrong.metric"))
  expect_error(heatmap.timecourse(DEprot.timecourse.object = tc, order.by = "wrong.order"))
  expect_error(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                       protein.id = "protein.1",
                                       values = "wrong.values"))
  expect_error(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                       protein.id = "not.a.protein"))
})



##########################################################################################
####                              ANALYZE.TIMECOURSE                                  ####
##########################################################################################

test_that("no error is returned when running a time-course analysis", {
  expect_no_error(analyze.timecourse(DEprot.object = dpo,
                                     time.column = "time.hours",
                                     n.clusters = 3,
                                     verbose = FALSE))
})


test_that("the output is an object of class DEprot.timecourse", {
  expect_s4_class(tc, "DEprot.timecourse")
  expect_true(is.data.frame(tc@results))
  expect_equal(nrow(tc@results), nrow(sim$counts))
  expect_equal(tc@timepoints, c(0, 1, 2, 6, 24))
  expect_equal(tc@data.used, "imputed")
})


test_that("the results table contains the expected columns", {
  expect_true(all(c("prot.id", "F.statistic", "p.value", "padj", "amplitude",
                    "initial.slope", "peak.time", "trend.shape", "trend.status",
                    "score", "cluster", "membership",
                    "rank.overall", "rank.in.cluster") %in% colnames(tc@results)))
})


test_that("the fitted curves and the grid are consistent", {
  expect_type(tc@fitted.curves, "list")
  expect_named(tc@fitted.curves, "all")
  expect_equal(ncol(tc@fitted.curves[["all"]]), tc@params$grid.n)
  expect_equal(length(tc@time.grid), tc@params$grid.n)

  ## the grid must be expressed on the ORIGINAL time scale, despite the log2 transform
  expect_equal(range(tc@time.grid), range(tc@timepoints))
})


test_that("the observed means have one column per timepoint", {
  expect_equal(ncol(tc@observed.means[["all"]]), length(tc@timepoints))
  expect_equal(colnames(tc@observed.means[["all"]]), as.character(tc@timepoints))
})


test_that("the trends are detected and the shapes match the simulated archetypes", {
  trending <- get.timecourse.results(tc)

  expect_gt(nrow(trending), 50)

  ## the responsive proteins must be strongly enriched among the trending ones
  archetypes <- sim$truth$archetype[match(trending$prot.id, sim$truth$prot.id)]
  expect_gt(mean(archetypes != "flat"), 0.9)

  ## a monotone.up protein cannot be called monotone.down
  up <- trending[trending$prot.id %in% sim$truth$prot.id[sim$truth$archetype == "monotone.up"],]
  expect_true(all(up$trend.shape != "monotone.down"))
})


test_that("the amplitude threshold moves proteins to 'unresponsive'", {
  tc.strict <- analyze.timecourse(DEprot.object = dpo,
                                  time.column = "time.hours",
                                  time.transform = "log2",
                                  log2.amplitude.th = 10,
                                  n.clusters = 0,
                                  verbose = FALSE)

  expect_equal(sum(tc.strict@results$trend.status == "trending"), 0)
  expect_gt(sum(tc.strict@results$trend.status == "unresponsive"), 0)
})


test_that("the spline degrees of freedom are capped on the number of timepoints", {
  ## 5 timepoints -> the maximum is 3
  expect_warning(tc.df <- analyze.timecourse(DEprot.object = dpo,
                                             time.column = "time.hours",
                                             spline.df = 8,
                                             n.clusters = 0,
                                             verbose = FALSE))
  expect_equal(tc.df@params$spline.df, 3)
})


test_that("a linear trend is fitted below four timepoints", {
  sim.3 <- DEprot:::.simulate.timecourse(n.proteins = 200,
                                         timepoints = c(0, 2, 24),
                                         n.replicates = 3,
                                         seed = 2)

  dpo.3 <- load.counts2(counts = sim.3$counts,
                        metadata = sim.3$metadata,
                        data.type = "imputed",
                        log.base = 2,
                        column.id = "column.id")

  expect_warning(tc.3 <- analyze.timecourse(DEprot.object = dpo.3,
                                            time.column = "time.hours",
                                            n.clusters = 0,
                                            verbose = FALSE))
  expect_equal(tc.3@params$spline.df, 1)
  expect_equal(tc.3@params$model, "linear")
})


test_that("the clustering can be skipped", {
  tc.noclust <- analyze.timecourse(DEprot.object = dpo,
                                   time.column = "time.hours",
                                   n.clusters = 0,
                                   verbose = FALSE)

  expect_null(tc.noclust@clusters)
  expect_true(all(is.na(tc.noclust@results$cluster)))
  expect_null(tc.noclust@profile.plot)
})


test_that("the clustering returns the requested number of clusters", {
  expect_equal(tc@clusters$k, 5)
  expect_equal(tc@clusters$method, "cmeans")
  expect_true(tc@clusters$fuzzifier > 1)
  expect_equal(length(unique(stats::na.omit(tc@results$cluster))), 5)

  ## the membership exists only for the trending proteins
  expect_true(all(!is.na(tc@results$membership[tc@results$trend.status == "trending"])))
  expect_true(all(is.na(tc@results$membership[tc@results$trend.status != "trending"])))
})


test_that("the hard clustering works and returns no membership", {
  tc.pam <- analyze.timecourse(DEprot.object = dpo,
                               time.column = "time.hours",
                               time.transform = "log2",
                               log2.amplitude.th = 0.5,
                               clustering.method = "pam",
                               n.clusters = 4,
                               verbose = FALSE)

  expect_equal(tc.pam@clusters$method, "pam")
  expect_true(all(is.na(tc.pam@results$membership)))
  expect_no_error(plot.timecourse.profiles(DEprot.timecourse.object = tc.pam))
})


test_that("the analysis is reproducible with the same seed", {
  tc.a <- analyze.timecourse(DEprot.object = dpo, time.column = "time.hours",
                             n.clusters = 4, seed = 99, verbose = FALSE)
  tc.b <- analyze.timecourse(DEprot.object = dpo, time.column = "time.hours",
                             n.clusters = 4, seed = 99, verbose = FALSE)

  expect_equal(tc.a@results$padj, tc.b@results$padj)
  expect_equal(tc.a@results$cluster, tc.b@results$cluster)
})


test_that("all the time transformations are accepted and reversible", {
  for (tr in c("none", "log2", "log10", "log1p", "sqrt", "LOG2")) {
    tc.tr <- analyze.timecourse(DEprot.object = dpo,
                                time.column = "time.hours",
                                time.transform = tr,
                                n.clusters = 0,
                                verbose = FALSE)

    ## the back-transformed grid must return to the original time range
    expect_false(is.null(tc.tr@time.grid))
    expect_equal(range(tc.tr@time.grid), range(tc@timepoints), tolerance = 1e-6)
  }
})


test_that("all the count types are accepted through their aliases", {
  expect_no_error(analyze.timecourse(DEprot.object = dpo, time.column = "time.hours",
                                     which.data = "imp", n.clusters = 0, verbose = FALSE))
  expect_error(analyze.timecourse(DEprot.object = dpo, time.column = "time.hours",
                                  which.data = "wrong.counts", n.clusters = 0, verbose = FALSE))
})



##########################################################################################
####                            GROUP x TIME INTERACTION                              ####
##########################################################################################

test_that("the interaction model runs and returns group-specific descriptors", {
  expect_s4_class(tc.groups, "DEprot.timecourse")
  expect_equal(sort(tc.groups@params$group.levels), c("control", "treated"))
  expect_named(tc.groups@fitted.curves, c("control", "treated"))

  expect_true(all(c("amplitude.treated", "amplitude.control",
                    "peak.time.treated", "peak.time.control") %in% colnames(tc.groups@results)))

  ## the global amplitude is the largest one across the groups
  expect_equal(tc.groups@results$amplitude,
               pmax(tc.groups@results$amplitude.treated, tc.groups@results$amplitude.control))
})


test_that("the interaction detects the proteins moving only in the treated arm", {
  trending <- get.timecourse.results(tc.groups)
  archetypes <- sim.groups$truth$archetype[match(trending$prot.id, sim.groups$truth$prot.id)]

  expect_gt(nrow(trending), 20)
  expect_gt(mean(archetypes != "flat"), 0.9)
})


test_that("the design is rejected when the samples are not enough", {
  sim.small <- DEprot:::.simulate.timecourse(n.proteins = 100,
                                             timepoints = c(0, 1, 2, 6),
                                             n.replicates = 1,
                                             groups = c("treated", "control"),
                                             seed = 3)

  dpo.small <- load.counts2(counts = sim.small$counts,
                            metadata = sim.small$metadata,
                            data.type = "imputed",
                            log.base = 2,
                            column.id = "column.id")

  expect_error(analyze.timecourse(DEprot.object = dpo.small,
                                  time.column = "time.hours",
                                  group.column = "group",
                                  spline.df = 2,
                                  verbose = FALSE))
})



##########################################################################################
####                            RANKING AND GETTERS                                   ####
##########################################################################################

test_that("the ranking is recomputed without refitting", {
  tc.memb <- rank.timecourse(DEprot.timecourse.object = tc, rank.by = "membership")

  expect_s4_class(tc.memb, "DEprot.timecourse")
  expect_equal(tc.memb@params$rank.by, "membership")

  ## only the ranking columns change
  expect_equal(tc.memb@results$padj, tc@results$padj)
  expect_equal(tc.memb@results$cluster, tc@results$cluster)
  expect_false(identical(tc.memb@results$rank.in.cluster, tc@results$rank.in.cluster))
})


test_that("the ranking falls back to the score when no membership is available", {
  tc.pam <- analyze.timecourse(DEprot.object = dpo, time.column = "time.hours",
                               log2.amplitude.th = 0.5, clustering.method = "pam",
                               n.clusters = 3, verbose = FALSE)

  expect_warning(tc.pam <- rank.timecourse(DEprot.timecourse.object = tc.pam,
                                           rank.by = "membership"))
  expect_equal(tc.pam@params$rank.by, "score")
})


test_that("the rankings start at 1 within each cluster and globally", {
  for (k in sort(unique(stats::na.omit(tc@results$cluster)))) {
    expect_equal(min(tc@results$rank.in.cluster[which(tc@results$cluster == k)]), 1)
  }
  expect_equal(min(tc@results$rank.overall, na.rm = TRUE), 1)
})


test_that("the getter subsets and sorts the results", {
  all.res <- get.timecourse.results(DEprot.timecourse.object = tc)
  expect_true(all(all.res$trend.status == "trending"))
  expect_equal(all.res$rank.overall, sort(all.res$rank.overall))

  cl2 <- get.timecourse.results(DEprot.timecourse.object = tc, cluster = 2)
  expect_true(all(cl2$cluster == 2))
  expect_equal(cl2$rank.in.cluster, sort(cl2$rank.in.cluster))

  top5 <- get.timecourse.results(DEprot.timecourse.object = tc, cluster = 2, top.n = 5)
  expect_equal(nrow(top5), 5)

  ## with trending.only = FALSE the whole table is returned
  expect_equal(nrow(get.timecourse.results(DEprot.timecourse.object = tc,
                                           trending.only = FALSE)),
               nrow(tc@results))
})



test_that("t.half is a pure timing parameter", {
  trending <- get.timecourse.results(tc)
  expect_true(all(trending$t.half >= min(tc@timepoints) &
                    trending$t.half <= max(tc@timepoints), na.rm = TRUE))
  ## early archetypes must be faster than late ones
  early <- trending$t.half[trending$trend.shape == "transient.up"]
  expect_lt(median(early, na.rm = TRUE), max(tc@timepoints) / 2)
})



##########################################################################################
####                              METHODS AND PLOTS                                   ####
##########################################################################################

test_that("the show and summary methods work", {
  expect_no_error(capture.output(show(tc)))
  expect_no_error(capture.output(show(tc.groups)))

  smr <- summary(tc)
  expect_true(is.data.frame(smr))
  expect_equal(nrow(smr), tc@clusters$k)
  expect_true(all(c("cluster", "n", "dominant.shape",
                    "median.amplitude", "median.peak.time", "top.protein") %in% colnames(smr)))
})


test_that("the profile plot is stored in the object and returned by plot()", {
  expect_s3_class(tc@profile.plot, "ggplot")
  expect_s3_class(plot(tc), "ggplot")
})


test_that("no error is returned by the individual protein plots", {
  best <- get.timecourse.results(DEprot.timecourse.object = tc, top.n = 1)$prot.id

  for (v in c("counts", "log2FC", "zscore")) {
    expect_s3_class(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                            protein.id = best,
                                            values = v),
                    "ggplot")
  }

  ## several proteins, custom color, shapes and log axis
  multi <- get.timecourse.results(DEprot.timecourse.object = tc, top.n = 4)$prot.id
  expect_s3_class(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                          protein.id = multi,
                                          color = "steelblue",
                                          shape.column = "replicate",
                                          log.x = TRUE,
                                          ncol = 2),
                  "ggplot")

  ## the groups get their own curve
  expect_s3_class(plot.timecourse.protein(DEprot.timecourse.object = tc.groups,
                                          protein.id = get.timecourse.results(tc.groups, top.n = 1)$prot.id,
                                          values = "log2FC",
                                          group.colors = c(treated = "red", control = "grey40")),
                  "ggplot")
})


test_that("the deprecated scale.expression argument still works", {
  best <- get.timecourse.results(DEprot.timecourse.object = tc, top.n = 1)$prot.id

  expect_warning(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                         protein.id = best,
                                         scale.expression = TRUE))
})


test_that("no error is returned by the cluster profile plots", {
  for (v in c("zscore", "log2FC", "counts")) {
    expect_s3_class(plot.timecourse.profiles(DEprot.timecourse.object = tc, values = v), "ggplot")
  }

  expect_s3_class(plot.timecourse.profiles(DEprot.timecourse.object = tc,
                                           clusters = c(1, 2),
                                           top.n = 20,
                                           line.color = "grey60",
                                           centroid.color = "red",
                                           log.x = TRUE),
                  "ggplot")

  expect_s3_class(plot.timecourse.profiles(DEprot.timecourse.object = tc.groups), "ggplot")

  ## nothing left to plot
  expect_error(plot.timecourse.profiles(DEprot.timecourse.object = tc, clusters = 99))
})


test_that("no error is returned by the heatmaps", {
  for (v in c("zscore", "log2FC", "counts")) {
    expect_s3_class(heatmap.timecourse(DEprot.timecourse.object = tc, values = v), "ggplot")
  }

  for (o in c("rank", "membership", "peak.time", "amplitude", "hclust")) {
    expect_s3_class(heatmap.timecourse(DEprot.timecourse.object = tc,
                                       top.n = 20,
                                       order.by = o),
                    "ggplot")
  }

  expect_s3_class(heatmap.timecourse(DEprot.timecourse.object = tc,
                                     use.fitted = TRUE,
                                     order.by = "peak.time",
                                     panel.border = FALSE,
                                     show.protein.names = TRUE,
                                     protein.names.pattern = "protein\\."),
                  "ggplot")

  expect_s3_class(heatmap.timecourse(DEprot.timecourse.object = tc.groups), "ggplot")
  expect_s3_class(heatmap.timecourse(DEprot.timecourse.object = tc.groups,
                                     group.subset = "treated"),
                  "ggplot")

  expect_error(heatmap.timecourse(DEprot.timecourse.object = tc.groups,
                                  group.subset = "not.a.group"))
})


test_that("the reference time is snapped to the closest timepoint", {
  best <- get.timecourse.results(DEprot.timecourse.object = tc, top.n = 1)$prot.id

  expect_warning(plot.timecourse.protein(DEprot.timecourse.object = tc,
                                         protein.id = best,
                                         values = "log2FC",
                                         reference.time = 3))
})



##########################################################################################
####                             ENRICHMENT PER CLUSTER                               ####
##########################################################################################

test_that("errors are returned by the enrichment with a wrong input", {
  expect_error(timecourse.enrichment(DEprot.timecourse.object = dpo,
                                     TERM2GENE = sim$TERM2GENE))

  expect_error(timecourse.enrichment(DEprot.timecourse.object = tc,
                                     TERM2GENE = sim$TERM2GENE,
                                     size.by = "wrong.metric"))

  expect_error(timecourse.enrichment(DEprot.timecourse.object = tc,
                                     TERM2GENE = sim$TERM2GENE,
                                     clusters = 99))
})


test_that("the enrichment returns a DEprot.timecourse.enrichment object", {
  ora <- timecourse.enrichment(DEprot.timecourse.object = tc,
                               TERM2GENE = sim$TERM2GENE,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               dotplot.n = 3)

  expect_s4_class(ora, "DEprot.timecourse.enrichment")
  expect_true(is.data.frame(ora@results))
  expect_true("cluster" %in% colnames(ora@results))
  expect_true("FoldEnrichment" %in% colnames(ora@results))

  ## the default universe is made of all the quantified proteins
  expect_equal(length(ora@universe), nrow(tc@counts.used))
  expect_equal(ora@parameters$universe.size, nrow(tc@counts.used))

  expect_s3_class(ora@dotplot, "ggplot")
  expect_s3_class(plot(ora), "ggplot")
  expect_no_error(capture.output(show(ora)))
})


test_that("the clusters too small are skipped", {
  expect_warning(timecourse.enrichment(DEprot.timecourse.object = tc,
                                       TERM2GENE = sim$TERM2GENE,
                                       min.cluster.size = 1e5))
})
