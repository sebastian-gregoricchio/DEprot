## ----setup, include = FALSE---------------------------------------------------
# knitr::opts_chunk$set(collapse = TRUE, comment = ">", dev = "svg",
#                       warning = F, message = F, fig.align = "center",
#                       rows.print=12)
#knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

if (Sys.info()[["sysname"]] == "Darwin") {
  knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                        warning = F, message = F, fig.align = "center",
                        rows.print = 12, dev = "png", dpi = 96)
} else {
  knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                        warning = F, message = F, fig.align = "center",
                        rows.print = 12, dev = "png", dpi = 96)
}

#options(tibble.print_min = 4L, tibble.print_max = 4L)

# Load libraries required
require(DEprot)
require(dplyr)
require(ggplot2)

## ----simulate-----------------------------------------------------------------
sim <- DEprot:::.simulate.timecourse(n.proteins = 1000,
                                     timepoints = c(0, 1, 2, 6, 24),
                                     n.replicates = 3,
                                     fraction.responsive = 0.3,
                                     seed = 1234)

names(sim)

## ----print_metadata, echo=FALSE-----------------------------------------------
knitr::kable(head(sim$metadata, 6), row.names = F,
             caption = "**Sample metadata table** (first 6 rows)")

## ----print_counts, echo=FALSE-------------------------------------------------
knitr::kable(data.frame(round(sim$counts[1:5, 1:6], 2)), row.names = T,
             caption = "**Simulated log2(LFQ) values**")

## ----make_dpo-----------------------------------------------------------------
dpo <- load.counts2(counts = sim$counts,
                    metadata = sim$metadata,
                    data.type = "imputed",
                    log.base = 2,
                    column.id = "column.id")

dpo

## ----analyze------------------------------------------------------------------
tc <- analyze.timecourse(DEprot.object = dpo,
                         time.column = "time.hours",
                         time.transform = "log2",
                         which.data = "imputed",
                         padj.th = 0.05,
                         log2.amplitude.th = 0.5,
                         n.clusters = 5,
                         seed = 1234)

tc

## ----results_table, eval=F----------------------------------------------------
# head(get.timecourse.results(tc, top.n = 8))

## ----print_results, echo=FALSE------------------------------------------------
res.head <- get.timecourse.results(tc, top.n = 8)
knitr::kable(res.head[,c("prot.id","F.statistic","padj","amplitude","initial.slope","t.half","peak.time","trend.shape","cluster","membership")],
             row.names = F, digits = 3, caption = "**Top 8 trending proteins**")

## ----truth_check, eval=F------------------------------------------------------
# trending <- get.timecourse.results(tc)
# 
# table(sim$truth$archetype[match(trending$prot.id, sim$truth$prot.id)],
#       trending$trend.shape)

## ----print_truth, echo=FALSE--------------------------------------------------
trending <- get.timecourse.results(tc)
knitr::kable(as.data.frame.matrix(table(simulated = sim$truth$archetype[match(trending$prot.id, sim$truth$prot.id)],
                                        detected = trending$trend.shape)),
             caption = "**Simulated archetype (rows) versus detected shape (columns)**")

## ----plot_protein, fig.width=5, fig.height=4.2--------------------------------
best.protein <- get.timecourse.results(tc, top.n = 1)$prot.id

plot.timecourse.protein(DEprot.timecourse.object = tc,
                        protein.id = best.protein,
                        shape.column = "replicate",
                        log.x = TRUE)

## ----plot_protein_fc, fig.width=5, fig.height=4.2-----------------------------
plot.timecourse.protein(DEprot.timecourse.object = tc,
                        protein.id = best.protein,
                        values = "log2FC",
                        reference.time = 0,
                        log.x = TRUE)

## ----plot_protein_multi, fig.width=7, fig.height=5----------------------------
top.per.cluster <-
  sapply(1:5, function(k) {get.timecourse.results(tc, cluster = k, top.n = 1)$prot.id})

plot.timecourse.protein(DEprot.timecourse.object = tc,
                        protein.id = top.per.cluster,
                        values = "log2FC",
                        show.points = FALSE,
                        log.x = TRUE,
                        ncol = 3)

## ----plot_profiles, fig.width=7.5, fig.height=5-------------------------------
plot.timecourse.profiles(DEprot.timecourse.object = tc,
                         log.x = TRUE)

## ----plot_profiles_fc, fig.width=7.5, fig.height=5----------------------------
plot.timecourse.profiles(DEprot.timecourse.object = tc,
                         values = "log2FC",
                         reference.time = 0,
                         top.n = 100,
                         log.x = TRUE)

## ----summary_tc, eval=F-------------------------------------------------------
# summary(tc)

## ----print_summary, echo=FALSE------------------------------------------------
knitr::kable(summary(tc), row.names = F, digits = 3, caption = "**Cluster summary**")

## ----heatmap_tc, fig.width=6, fig.height=7------------------------------------
heatmap.timecourse(DEprot.timecourse.object = tc,
                   values = "zscore",
                   order.by = "peak.time")

## ----heatmap_fitted, fig.width=6, fig.height=7--------------------------------
heatmap.timecourse(DEprot.timecourse.object = tc,
                   values = "log2FC",
                   reference.time = 0,
                   use.fitted = TRUE,
                   order.by = "peak.time")

## ----ranking, eval=F----------------------------------------------------------
# # strongest movers of cluster 1
# get.timecourse.results(tc, cluster = 1, top.n = 5)
# 
# # most prototypical proteins of cluster 1
# tc.membership <- rank.timecourse(tc, rank.by = "membership")
# get.timecourse.results(tc.membership, cluster = 1, top.n = 5)

## ----print_ranking, echo=FALSE------------------------------------------------
tc.membership <- rank.timecourse(tc, rank.by = "membership")

knitr::kable(get.timecourse.results(tc, cluster = 1, top.n = 5)[,c("prot.id","amplitude","padj","score","membership")],
             row.names = F, digits = 3, caption = "**Cluster 1 ranked by score**")

knitr::kable(get.timecourse.results(tc.membership, cluster = 1, top.n = 5)[,c("prot.id","amplitude","padj","score","membership")],
             row.names = F, digits = 3, caption = "**Cluster 1 ranked by membership**")

## ----ordering, eval=F---------------------------------------------------------
# # the clusters, in the order in which they respond
# summary(tc)[,c("cluster", "n", "median.t.half", "dominant.shape")]
# 
# # within a wave, the proteins ordered by when they respond
# tc.timing <- rank.timecourse(tc, rank.by = "t.half")
# get.timecourse.results(tc.timing, cluster = 1, top.n = 10)

## ----print_ordering, echo=FALSE-----------------------------------------------
tc.timing <- rank.timecourse(tc, rank.by = "t.half")

knitr::kable(summary(tc)[,c("cluster", "n", "median.t.half", "median.amplitude", "dominant.shape")],
             row.names = F, digits = 3,
             caption = "**Clusters ordered by median time to half-maximum**")

knitr::kable(get.timecourse.results(tc.timing, cluster = 1, top.n = 8)[,c("prot.id","t.half","initial.slope","amplitude","padj")],
             row.names = F, digits = 3, caption = "**Cluster 1 ranked by t.half**")

## ----enrichment---------------------------------------------------------------
ora <- timecourse.enrichment(DEprot.timecourse.object = tc,
                             TERM2GENE = sim$TERM2GENE,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             dotplot.n = 3)

ora

## ----enrichment_dotplot, fig.width=7, fig.height=5----------------------------
plot(ora)

## ----enrichment_results, eval=F-----------------------------------------------
# head(ora@results[,c("cluster","ID","GeneRatio","FoldEnrichment","p.adjust","Count")])

## ----print_enrichment, echo=FALSE---------------------------------------------
knitr::kable(head(ora@results[,c("cluster","ID","GeneRatio","FoldEnrichment","p.adjust","Count")], 8),
             row.names = F, digits = 4, caption = "**Enrichment results (first 8 rows)**")

## ----simulate_groups----------------------------------------------------------
sim2 <- DEprot:::.simulate.timecourse(n.proteins = 800,
                                      timepoints = c(0, 1, 2, 6, 24),
                                      n.replicates = 3,
                                      groups = c("treated", "control"),
                                      fraction.responsive = 0.25,
                                      seed = 42)

dpo2 <- load.counts2(counts = sim2$counts,
                     metadata = sim2$metadata,
                     data.type = "imputed",
                     log.base = 2,
                     column.id = "column.id")

tc2 <- analyze.timecourse(DEprot.object = dpo2,
                          time.column = "time.hours",
                          group.column = "group",
                          time.transform = "log2",
                          padj.th = 0.05,
                          log2.amplitude.th = 0.5,
                          n.clusters = 4,
                          seed = 1234)

tc2

## ----plot_protein_groups, fig.width=5.5, fig.height=4.2-----------------------
plot.timecourse.protein(DEprot.timecourse.object = tc2,
                        protein.id = get.timecourse.results(tc2, top.n = 1)$prot.id,
                        values = "log2FC",
                        reference.time = 0,
                        log.x = TRUE)

## ----session_info-------------------------------------------------------------
sessionInfo()

