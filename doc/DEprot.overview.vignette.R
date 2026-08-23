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
require(legendry)
require(ggpubr)

## ----citation, message=FALSE, warning=FALSE-----------------------------------
citation("DEprot")

## ----read_metadata, eval=F----------------------------------------------------
# # Metadata
# data("sample.config", package = "DEprot")
# sample.config

## ----print_metadata, echo=FALSE-----------------------------------------------
data("sample.config", package = "DEprot")
knitr::kable(sample.config, row.names = F, caption = "**Sample metadata table**")

## ----read_counts, eval=F------------------------------------------------------
# # log2(LFQ) values (not imputed)
# data("unimputed.counts", package = "DEprot")
# head(unimputed.counts[,1:6])

## ----print_counts, echo=FALSE-------------------------------------------------
data("unimputed.counts", package = "DEprot")
knitr::kable(data.frame(head(unimputed.counts[,1:6])), row.names = T, caption = "**Unimputed log2(LFQ) values**")

## ----read_protein_info, eval=F------------------------------------------------
# # Protein annotation
# ## here a dummy example is generated using the same names of the counts row names
# protein.annotation <- data.frame(protein.name = rownames(unimputed.counts))
# rownames(protein.annotation) = rownames(unimputed.counts)
# 
# head(protein.annotation)

## ----print_protein_info, echo=FALSE-------------------------------------------
protein.annotation <- data.frame(protein.name = rownames(unimputed.counts))
rownames(protein.annotation) = rownames(unimputed.counts)
knitr::kable(head(protein.annotation), row.names = TRUE, caption = "**Protein annotation table**")

## ----make_dpo_info------------------------------------------------------------
dpo <- load.counts2(counts = unimputed.counts,
                    metadata = sample.config,
                    data.type = "raw",
                    log.base = 2,
                    column.id = "column.id",
                    protein.info = protein.annotation) # set as NULL if not available

dpo

## ----add_protein_info, eval=F-------------------------------------------------
# dpo <- add.protein.info(DEprot.object = dpo,
#                         protein.info = protein.annotation)

## ----get_protein_info, eval=F-------------------------------------------------
# head(get.protein.info(dpo))

## ----print_get_protein_info, echo=FALSE---------------------------------------
knitr::kable(head(get.protein.info(dpo)), row.names = TRUE, caption = "**Protein annotation stored in the object**")

## ----rename_samples-----------------------------------------------------------
dpo <- rename.samples(DEprot.object = dpo,
                      metadata.column = "sample.id")

get.metadata(dpo)

## ----display_rename_samples---------------------------------------------------
head(dpo@raw.counts[,1:6])

## ----display_boxplot_raw, fig.width=5-----------------------------------------
dpo@boxplot.raw

## ----normalize----------------------------------------------------------------
dpo <- normalize.counts(DEprot.object = dpo,
                        NRI.RI.ratio.threshold = 0.5,
                        balancing.function = "median")

dpo

## ----normalize2---------------------------------------------------------------
dpo@normalization.method

## ----normalize3---------------------------------------------------------------
head(dpo@raw.counts[,1:6])

## ----plot_norm, fig.width=8---------------------------------------------------
plot(dpo, nrow = 1)

## ----eval = F-----------------------------------------------------------------
# ## Adding batch column to the metadata table
# dpo@metadata$batch = rep(c("batch_1","batch_2"), each = 6)
# 
# get.metadata(dpo)

## ----echo=FALSE---------------------------------------------------------------
new_meta = dpo@metadata
new_meta$batch = rep(c("batch_1","batch_2"), each = 6)
new_meta

## ----eval = F-----------------------------------------------------------------
# ## batch correction
# dpo <- harmonize.batches(DEprot.object = dpo,
#                          batch.column = "batch",
#                          cores = 1)

## ----eval = FALSE-------------------------------------------------------------
# dpo <- harmonize.batches(DEprot.object = dpo,
#                          batch.column = "batch",
#                          algorithm = "limma",
#                          cores = 1)

## ----randomize_data-----------------------------------------------------------
dpo <- randomize.missing.values(DEprot.object = dpo,
                                group.column = "combined.id",
                                percentage.missing = 100, # completely missing
                                tail.percentage = 3,
                                seed = 1234,
                                verbose = FALSE)

## ----missingness_diagnostic---------------------------------------------------
miss <- missingness.diagnostic(DEprot.object = dpo,
                               verbose = FALSE)

miss

## ----missingness_summary_plots, fig.width = 12, fig.height = 9----------------
plot(miss, ncol = 2)

## ----missingness_heatmap, fig.width = 8, fig.height = 9-----------------------
plot(miss, plot.type = "heatmap")

## ----choice_of_imp, fig.width = 15, fig.height = 12, message=FALSE, warning=FALSE----
imp.comparison <- compare.imp.methods(DEprot.object = dpo,
                                      percentage.test = 30,
                                      sample.group.column = "combined.id",
                                      which.data = "normalized",
                                      seed = 1234,
                                      run.kNN = FALSE, # time consuming
                                      verbose = FALSE)

patchwork::wrap_plots(c(imp.comparison@correlation.plots,
                        imp.comparison@density.residuals))

## ----summary_impute_comparison, echo=FALSE------------------------------------
summary(imp.comparison)

## ----imputation_example, eval = F---------------------------------------------
# ## Without parallelization
# dpo <- impute.counts(DEprot.object = dpo,
#                      method = "missForest",
#                      which.data = "randomized",
#                      missForest.max.iterations = 100,
#                      missForest.variable.wise.OOBerror = TRUE,
#                      seed = 1234)
# 
# 
# ## With parallelization
# dpo <- impute.counts(DEprot.object = dpo,
#                      method = "missForest",
#                      which.data = "randomized",
#                      missForest.max.iterations = 100,
#                      missForest.variable.wise.OOBerror = TRUE,
#                      missForest.cores = 20,
#                      missForest.parallel.mode = "variables",
#                      seed = 1234)
# 
# dpo
# 
# head(dpo@imputation.method$OOBerror)
# 
# data.frame(dpo@imputation.method[-3])

## ----load_imputation, echo=FALSE----------------------------------------------
data("dpo.imputed.counts", package = "DEprot")
dpo = dpo.imputed.counts
dpo

## ----load_imputation2, echo=FALSE---------------------------------------------
error = dpo@imputation.method$OOBerror
names(error) = colnames(dpo@imputed.counts)
head(error)

knitr::kable(data.frame(dpo@imputation.method[-3]), row.names = F)

## ----plot_imputed, fig.width=7, fig.height=8, fig.align='center'--------------
plot(dpo, ncol = 2)

## ----make_PCA, fig.width=8----------------------------------------------------
## Perform the analyses (DEprot.PCA object)
PCA <- perform.PCA(DEprot.object = dpo,
                   which.data = "imputed") # possible: raw, normalized, imputed

## ----run_cumulative-----------------------------------------------------------
## Plot cumulative variance of all PCs
#### equivalent to `PCA@cumulative.PC.plot`
plot.PC.cumulative(DEprot.PCA.object = PCA,
                   bar.color = "steelblue",
                   line.color = "navyblue")

## ----run_PCA------------------------------------------------------------------
## Plot PC scatters
PC_1.2 <-
  plot.PC.scatter(DEprot.PCA.object = PCA,
                  PC.x = 1,
                  PC.y = 2,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.line.x = TRUE,
                  plot.zero.line.y = TRUE) +
  theme(legend.position = "none")


PC_2.3 <-
  plot.PC.scatter(DEprot.PCA.object = PCA,
                  PC.x = 2,
                  PC.y = 3,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.line.x = TRUE,
                  plot.zero.line.y = TRUE)

patchwork::wrap_plots(PC_1.2, PC_2.3, nrow = 1)

## ----run_plotPCscatter123-----------------------------------------------------
plot.PC.scatter.123(DEprot.PCA.object = PCA,
                    color.column = "condition",
                    shape.column = "replicate",
                    label.column = "replicate",
                    dot.colors = c("6h.10nM.E2" = "indianred",
                                   "6h.DMSO" = "steelblue",
                                   "FBS" =  "forestgreen"),
                    plot.zero.line.y.12 = TRUE,
                    plot.zero.line.x.12 = TRUE,
                    plot.zero.line.y.23 = TRUE,
                    plot.zero.line.x.23 = TRUE)

## ----run_biplot, fig.height=4-------------------------------------------------
PC_biplot_1.2 <-
  plot.PC.biplot(DEprot.PCA.object = PCA,
                 PC.x = 1,
                 PC.y = 2,
                 color.column = "condition",
                 shape.column = "replicate",
                 label.column = NULL,
                 n.loadings = 5,
                 plot.zero.line.x = TRUE,
                 plot.zero.line.y = TRUE)

PC_biplot_1.2

## ----make_PCA_subset, fig.width=8---------------------------------------------
## Perform the analyses (DEprot.PCA object)
PCA.fbs.e2 <-
  perform.PCA(DEprot.object = dpo,
              sample.subset = dpo@metadata$column.id[grepl("E2|FBS",
                                                           dpo@metadata$column.id)],
              which.data = "imputed")


## Plot cumulative variance of all PCs
plot.PC.cumulative(DEprot.PCA.object = PCA.fbs.e2,
                   bar.color = "indianred",
                   line.color = "firebrick4",
                   title = "**Only ERa active**")

## ----make_PCA_scatters_subset, fig.width=8------------------------------------
## Plot PC scatters
PC.fbs.e2_1.2 <-
  plot.PC.scatter(DEprot.PCA.object = PCA.fbs.e2,
                  PC.x = 1,
                  PC.y = 2,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.line.x = TRUE,
                  plot.zero.line.y = TRUE) +
  theme(legend.position = "none")
  
  
PC.fbs.e2_2.3 <-
  plot.PC.scatter(DEprot.PCA.object = PCA.fbs.e2,
                  PC.x = 2,
                  PC.y = 3,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.line.x = TRUE,
                  plot.zero.line.y = TRUE)

patchwork::wrap_plots(PC.fbs.e2_1.2, PC.fbs.e2_2.3, nrow = 1)

## ----make_correlation_all, fig.width=9----------------------------------------
corr.all.samples <-
  plot.correlation.heatmap(DEprot.object = dpo,
                           which.data = "imputed",
                           palette = viridis::mako(n = 10, direction = -1, begin = 0.25),
                           correlation.scale.limits = c(0.9,1),
                           correlation.method = "pearson",
                           plot.subtitle = "All samples",
                           display.values = TRUE)
corr.all.samples

## ----make_correlation_subset--------------------------------------------------
corr.ERa.active <-
  plot.correlation.heatmap(DEprot.object = dpo,
                           which.data = "imputed",
                           sample.subset = dpo@metadata$column.id[grepl("E2|FBS",
                                                                        dpo@metadata$column.id)],
                           palette = viridis::magma(n = 10, direction = -1, begin = 0.25),
                           correlation.scale.limits = c(0.9,1),
                           correlation.method = "pearson",
                           plot.subtitle = "Only ERa active",
                           clustering.method = "complete",
                           display.values = TRUE)

corr.ERa.active

## ----pcoa.compute-------------------------------------------------------------
## Perform the analyses (DEprot.PCoA object)
PCoA.ERa.active <- perform.PCoA(DEprot.correlation.object = corr.ERa.active,
                                DEprot.object = dpo) # needed only for the biplot

## ----pcoa.shepard-------------------------------------------------------------
PCoA.ERa.active@shepard.plot

## ----pcoa.cumulative----------------------------------------------------------
## Plot cumulative variance of all PCos
#### equivalent to `PCoA.ERa.active@cumulative.PCo.plot`
plot.PCoA.cumulative(DEprot.PCoA.object = PCoA.ERa.active,
                     bar.color = "indianred",
                     line.color = "firebrick4",
                     title = "**Only ERa active**")

## ----pcoa.scatter-------------------------------------------------------------
## Combined PCo1-vs-PCo2 and PCo3-vs-PCo2 scatters
#### equivalent to `PCoA.ERa.active@scatter.plot.123`
plot.PCoA.scatter.123(DEprot.PCoA.object = PCoA.ERa.active,
                      color.column = "condition",
                      shape.column = "replicate",
                      label.column = "replicate",
                      dot.colors = c("6h.10nM.E2" = "indianred",
                                     "FBS" = "forestgreen"))

## ----pcoa.biplot--------------------------------------------------------------
## Top proteins associated with the separation (projection onto the PCos)
PCoA_biplot_1.2 <- plot.PCoA.biplot(DEprot.PCoA.object = PCoA.ERa.active,
                                    PCo.x = 1,
                                    PCo.y = 2,
                                    color.column = "condition",
                                    shape.column = "replicate",
                                    label.column = NULL,
                                    n.loadings = 5)
PCoA_biplot_1.2

## ----detect_outliers----------------------------------------------------------
outliers <- detect.outliers(DEprot.object = dpo,
                            which.data = "imputed",
                            group.column = "condition",
                            correlation.method = "pearson",
                            n.PCs = 3,
                            min.flags = 2)

outliers

## ----plot_outliers, fig.height=9, fig.width=7---------------------------------
plot(outliers)

## ----filter_outliers, eval=FALSE----------------------------------------------
# dpo.clean <- filter.samples(DEprot.object = dpo,
#                             samples = outliers@outliers,
#                             mode = "remove")

## ----compute_diff_exp_examples_Ttest, eval=F----------------------------------
# ## Unpaired test
# dpo_analyses <- diff.analyses(DEprot.object = dpo,
#                               contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                    c("condition", "6h.10nM.E2", "FBS")),
#                               linear.FC.th = 2,
#                               padj.th = 0.05,
#                               padj.method = "BH",
#                               stat.test = "t.test",
#                               which.data = "imputed")
# 
# ## Paired test
# dpo_analyses <- diff.analyses(DEprot.object = dpo,
#                               contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                    c("condition", "6h.10nM.E2", "FBS")),
#                               replicate.column = "replicate",
#                               paired.test = TRUE,
#                               linear.FC.th = 2,
#                               padj.th = 0.05,
#                               padj.method = "BH",
#                               stat.test = "t.test",
#                               which.data = "imputed")
# 
# dpo_analyses

## ----compute_diff_exp_paired, echo=FALSE--------------------------------------
## Paired test
dpo_analyses <- diff.analyses(DEprot.object = dpo,
                              contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
                                                   c("condition", "6h.10nM.E2", "FBS")),
                              replicate.column = "replicate",
                              paired.test = TRUE,
                              linear.FC.th = 2,
                              padj.th = 0.05,
                              padj.method = "BH",
                              stat.test = "t.test",
                              which.data = "imputed")

dpo_analyses

## ----analyses_summary, eval=F-------------------------------------------------
# diff.analyses_summary <- summary(dpo)

## ----normality_test, message=T, fig.width=10, fig.height=10-------------------
normality <- check.normality(DEprot.object = dpo,
                             p.threshold = 0.05,
                             which.data = "imputed",
                             verbose = TRUE)

## ----normality_test_plots-----------------------------------------------------
## example of Q-Q and density plots
plot(normality, n.samples = 1)

## ----compute_diff_exp_examples_limma, eval=F----------------------------------
# ## Unpaired test
# dpo_analyses <- diff.analyses.limma(DEprot.object = dpo,
#                                     contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                          c("condition", "6h.10nM.E2", "FBS")),
#                                     linear.FC.th = 2,
#                                     padj.th = 0.05,
#                                     padj.method = "BH",
#                                     fitting.method = "ls",
#                                     which.data = "imputed")
# 
# ## Paired test
# dpo_analyses <- diff.analyses.limma(DEprot.object = dpo,
#                                     contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                          c("condition", "6h.10nM.E2", "FBS")),
#                                     replicate.column = "replicate",
#                                     include.rep.model = TRUE,
#                                     linear.FC.th = 2,
#                                     padj.th = 0.05,
#                                     padj.method = "BH",
#                                     fitting.method = "ls",
#                                     which.data = "imputed")

## ----prolfqua_analyses, eval = FALSE------------------------------------------
# dpo_prolfqua <-
#   diff.analyses.prolfqua(DEprot.object = dpo_imputed,
#                          contrast.list = list(c("condition", "FBS", "6h.DMSO"),
#                                               c("condition", "6h.10nM.E2", "6h.DMSO")),
#                          strategy = "lm",
#                          linear.FC.th = 1.5,
#                          FDR.th = 0.05,
#                          which.data = "imputed")

## ----prolfqua_lmer, eval = FALSE----------------------------------------------
# dpo_prolfqua <-
#   diff.analyses.prolfqua(DEprot.object = dpo_imputed,
#                          contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO")),
#                          strategy = "lmer",
#                          replicate.column = "replicate",
#                          linear.FC.th = 1.5)

## ----prolfqua_scaling, eval = FALSE-------------------------------------------
# dpo_prolfqua <-
#   diff.analyses.prolfqua(DEprot.object = dpo_norm,
#                          contrast.list = list(c("condition", "FBS", "6h.DMSO")),
#                          strategy = "lm",
#                          robust.scaling = FALSE,
#                          which.data = "normalized")

## ----prolfqua_scaling_factors, eval = FALSE-----------------------------------
# dpo_prolfqua@analyses.result.list$condition_FBS.vs.6h.DMSO$prolfqua.out$scaling.factors

## ----prolfqua_unimputed, eval = FALSE-----------------------------------------
# dpo_prolfqua <-
#   diff.analyses.prolfqua(DEprot.object = dpo_norm,
#                          contrast.list = list(c("condition", "FBS", "6h.DMSO")),
#                          strategy = "lm",
#                          robust.scaling = FALSE,
#                          which.data = "normalized")

## ----proDA_analyses, eval = FALSE---------------------------------------------
# dpo_proDA <-
#   diff.analyses.proDA(DEprot.object = dpo_norm,
#                       contrast.list = list(c("condition", "FBS", "6h.DMSO"),
#                                            c("condition", "6h.10nM.E2", "6h.DMSO")),
#                       linear.FC.th = 1.5,
#                       padj.th = 0.05,
#                       which.data = "normalized")

## ----proDA_missingness, eval = FALSE------------------------------------------
# miss <- missingness.diagnostic(DEprot.object = dpo_norm)
# 
# dpo_proDA <-
#   diff.analyses.proDA(DEprot.object = dpo_norm,
#                       contrast.list = list(c("condition", "FBS", "6h.DMSO")),
#                       missingness.object = miss,
#                       linear.FC.th = 1.5)

## ----proDA_paired, eval = FALSE-----------------------------------------------
# dpo_proDA <-
#   diff.analyses.proDA(DEprot.object = dpo_norm,
#                       contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO")),
#                       include.rep.model = TRUE,
#                       replicate.column = "replicate",
#                       linear.FC.th = 1.5)

## ----proDA_dropout, eval = FALSE----------------------------------------------
# dpo_proDA@analyses.result.list$condition_FBS.vs.6h.DMSO$proDA.fit$dropout.curves

## ----get_results, eval = F----------------------------------------------------
# ## Direct access
# results = dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$results
# 
# ## Function
# results = get.results(dpo_analyses, contrast = 1, protein.info.columns = "all")
# 
# head(results)

## ----get_results2, echo=FALSE-------------------------------------------------
knitr::kable(get.results(dpo_analyses, contrast = 1, protein.info.columns = "all")[1:6,], row.names = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# res <- get.results(dpo, contrast = 1)
# 
# confidence <- 0.95
# alpha <- 1 - confidence
# t.crit <- qt(1 - alpha/2, df = res$df)   # upper-tail critical value, per protein
# 
# res$CI.lower <- res$log2.Fold_FBS.vs.6h.DMSO - t.crit * res$lfcSE
# res$CI.upper <- res$log2.Fold_FBS.vs.6h.DMSO + t.crit * res$lfcSE

## ----DE_PCA_scatters, fig.width=8, eval=F-------------------------------------
# dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.plots

## ----DE_PCA_scatters_replotting, echo=FALSE, fig.width=8----------------------
scatter_PC12 = plot.PC.scatter(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data, 1,2, color.column = "condition") + theme(legend.position = "none")
scatter_PC23 = plot.PC.scatter(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data, 2,3, color.column = "condition")
cumulative = plot.PC.cumulative(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data)

scatters = cowplot::plot_grid(scatter_PC12, scatter_PC23, nrow = 1, align = "hv", axis = "tblr")
cowplot::plot_grid(scatters, cumulative, ncol = 1, axis = "tblr")

## ----DE_correlations, fig.width=10--------------------------------------------
dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$correlations

## ----check_pvalues, fig.width = 8, fig.height = 4-----------------------------
pval.distribution <- check.pvalues(DEprot.analyses.object = dpo_analyses,
                                   contrast = 2,
                                   histogram.binwidth = 0.025)

pval.distribution

## ----DE_volcano_MA, fig.width=9-----------------------------------------------
volcano <- plot.volcano(dpo_analyses,
                        contrast = 1,
                        label.top.n = 10,
                        label.font.size = 4,
                        use.uncorrected.pvalue = TRUE)

MAplot <- plot.MA(dpo_analyses,
                  contrast = 1,
                  use.uncorrected.pvalue = TRUE)

patchwork::wrap_plots(volcano, MAplot)

## ----contrast_scatter, fig.width=5, fig.height=5------------------------------
contrast.scatter <-
  contrast.scatter(DEprot.analyses.object = dpo_analyses,
                   contrast.x = 1,
                   contrast.y = 2,
                   regression.line.color = "firebrick",
                   correlation.method = "pearson",
                   add.foldchange.threshold = TRUE,
                   symmetric.axes = TRUE)

contrast.scatter

## ----contrast_LFQ, fig.width=5, fig.height=5----------------------------------
contrast_LFQ <-
  contrast.LFQ(DEprot.analyses.object = dpo_analyses,
               contrast = 2,
               dot.labels = "protein.3081")

contrast_LFQ

## ----heatmap_counts, fig.width=13, fig.height=4.5-----------------------------
## Plotting from a DEprot object
imputed_counts_heatmap <- 
  heatmap.counts(DEprot.object = dpo,
                 which.data = "imputed",
                 sample.subset = dpo@metadata$column.id[grep("6h", dpo@metadata$column.id)],
                 show.protein.names = TRUE,
                 protein.subset = c("protein.2295", "protein.304", "protein.657",
                                    "protein.2819", "protein.2168", "protein.10594"),
                 title = "Imputed counts | protein and sample selection") 



## Plotting from a DEprot.analyses object
## top 15 differential proteins from contrast 1
imputed_counts_heatmap_diffProteins <- 
  heatmap.counts(DEprot.object = dpo_analyses,
                 which.data = "imputed",
                 contrast = 1,
                 top.n = 15,
                 palette = viridis::mako(n = 100, direction = -1),
                 cell.border.color = "white",
                 show.protein.names = TRUE,
                 sample.subset = dpo@metadata$column.id[grep("6h", dpo@metadata$column.id)],
                 use.uncorrected.pvalue = TRUE,
                 protein.names.pattern = "protein[.]",
                 title = "condition: **6h.10nM.E2** *vs* **6h.DMSO** (top 15) | Imputed counts")


## Combine heatmaps
patchwork::wrap_plots(imputed_counts_heatmap@heatmap,
                      imputed_counts_heatmap_diffProteins@heatmap)

## ----heatmap_counts_zscores, fig.width=10, fig.height=4.5---------------------
## Z-score by row
imputed_counts_heatmap_diffProteins_rowScaled <- 
  heatmap.counts(DEprot.object = dpo_analyses,
                 which.data = "imputed",
                 contrast = 1,
                 top.n = 15,
                 high.color = "purple4",
                 low.color = "darkorange",
                 mid.color = "white",
                 cell.border.color = "white",
                 show.protein.names = TRUE,
                 sample.subset = dpo@metadata$column.id[grep("6h", dpo@metadata$column.id)],
                 use.uncorrected.pvalue = TRUE,
                 scale = "rows",
                 title = "condition: **6h.10nM.E2** *vs* **6h.DMSO** (top 15)<br>Imputed counts Z-score")


## Z-score by column
imputed_counts_heatmap_diffProteins_columnScaled <- 
  heatmap.counts(DEprot.object = dpo_analyses,
                 which.data = "imputed",
                 contrast = 1,
                 top.n = 15,
                 high.color = "firebrick",
                 low.color = "steelblue",
                 mid.color = "white",
                 cell.border.color = "white",
                 show.protein.names = TRUE,
                 sample.subset = dpo@metadata$column.id[grep("6h", dpo@metadata$column.id)],
                 use.uncorrected.pvalue = TRUE,
                 scale = "columns",
                 title = "condition: **6h.10nM.E2** *vs* **6h.DMSO** (top 15)<br>Imputed counts Z-score")


## Combine heatmaps
patchwork::wrap_plots(imputed_counts_heatmap_diffProteins_rowScaled@heatmap,
                      imputed_counts_heatmap_diffProteins_columnScaled@heatmap)


## ----heatmap_counts_groups, fig.width = 4, fig.height = 6---------------------
imputed_counts_heatmap_diffProteins_rowScaled_grouped.by.condition <- 
  heatmap.counts(DEprot.object = dpo_analyses,
                 group.by.metadata.column = "combined.id",
                 which.data = "imputed",
                 contrast = 1,
                 high.color = "firebrick",
                 low.color = "steelblue",
                 mid.color = "white",
                 cell.border.color = "white",
                 show.protein.names = TRUE,
                 use.uncorrected.pvalue = TRUE,
                 scale = "rows",
                 title = "condition: **6h.10nM.E2** *vs* **6h.DMSO** (all)<br>Imputed counts Z-score")

imputed_counts_heatmap_diffProteins_rowScaled_grouped.by.condition@heatmap

## ----heatmap_counts_anno, fig.width=7, fig.height=6---------------------------
counts_heatmap_annotated <-
  heatmap.counts.anno(DEprot.object = dpo_analyses,
                      which.data = "imputed",
                      contrast = 1,
                      top.n = 15,
                      scale = "row",
                      use.uncorrected.pvalue = TRUE,
                      column.annotation = c("condition", "replicate"),
                      show.protein.names = TRUE,
                      protein.names.pattern = "protein[.]",
                      title = "condition: 6h.10nM.E2 vs 6h.DMSO (top 15)")

counts_heatmap_annotated

## ----heatmap_counts_anno_colors, fig.width=7, fig.height=6--------------------
counts_heatmap_annotated_split <-
  heatmap.counts.anno(DEprot.object = dpo_analyses,
                      which.data = "imputed",
                      contrast = 1,
                      top.n = 15,
                      scale = "row",
                      use.uncorrected.pvalue = TRUE,
                      column.annotation = c("condition", "replicate"),
                      column.split = "condition",
                      annotation.colors = list(condition = c("6h.10nM.E2" = "indianred",
                                                             "6h.DMSO" = "steelblue",
                                                             "FBS" = "forestgreen")),
                      high.color = "purple4",
                      low.color = "darkorange",
                      mid.color = "white",
                      show.protein.names = TRUE,
                      protein.names.pattern = "protein[.]")

counts_heatmap_annotated_split

## ----heatmap_counts_anno_rows, fig.width=7, fig.height=6.5--------------------
## Protein annotation built on the results of the first contrast
results_contrast.1 <- get.results(DEprot.analyses.object = dpo_analyses, contrast = 1)

protein.annotation_contrast.1 <-
  data.frame(diff.status = results_contrast.1$diff.status,
             log2FC = results_contrast.1[,grep("^log2.Fold", colnames(results_contrast.1))],
             row.names = results_contrast.1$prot.id)

dpo_analyses_annotated <-
  add.protein.info(DEprot.object = dpo_analyses,
                   protein.info = protein.annotation_contrast.1,
                   overwrite = TRUE)


counts_heatmap_annotated_rows <-
  heatmap.counts.anno(DEprot.object = dpo_analyses_annotated,
                      which.data = "imputed",
                      contrast = 1,
                      top.n = 20,
                      scale = "row",
                      use.uncorrected.pvalue = TRUE,
                      column.annotation = "condition",
                      row.annotation = c("diff.status", "log2FC"),
                      row.split = "diff.status",
                      annotation.colors = list(condition = c("6h.10nM.E2" = "indianred",
                                                             "6h.DMSO" = "steelblue",
                                                             "FBS" = "forestgreen"),
                                               diff.status = c("6h.10nM.E2" = "indianred",
                                                               "6h.DMSO" = "steelblue",
                                                               "unresponsive" = "gray70",
                                                               "null" = "gray90"),
                                               log2FC = c("darkorange", "white", "purple4")),
                      show.protein.names = TRUE,
                      protein.names.pattern = "protein[.]")

counts_heatmap_annotated_rows

## ----heatmap_counts_anno_colorRamp, eval=FALSE--------------------------------
# annotation.colors <-
#   list(log2FC = circlize::colorRamp2(breaks = c(-2, 0, 2),
#                                      colors = c("darkorange", "white", "purple4")))

## ----heatmap_counts_anno_groups, fig.width=5.5, fig.height=6.5----------------
counts_heatmap_annotated_grouped <-
  heatmap.counts.anno(DEprot.object = dpo_analyses_annotated,
                      group.by.metadata.column = "combined.id",
                      which.data = "imputed",
                      contrast = 1,
                      top.n = 20,
                      scale = "row",
                      use.uncorrected.pvalue = TRUE,
                      column.annotation = c("condition", "cell"),
                      row.annotation = "diff.status",
                      show.protein.names = TRUE,
                      protein.names.pattern = "protein[.]")

counts_heatmap_annotated_grouped

## ----heatmap_counts_anno_draw, eval=FALSE-------------------------------------
# ## Legends on the left, all of them merged in a single block
# ComplexHeatmap::draw(counts_heatmap_annotated_rows@heatmap,
#                      heatmap_legend_side = "left",
#                      merge_legend = TRUE)
# 
# ## Two heatmaps sharing the rows, one next to the other
# ComplexHeatmap::draw(counts_heatmap_annotated_rows@heatmap + counts_heatmap_annotated_grouped@heatmap)

## ----heatmap_counts_anno_grab, eval=FALSE-------------------------------------
# heatmap_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(counts_heatmap_annotated@heatmap))
# 
# patchwork::wrap_plots(heatmap_grob, imputed_counts_heatmap@heatmap)

## ----heatmap_foldchanges, fig.width = 4, fig.height = 6-----------------------
FC_heatmap <-
  heatmap.contrasts(DEprot.analyses.object = dpo_analyses,
                    contrasts = c(1:2),
                    top.n = 20,
                    high.color = "#35978F",
                    low.color = "#BF812D",
                    mid.color = "white",
                    show.protein.names = TRUE,
                    use.uncorrected.pvalue = TRUE,
                    protein.names.pattern = "protein[.]")

FC_heatmap@heatmap

## ----upset_plot, message=F, warning=F-----------------------------------------
upset.plot <- plot.upset(DEprot.analyses.object = dpo_analyses,
                         contrast.subset = c(1,2),
                         title = "**My upset plot**",
                         sort.intersections = "cardinality",
                         sort.sets = "descending",
                         intersection.bar.color = "navy",
                         setsize.bar.color = "black",
                         show.counts = T,
                         height.ratio = 0.5,
                         width.ratio = 0.4,
                         use.uncorrected.pvalue = TRUE)

upset.plot  # or upset.plot@upset

## ----eval=FALSE---------------------------------------------------------------
# upset.plot@obs.matrix

## ----upset_tb, echo=FALSE-----------------------------------------------------
knitr::kable(upset.plot@obs.matrix[1:5,], row.names = F, caption = "**Upset observations matrix**")

## ----export_analyses, eval = F------------------------------------------------
# export.analyses(DEprot.analyses.object = dpo_analyses, output.folder = "./export")

## ----export_qc_report, eval = F-----------------------------------------------
# export.report(DEprot.object = dpo_analyses,
#               output.file = "QC_report.html",
#               report.title = "QC Report",
#               author.name = "Your Name",
#               protein.summary.group.column = "combined.id",
#               PCA.color.column = "combined.id",
#               PCA.shape.column = "replicate")

## ----run_saint, fig.width = 8, fig.height=5-----------------------------------
saint_deprot <-
  SAINTq(DEprot.object = DEprot::rime.dpo,
         metadata.column = "group",
         control = "LNCaP_TRIM33-5#MC-C2_FBS_IgG",
         bait = "LNCaP_TRIM33-5#MC-C2_FBS_AR",
         which.data = "imputed", # use 'raw' to be closer to the original SAINTq
         fold = 5)
saint_deprot

## ----show_saint_results, eval = F---------------------------------------------
# head(summary(saint_deprot))

## ----display_saint_tb_head, echo=FALSE----------------------------------------
knitr::kable(saint_deprot@scores$`LNCaP_TRIM33-5#MC-C2_FBS_AR`[1:5,], row.names = FALSE, caption = "**SAINTq results**")

## ----compare_saintq_to_deprot, fig.height=5, fig.width=6----------------------
sq_combo = dplyr::left_join(x = saint_deprot@scores$`LNCaP_TRIM33-5#MC-C2_FBS_AR`,
                            y = DEprot::rime.saintq, # original tool
                            by = "Prey")

ggpubr::ggscatter(data = sq_combo,
                  x = "AvgP.x",
                  y = "AvgP.y",
                  alpha = 0.5,
                  stroke = NA,
                  xlab = "DEprot-computed SAINT",
                  ylab = "SAINTq",
                  title = "SAINT score computation benchmarking") +
  geom_smooth(formula = y ~ x, method = "lm", color = "steelblue", fill = "steelblue") +
  ggpubr::stat_cor(method = "pearson", r.digits = 3) +
  geom_vline(xintercept = 0.95, color = 'gray', linetype = 2) +
  geom_hline(yintercept = 0.95, color = 'gray', linetype = 2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5),
        panel.border = element_rect(fill = NA, colour = "black"))

## ----power_estimate, fig.width = 10, fig.height = 4---------------------------
power.estimation <- estimate.power(DEprot.analyses.object = dpo_analyses,
                                   contrast = 1,
                                   sample.size.range = c(2, 30),
                                   target.power = 0.8)

power.estimation

## ----power_plots, fig.width = 12, fig.height = 4------------------------------
plot(power.estimation, nrow = 1)

## ----power_table, eval = FALSE------------------------------------------------
# head(power.estimation@power.table)

## ----show_power_table, echo=FALSE---------------------------------------------
knitr::kable(head(power.estimation@power.table), row.names = FALSE, caption = "**Power table results**")

## ----power_desiredFC, fig.width = 10, fig.height = 4--------------------------
power.FC1.5 <- estimate.power(DEprot.analyses.object = dpo_analyses,
                              contrast = 1,
                              desired.FC = 1.5,
                              sd.quantile = 0.75,
                              target.power = 0.8)

power.FC1.5

## ----power_complete_cases, eval = FALSE---------------------------------------
# complete.proteins <- rownames(dpo_norm@norm.counts)[rowSums(is.na(dpo_norm@norm.counts)) == 0]
# 
# dpo_complete <- filter.proteins(DEprot.object = dpo_norm,
#                                 proteins = complete.proteins,
#                                 mode = "keep")
# 
# dpo_complete <- diff.analyses.limma(DEprot.object = dpo_complete,
#                                     contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO")),
#                                     which.data = "normalized")
# 
# estimate.power(dpo_complete, contrast = 1)

## ----power4peaks_conversion, eval = FALSE-------------------------------------
# p4p.stats <- power4peaks::as.power4peaks(object = dpo_analyses,
#                                          contrast = 1)
# 
# p4p.power <- power4peaks::compute.power(power4peaks.stats = p4p.stats,
#                                         sample.size.range = c(2, 30),
#                                         power.threshold = 0.8)
# 
# p4p.power

## ----protein_summary, fig.height=3, fig.width=10------------------------------
protein.counts <-
  protein.summary(DEprot.object = dpo_analyses,
                  n.labels = "counts",
                  show.frequency = FALSE,
                  colors = c("gray", "steelblue4"),
                  title = "**# Proteins identified in each sample**")

protein.counts

## ----protein_summary_by_condition, fig.width=4, fig.height=3------------------
protein.counts.byCondition <- 
  protein.summary(DEprot.object = dpo_analyses,
                  group.column = "condition",
                  n.labels = "percentage",
                  show.frequency = TRUE,
                  x.label.angle = 0,
                  title = "**# Proteins identified per _Condition_**")

protein.counts.byCondition

## ----expression_boxplot, fig.width=10.5, fig.height=4-------------------------
### raw expression
protein.1733_raw <-
  expression.boxplot(DEprot.object = dpo,
                     protein.id = "protein.1733",
                     which.data = "imputed",
                     shape.column = "replicate",
                     group.by.metadata.column = "condition",
                     group.levels = c("6h.DMSO", "6h.10nM.E2", "FBS"),
                     scale.expression = FALSE,
                     x.label.angle = 90)


### scaled expression
protein.1733_scaled <-
  expression.boxplot(DEprot.object = dpo,
                     protein.id = "protein.1733",
                     which.data = "imputed",
                     shape.column = "replicate",
                     group.by.metadata.column = "condition",
                     group.levels = c("6h.DMSO", "6h.10nM.E2", "FBS"),
                     scale.expression = TRUE,
                     x.label.angle = 90)


### scaled expression + p-values
protein.1733_scaled_pairwise <-
  expression.boxplot(DEprot.object = dpo,
                     protein.id = "protein.1733",
                     which.data = "imputed",
                     shape.column = "replicate",
                     group.by.metadata.column = "condition",
                     group.levels = c("6h.DMSO", "6h.10nM.E2", "FBS"),
                     scale.expression = TRUE,
                     x.label.angle = 90,
                     pairwise.comparisons = TRUE,
                     pairwise.test.type = "wilcox",
                     pairwise.p.label = "p.value",
                     pairwise.p.decimals = 3,
                     pairwise.include.ns = FALSE)

patchwork::wrap_plots(protein.1733_raw,
                      protein.1733_scaled,
                      protein.1733_scaled_pairwise,
                      nrow = 1)

## ----keep_nuclear_prot, eval = F----------------------------------------------
# nucleus <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
#                                  keytype = "GOALL",
#                                  keys = "GO:0005634", #nucleus
#                                  columns = c("SYMBOL", "UNIPROT"))
# 
# dpo_analyses_nuclear <- filter.proteins(DEprot.object = dpo_analyses,
#                                         proteins = nucleus$SYMBOL,
#                                         mode = "keep")

## ----remove_cytoplasmic_proteins, eval = F------------------------------------
# cytoplasm <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
#                                    keytype = "GOALL",
#                                    keys = "GO:0005737", #cytoplasm
#                                    columns = c("SYMBOL", "UNIPROT"))
# 
# dpo_analyses_nuclear <- filter.proteins(DEprot.object = dpo_analyses,
#                                         proteins = cytoplasm$SYMBOL,
#                                         mode = "remove")

## ----filter_samples, eval = FALSE---------------------------------------------
# dpo.clean <- filter.samples(DEprot.object = dpo.imputed,
#                             samples = c("Sample_A", "Sample_B"),
#                             mode = "remove",
#                             verbose = TRUE)

## -----------------------------------------------------------------------------
dpo_analyses_fc1.5 <-
  reapply.thresholds(dpo_analyses,
                     linear.FC = 1.5,
                     p.adjusted = 0.05,
                     linear.FC.unresp.range = c(1/1.1, 1.1),
                     up.color = "indianred",
                     down.color = "steelblue",
                     unresponsive.color = "purple",
                     null.color = "gray")

summary(dpo_analyses_fc1.5)

## ----define_corum, eval = F---------------------------------------------------
# ## GeneSet Enrichment Analyses
# data("corum_v5.0", package = "DEprot")
# 
# corum_geneSet <-
#   corum_v5.0 %>%
#   dplyr::filter(organism == "Human") %>%
#   dplyr::rename(gs_name = complex.name,
#                 gene_symbol = protein.members) %>%
#   dplyr::select(gs_name, gene_symbol)
# 
# corum_geneSet

## ----print_corum_geneset, echo=FALSE------------------------------------------
data("corum_v5.0", package = "DEprot")

corum_geneSet <-
  corum_v5.0 %>%
  dplyr::filter(organism == "Human") %>%
  dplyr::rename(gs_name = complex.name,
                gene_symbol = protein.members) %>%
  dplyr::select(gs_name, gene_symbol)

knitr::kable(corum_geneSet[1:10,], row.names = FALSE, caption = "**CORUM protein complexes (v5.0)**")

## ----compare_ranking, fig.width=15, fig.height=10-----------------------------
compare.ranking(DEprot.analyses.object = dpo_analyses,
                contrast = 2)

## ----GSEA, eval = F-----------------------------------------------------------
# GSEA.results <-
#   geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                      contrast = 1,
#                      TERM2GENE = corum_geneSet,
#                      enrichment.type = "GSEA",
#                      gsea.rank.method = "foldchange", # or correlation
#                      gsub.pattern.prot.id = "_HUMAN|;.*",
#                      pvalueCutoff = 0.05,
#                      pAdjustMethod = "BH",
#                      dotplot.n = 10)

## ----ORA, eval = FALSE--------------------------------------------------------
# ORA.results <-
#   geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                      contrast = 1,
#                      TERM2GENE = corum_geneSet,
#                      enrichment.type = "ORA",
#                      gsub.pattern.prot.id = "_HUMAN|;.*",
#                      pvalueCutoff = 0.05,
#                      qvalueCutoff = 0.05,
#                      pAdjustMethod = "BH",
#                      diff.status.category = "6h.10nM.E2",
#                      dotplot.n = 10)

## ----simplify_enrichment, eval = FALSE----------------------------------------
# GSEA.results.simplified <- simplify.enrichment(GSEA.results)
# ORA.results.simplified <- simplify.enrichment(ORA.results)

## ----combine_enrichments, eval = FALSE----------------------------------------
# combined <-
#   combine.enrichments(enrichment.list = list(`E2 6h` = ORA.results.6h,
#                                              `E2 24h` = ORA.results.24h,
#                                              `E2 48h` = ORA.results.48h),
#                       dotplot.n = 5,
#                       padj.cutoff = 0.05,
#                       size.by = "FoldEnrichment",
#                       order.by = "discovery")
# 
# combined$dotplot

## ----divergent_enrichment, eval = FALSE---------------------------------------
# ORA.E2 <-
#   geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                      contrast = 1,
#                      TERM2GENE = corum_geneSet,
#                      enrichment.type = "ORA",
#                      gsub.pattern.prot.id = "_HUMAN|;.*",
#                      diff.status.category = "6h.10nM.E2")
# 
# ORA.DMSO <-
#   geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                      contrast = 1,
#                      TERM2GENE = corum_geneSet,
#                      enrichment.type = "ORA",
#                      gsub.pattern.prot.id = "_HUMAN|;.*",
#                      diff.status.category = "6h.DMSO")
# 
# 
# divergent <-
#   divergent.enrichment(enrichment.list = list(`6h E2` = ORA.E2,
#                                               `6h DMSO` = ORA.DMSO),
#                        value = "FoldEnrichment",
#                        top.n = 10,
#                        padj.cutoff = 0.05)
# 
# divergent$divergent.plot

## ----import_external, eval=F--------------------------------------------------
# dpo <- import.external(file = "report.pg_matrix.tsv",
#                        metadata = sample.config,
#                        source = "diann.matrix")

## ----import_diann_maxlfq, eval=F----------------------------------------------
# dpo <- import.external(file = "report.parquet",
#                        metadata = sample.config,
#                        source = "diann",
#                        summarization = "maxlfq",
#                        q.value = 0.01,
#                        pg.q.value = 0.01)

## ----import_install, eval=F---------------------------------------------------
# # install.missing: "ask" (default), "always" or "never"
# dpo <- import.external(file = "report.parquet",
#                        metadata = sample.config,
#                        install.missing = "ask")

## ----import_shortcut, eval=F--------------------------------------------------
# dpo <- read.fragpipe(file = "combined_protein.tsv",
#                      metadata = sample.config,
#                      quantity = "MaxLFQ Intensity")

## ----import_msstats, eval=F---------------------------------------------------
# # label-free
# summarized <- MSstats::dataProcess(raw)
# dpo <- import.msstats(object = summarized)
# 
# # isobaric labelling (TMT)
# summarized.tmt <- MSstatsTMT::proteinSummarization(input.pd)
# dpo <- import.msstats(object = summarized.tmt)

## ----import_msstats_metadata, eval=F------------------------------------------
# dpo <- import.msstats(object = summarized.tmt)
# get.metadata(dpo)

## ----import_msstats_harmonize, eval=F-----------------------------------------
# dpo <- harmonize.batches(DEprot.object = dpo,
#                          batch.column = "Mixture")

## ----export_external, eval=F--------------------------------------------------
# se <- export.external(DEprot.object = dpo,
#                       format = "SummarizedExperiment")

## ----export_counts_type, eval=F-----------------------------------------------
# # imputed counts as primary assay, all the others still accessible
# se <- export.external(DEprot.object = dpo, counts.type = "imputed")
# 
# # only two matrices, normalized as primary
# se <- export.external(DEprot.object = dpo, assays = c("raw", "normalized"))

## ----export_results, eval=F---------------------------------------------------
# se <- export.external(DEprot.object = dpo.analyses,
#                       add.results = TRUE,
#                       contrast.subset = c(1, 3))
# 
# SummarizedExperiment::rowData(se)

## ----export_metadata, eval=F--------------------------------------------------
# se <- export.external(DEprot.object = dpo, keep.object = TRUE)
# 
# S4Vectors::metadata(se)$imputation.method
# dpo <- S4Vectors::metadata(se)$DEprot.object

## ----export_shortcut, eval=F--------------------------------------------------
# se <- as.SummarizedExperiment(dpo.analyses)
# qf <- as.QFeatures(dpo, assay.name = "proteins")

## ----session_info-------------------------------------------------------------
sessionInfo()

