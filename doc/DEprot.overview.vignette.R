## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">", dev = "svg",
                      warning = F, message = F, fig.align = "center",
                      rows.print=12)
#knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
#options(tibble.print_min = 4L, tibble.print_max = 4L)


# Load libraries required
require(DEprot)
require(dplyr)

## ----citation, message=FALSE, warning=FALSE-----------------------------------
citation("DEprot")

## ----read_metadata, eval=F----------------------------------------------------
#  # Metadata
#  data("sample.config", package = "DEprot")
#  sample.config

## ----print_metadata, echo=FALSE-----------------------------------------------
data("sample.config", package = "DEprot")
knitr::kable(sample.config, row.names = F, caption = "**Sample metadata table**")

## ----read_counts, eval=F------------------------------------------------------
#  # log2(LFQ) values (not imputed)
#  data("unimputed.counts", package = "DEprot")
#  head(unimputed.counts[,1:6])

## ----print_counts, echo=FALSE-------------------------------------------------
data("unimputed.counts", package = "DEprot")
knitr::kable(data.frame(head(unimputed.counts[,1:6])), row.names = T, caption = "**Unimputed log2(LFQ) values**")

## ----make_dpo-----------------------------------------------------------------
dpo <- load.counts(counts = unimputed.counts,
                   metadata = sample.config,
                   log.base = 2,
                   imputation = NA,
                   normalization.method = NA,
                   column.id = "column.id")

dpo

## ----rename_samples-----------------------------------------------------------
dpo <- rename.samples(DEprot.object = dpo,
                      metadata.column = "sample.id")

get.metadata(dpo)

## ----display_rename_samples---------------------------------------------------
head(dpo@raw.counts[,1:6])

## ----display_boxplot_raw, fig.width=5, fig.cap='<div style="text-align: justify"> <font size="-0.5"> **Violin/boxplot LFQ intensities of unnormalized data** <br> Boxplots display the quantiles of the LFQ intensities, while red and blue dahsed lines correspond to maximum and minimum LFQ value for each sample. </font> </div> <br>'----
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
patchwork::wrap_plots(dpo@boxplot.raw, dpo@boxplot.norm, nrow = 1)

## ----eval = F-----------------------------------------------------------------
#  ## Adding batch column to the metadata table
#  dpo@metadata$batch = rep(c("batch_1","batch_2"), each = 6)
#  
#  get.metadata(dpo)

## ----echo=FALSE---------------------------------------------------------------
new_meta = dpo@metadata
new_meta$batch = rep(c("batch_1","batch_2"), each = 6)
new_meta

## ----eval = F-----------------------------------------------------------------
#  ## batch correction
#  dpo <- harmonize.batches(DEprot.object = dpo,
#                           batch.column = "batch",
#                           cores = 1)

## ----randomize_data-----------------------------------------------------------
dpo <- randomize.missing.values(DEprot.object = dpo,
                                group.column = "combined.id",
                                tail.percentage = 5)

## ----imputation_example, eval = F---------------------------------------------
#  ## Without parallelization
#  dpo <- impute.counts(DEprot.object = dpo,
#                       method = "missForest",
#                       use.normalized.data = TRUE,
#                       missForest.max.iterations = 100,
#                       missForest.variable.wise.OOBerror = T)
#  
#  
#  ## With parallelization
#  dpo <- impute.counts(DEprot.object = dpo,
#                       method = "missForest",
#                       use.normalized.data = TRUE,
#                       missForest.max.iterations = 100,
#                       missForest.variable.wise.OOBerror = TRUE,
#                       missForest.cores = 10,
#                       missForest.parallel.mode = "variables")
#  
#  dpo
#  
#  dpo@imputation$OOBerror
#  
#  data.frame(dpo@imputation[-3])

## ----load_imputation, echo=FALSE----------------------------------------------
dpo = readRDS(url("https://data.cyverse.org/dav-anon/iplant/home/sgregoricchio/deprot/dpo.imputed.rds"))
dpo

## ----load_imputation2, echo=FALSE---------------------------------------------
error = dpo@imputation$OOBerror
names(error) = colnames(dpo@imputed.counts)
error

knitr::kable(data.frame(dpo@imputation[-3]), row.names = F)

## ----plot_imputed, fig.width=9, fig.align='center'----------------------------
patchwork::wrap_plots(dpo@boxplot.raw, dpo@boxplot.norm, dpo@boxplot.imputed, nrow = 1)

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
                  plot.zero.lines = F) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  theme(legend.position = "none")


PC_2.3 <-
  plot.PC.scatter(DEprot.PCA.object = PCA,
                  PC.x = 2,
                  PC.y = 3,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.lines = TRUE)

patchwork::wrap_plots(PC_1.2, PC_2.3, nrow = 1)

## -----------------------------------------------------------------------------
plot.PC.scatter.123(DEprot.PCA.object = PCA,
                    color.column = "condition",
                    shape.column = "replicate",
                    label.column = "replicate",
                    dot.colors = c("6h.10nM.E2" = "indianred",
                                   "6h.DMSO" = "steelblue",
                                   "FBS" =  "forestgreen"),
                    plot.zero.line.y.12 = TRUE,
                    plot.zero.line.x.12 = FALSE,
                    plot.zero.line.y.23 = TRUE,
                    plot.zero.line.x.23 = TRUE)

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
                  plot.zero.lines = F) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  theme(legend.position = "none")
  
  
PC.fbs.e2_2.3 <-
  plot.PC.scatter(DEprot.PCA.object = PCA.fbs.e2,
                  PC.x = 2,
                  PC.y = 3,
                  color.column = "condition",
                  shape.column = "replicate",
                  label.column = NULL,
                  plot.zero.lines = T)

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

## ----compute_diff_exp_examples_limma, eval=F----------------------------------
#  ## Unpaired test
#  dpo_analyses <- diff.analyses.limma(DEprot.object = dpo,
#                                      contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                           c("condition", "6h.10nM.E2", "FBS")),
#                                      linear.FC.th = 2,
#                                      padj.th = 0.05,
#                                      padj.method = "BH",
#                                      fitting.method = "ls",
#                                      which.data = "imputed")
#  
#  ## Paired test
#  dpo_analyses <- diff.analyses.limma(DEprot.object = dpo,
#                                      contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                           c("condition", "6h.10nM.E2", "FBS")),
#                                      replicate.column = "replicate",
#                                      include.rep.model = TRUE,
#                                      linear.FC.th = 2,
#                                      padj.th = 0.05,
#                                      padj.method = "BH",
#                                      fitting.method = "ls",
#                                      which.data = "imputed")

## ----compute_diff_exp_examples_Ttest, eval=F----------------------------------
#  ## Unpaired test
#  dpo_analyses <- diff.analyses(DEprot.object = dpo,
#                                contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                     c("condition", "6h.10nM.E2", "FBS")),
#                                linear.FC.th = 2,
#                                padj.th = 0.05,
#                                padj.method = "bonferroni",
#                                stat.test = "t.test",
#                                which.data = "imputed")
#  
#  ## Paired test
#  dpo_analyses <- diff.analyses(DEprot.object = dpo,
#                                contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
#                                                     c("condition", "6h.10nM.E2", "FBS")),
#                                replicate.column = "replicate",
#                                paired.test = TRUE,
#                                linear.FC.th = 2,
#                                padj.th = 0.05,
#                                padj.method = "bonferroni",
#                                stat.test = "t.test",
#                                which.data = "imputed")
#  
#  dpo_analyses

## ----compute_diff_exp_paired, echo=FALSE--------------------------------------
## Paired test
dpo_analyses <- diff.analyses(DEprot.object = dpo,
                              contrast.list = list(c("condition", "6h.10nM.E2", "6h.DMSO"),
                                                   c("condition", "6h.10nM.E2", "FBS")),
                              replicate.column = "replicate",
                              paired.test = TRUE,
                              linear.FC.th = 2,
                              padj.th = 0.05,
                              padj.method = "bonferroni",
                              stat.test = "t.test",
                              which.data = "imputed")

dpo_analyses

## ----analyses_summary, eval=F-------------------------------------------------
#  diff.analyses_summary = summary(dpo)

## ----get_results, eval = F----------------------------------------------------
#  ## Direct access
#  results = dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$results
#  
#  ## Function
#  results = get.results(dpo_analyses, contrast = 1)
#  
#  head(results)

## ----get_results2, echo=FALSE-------------------------------------------------
knitr::kable(get.results(dpo_analyses, contrast = 1)[1:6,], row.names = F)

## ----DE_PCA_scatters, fig.width=8, eval=F-------------------------------------
#  dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.plots

## ----DE_PCA_scatters_replotting, echo=FALSE, fig.width=8----------------------
scatter_PC12 = plot.PC.scatter(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data, 1,2, color.column = "condition") + theme(legend.position = "none")
scatter_PC23 = plot.PC.scatter(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data, 2,3, color.column = "condition")
cumulative = plot.PC.cumulative(dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$PCA.data)

scatters = cowplot::plot_grid(scatter_PC12, scatter_PC23, nrow = 1, align = "hv", axis = "tblr")
cowplot::plot_grid(scatters, cumulative, ncol = 1, axis = "tblr")

## ----DE_correlations, fig.width=10--------------------------------------------
dpo_analyses@analyses.result.list$condition_6h.10nM.E2.vs.6h.DMSO$correlations

## ----DE_volcano_MA, fig.width=9-----------------------------------------------
volcano = plot.volcano(dpo_analyses, contrast = 1, use.uncorrected.pvalue = TRUE)
MAplot = plot.MA(dpo_analyses, contrast = 1, use.uncorrected.pvalue = TRUE)

patchwork::wrap_plots(volcano, MAplot)

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

## ----heatmap_foldchanges, fig.width = 4, fig.height = 6-----------------------
FC_heatmap <-
  heatmap.contrasts(DEprot.analyses.object = dpo_analyses,
                    contrasts = c(1:2),
                    top.n = 20,
                    high.color = "#35978F",
                    low.color = "#BF812D",
                    mid.color = "white",
                    show.protein.names = TRUE,
                    use.uncorrected.pvalue = TRUE)

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
#  upset.plot@obs.matrix

## ----upset_tb, echo=FALSE-----------------------------------------------------
knitr::kable(upset.plot@obs.matrix[1:5,], row.names = F, caption = "**Upset observations matrix**")

## ----protein_summary, fig.height=3, fig.width=10------------------------------
protein.counts <-
  protein.summary(DEprot.object = dpo_analyses,
                  n.lables = "counts",
                  show.frequency = F,
                  colors = c("gray", "steelblue4"),
                  title = "**# Proteins identified in each sample**")

protein.counts

## ----protein_summary_by_condition, fig.width=4, fig.height=3------------------
protein.counts.byCondition <- 
  protein.summary(DEprot.object = dpo_analyses,
                  group.column = "condition",
                  n.lables = "percentage",
                  show.frequency = T,
                  x.label.angle = 0,
                  title = "**# Proteins identified per _Condition_**")

protein.counts.byCondition

## ----keep_nuclear_prot, eval = F----------------------------------------------
#  nucleus <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
#                                   keytype = "GOALL",
#                                   keys = "GO:0005634", #nucleus
#                                   columns = c("SYMBOL", "UNIPROT"))
#  
#  dpo_analyses_nuclear <- filter.proteins(DEprot.object = dpo_analyses,
#                                          proteins = nucleus$SYMBOL,
#                                          mode = "keep")

## ----remove_cytoplasmic_proteins, eval = F------------------------------------
#  cytoplasm <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
#                                     keytype = "GOALL",
#                                     keys = "GO:0005737", #cytoplasm
#                                     columns = c("SYMBOL", "UNIPROT"))
#  
#  dpo_analyses_nuclear <- filter.proteins(DEprot.object = dpo_analyses,
#                                          proteins = cytoplasm$SYMBOL,
#                                          mode = "remove")

## ----define_corum, eval = F---------------------------------------------------
#  ## GeneSet Enrichment Analyses
#  data("corum_v5.0", package = "DEprot")
#  
#  corum_geneSet =
#    corum_v5.0 %>%
#    dplyr::filter(organism == "Human") %>%
#    dplyr::rename(gs_name = complex.name,
#                  gene_symbol = protein.members) %>%
#    dplyr::select(gs_name, gene_symbol)
#  
#  corum_geneSet

## ----print_corum_geneset, echo=FALSE------------------------------------------
## GeneSet Enrichment Analyses
data("corum_v5.0", package = "DEprot")

corum_geneSet =
  corum_v5.0 %>%
  dplyr::filter(organism == "Human") %>%
  dplyr::rename(gs_name = complex.name,
                gene_symbol = protein.members) %>%
  dplyr::select(gs_name, gene_symbol)

knitr::kable(corum_geneSet[1:10,], row.names = F, caption = "**CORUM protein complexes (v5.0)**")

## ----enrichment_analyses, eval = F--------------------------------------------
#  ### GSEA
#  GSEA.results =
#    geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                       contrast = 1,
#                       TERM2GENE = corum_geneSet,
#                       enrichment.type = "GSEA",
#                       gsub.pattern.prot.id = "_HUMAN|;.*",
#                       pvalueCutoff = 0.05,
#                       pAdjustMethod = "BH",
#                       dotplot.n = 10)
#  
#  
#  
#  ### OverRepresentation Analyses (ORA)
#  ORA.results =
#    geneset.enrichment(DEprot.analyses.object = dpo_analyses,
#                       contrast = 1,
#                       TERM2GENE = corum_geneSet,
#                       enrichment.type = "GSEA",
#                       gsub.pattern.prot.id = "_HUMAN|;.*",
#                       pvalueCutoff = 0.05,
#                       qvalueCutoff = 0.05,
#                       pAdjustMethod = "BH",
#                       diff.status.category = "6h.10nM.E2",
#                       dotplot.n = 10)

