#' @title DEprot class
#' @name DEprot
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{metadata}:}{The data.frame corresponding to the metadata table describing the samples. Class: \code{"ANY"}.}
#'    \item{\code{raw.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the raw counts. Class: \code{"ANY"}.}
#'    \item{\code{norm.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the normalized counts. Class: \code{"ANY"}.}
#'    \item{\code{imputed.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the imputed counts. Class: \code{"ANY"}.}
#'    \item{\code{log.base}:}{Numeric value indicating the base of the logarithm expressing the values in the loaded data. Class: \code{"ANY"}.}
#'    \item{\code{log.transformed}:}{Logical value indicating whether the data are log-transformed or not. Class: \code{"logical"}.}
#'    \item{\code{imputed}:}{Logical value indicating whether the data are imputed. Class: \code{"logical"}.}
#'    \item{\code{imputation}:}{String (or any other class) value indicating the imputation method. Class: \code{"ANY"}.}
#'    \item{\code{normalized}:}{Logical value indicating whether the data are normalized. Class: \code{"logical"}.}
#'    \item{\code{normalization.method}:}{String (or any other class) value indicating the normalization method. Class: \code{"ANY"}. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.raw}:}{Ggplot object showing the distribution of the raw values per sample. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.norm}:}{Ggplot object showing the distribution of the normalized values per sample. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.imputed}:}{Ggplot object showing the distribution of the imputed values per sample. Class: \code{"ANY"}.}
#'    \item{\code{analyses.result.list}:}{For this type of object the value is \code{NULL}. Class: \code{"ANY"}.}
#'    \item{\code{contrasts}:}{For this type of object the value is \code{NULL}. Class: \code{"ANY"}.}
#'    \item{\code{differential.analyses.params}:}{For this type of object the value is \code{NULL}. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot

setClass(Class = "DEprot",
         slots = list(metadata = "ANY",
                      raw.counts = "ANY",
                      norm.counts = "ANY",
                      imputed.counts = "ANY",
                      log.base = "numeric",
                      log.transformed = "logical",
                      imputed = "logical",
                      imputation = "ANY",
                      normalized = "logical",
                      normalization.method = "ANY",
                      boxplot.raw = "ANY",
                      boxplot.norm = "ANY",
                      boxplot.imputed = "ANY",
                      analyses.result.list = "ANY",
                      contrasts = "ANY",
                      differential.analyses.params = "ANY"))






#' @title DEprot.analyses class
#' @name DEprot.analyses
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{metadata}:}{The data.frame corresponding to the metadata table describing the samples. Class: \code{"ANY"}.}
#'    \item{\code{raw.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the raw counts. Class: \code{"ANY"}.}
#'    \item{\code{norm.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the normalized counts. Class: \code{"ANY"}.}
#'    \item{\code{imputed.counts}:}{Numeric matrix (rows: proteins, columns: samples) of the imputed counts. Class: \code{"ANY"}.}
#'    \item{\code{log.base}:}{Numeric value indicating the base of the logarithm expressing the values in the loaded data. Class: \code{"ANY"}.}
#'    \item{\code{log.transformed}:}{Logical value indicating whether the data are log-transformed or not. Class: \code{"logical"}.}
#'    \item{\code{imputed}:}{Logical value indicating whether the data are imputed. Class: \code{"logical"}.}
#'    \item{\code{imputation}:}{String (or any other class) value indicating the imputation method. Class: \code{"ANY"}.}
#'    \item{\code{normalized}:}{Logical value indicating whether the data are normalized. Class: \code{"logical"}.}
#'    \item{\code{normalization.method}:}{String (or any other class) value indicating the normalization method. Class: \code{"ANY"}. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.raw}:}{Ggplot object showing the distribution of the raw values per sample. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.norm}:}{Ggplot object showing the distribution of the normalized values per sample. Class: \code{"ANY"}.}
#'    \item{\code{boxplot.imputed}:}{Ggplot object showing the distribution of the imputed values per sample. Class: \code{"ANY"}.}
#'    \item{\code{analyses.result.list}:}{List containing the differential results for each contrast. Class: \code{"ANY"}. The list contains the following elements:
#'      \itemize{
#'        \item{\code{results}: }{a data.frame containing the results of the analyses; includes average expression of each group, basemean, foldchange, pvalue and p.adj, differential.status}
#'        \item{\code{n.diff}: }{a summary table showing the number of proteins in each differential expression status (up/down/unresponsive, null)}
#'        \item{\code{PCA.data}: }{output of \link{perform.PCA} for the subset of samples analyzed in a specific contrast}
#'        \item{\code{PCA.plots}: }{combination of 3 plots: scatter PC1-vs-PC2, scatter PC2-vs-PC3, and cumulative bar plot}
#'        \item{\code{correlations}: }{combination of Pearson and Spearman correlation heatmaps (obtained by \link{plot.correlation.heatmap}) for the subset of samples analyzed in a specific contrast}
#'        \item{\code{volcano}: }{volcano plot showing the log2(FoldChange) x -log10(p.adjusted) of differential expression results; it can be regenerated using \link{plot.volcano}}
#'        \item{\code{MA.plot}: }{MA-plot showing the log2(basemean) x log2(FoldChange) of differential expression results; it can be regenerated using \link{plot.MA}}}}
#'    \item{\code{contrasts}:}{List of contrasts. each contrast is a vector indicating, in the order: metadata.table.column - groupA - groupB; (groupA / group B). Class: \code{"ANY"}.}
#'    \item{\code{differential.analyses.params}:}{List of parameters used to run the differential analyses (fold change thresholds, p-value threshold, p-adjustement method, etc.). Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.analyses

setClass(Class = "DEprot.analyses",
         slots = list(metadata = "ANY",
                      raw.counts = "ANY",
                      norm.counts = "ANY",
                      imputed.counts = "ANY",
                      log.base = "numeric",
                      log.transformed = "logical",
                      imputed = "logical",
                      imputation = "ANY",
                      normalized = "logical",
                      normalization.method = "ANY",
                      boxplot.raw = "ANY",
                      boxplot.norm = "ANY",
                      boxplot.imputed = "ANY",
                      analyses.result.list = "ANY",
                      contrasts = "ANY",
                      differential.analyses.params = "ANY"))







#' @title DEprot.PCA class
#' @name DEprot.PCA
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{PCA.metadata}:}{metadata of the samples used in the PCA (subset of the original \code{DEprot@@metadata}. Class: \code{"ANY"}.}
#'    \item{\code{sample.subset}:}{vector containing the list of samples analyzed. Class: \code{"ANY"}.}
#'    \item{\code{data.used}:}{vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.}
#'    \item{\code{prcomp}:}{object of class \code{prcomp} (or output from \code{pcaMethods::pca}, method = "nipals") corresponding to the full PCA output. Class: \code{"ANY"}.}
#'    \item{\code{PCs}:}{data.frame combining the PC scores and the metadata table, useful for replotting. Class: \code{"ANY"}.}
#'    \item{\code{importance}:}{statistical summary table for the PCA analyses per each PC. Class: \code{"ANY"}.}
#'    \item{\code{cumulative.PC.plot}:}{ggplot object corresponding to out put of \code{plot.PC.cumulative} for this object. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.PCA

setClass(Class = "DEprot.PCA",
         slots = list(PCA.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      prcomp = "ANY",
                      PCs = "ANY",
                      importance = "ANY",
                      cumulative.PC.plot = "ANY"))






#' @title DEprot.correlation class
#' @name DEprot.correlation
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{heatmap}:}{ggplot object corresponding to the correlation heatmap. Class: \code{"ANY"}.}
#'    \item{\code{corr.metadata}:}{metadata of the samples used in the correlation (subset of the original \code{DEprot@@metadata}. Class: \code{"ANY"}.}
#'    \item{\code{sample.subset}:}{vector containing the list of samples analyzed. Class: \code{"ANY"}.}
#'    \item{\code{data.used}:}{vector indicating the type of counts used (imputed, normalized, raw). Class: \code{"ANY"}.}
#'    \item{\code{corr.matrix}:}{the correlation matrix on which the heatmap is base on. Class: \code{"ANY"}.}
#'    \item{\code{distance}:}{object of class \code{dist} corresponding to the output of \code{as.dist(1 - correlation.matrix)}. Class: \code{"ANY"}.}
#'    \item{\code{cluster}:}{\code{hclust} object generated by \code{hclust(d = as.dist(1 - correlation.matrix), method = clustering.method)}. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.correlation

setClass(Class = "DEprot.correlation",
         slots = list(heatmap = "ANY",
                      corr.metadata = "ANY",
                      sample.subset = "ANY",
                      data.used = "ANY",
                      corr.matrix = "ANY",
                      distance = "ANY",
                      cluster = "ANY"))



#' @title DEprot.upset class
#' @name DEprot.upset
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{upset}:}{Ggplot object corresponding to upset plot displaying the overlaps of differential proteins between cointrasts. Class: \code{"ANY"}.}
#'    \item{\code{obs.matrix}:}{Logical matrix indicating all the proteins that are differentially expressed at least in a contrast (rows). Columns indicate a specific contrast. The logical values indicate whether a protein is found differential in a specific contrast (column). Therefore, this table can be used to extract the proteins included in a specific overlap. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.upset

setClass(Class = "DEprot.upset",
         slots = list(upset = "ANY",
                      obs.matrix = "ANY"))



#' @title DEprot.contrast.heatmap class
#' @name DEprot.contrast.heatmap
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{heatmap}:}{Ggplot object corresponding to any heatmap generate by either \link{heatmap.contrasts}. Class: \code{"ANY"}.}
#'    \item{\code{cluster}:}{The \code{hclust} object of the rows (proteins). Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.contrast.heatmap

setClass(Class = "DEprot.contrast.heatmap",
         slots = list(heatmap = "ANY",
                      cluster = "ANY"))


#' @title DEprot.counts.heatmap class
#' @name DEprot.counts.heatmap
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{heatmap}:}{Ggplot object corresponding to any heatmap generate by either \link{heatmap.counts}. Class: \code{"ANY"}.}
#'    \item{\code{row.cluster}:}{The \code{hclust} object of the rows (proteins). Class: \code{"ANY"}.}
#'    \item{\code{column.cluster}:}{The \code{hclust} object of the columns (samples). Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.counts.heatmap

setClass(Class = "DEprot.counts.heatmap",
         slots = list(heatmap = "ANY",
                      row.cluster = "ANY",
                      column.cluster = "ANY"))



#' @title DEprot.enrichResult class
#' @name DEprot.enrichResult
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{enrichment.discovery}:}{the direct output from \code{clusterProfiler::GSEA} or \code{clusterProfiler::enricher} (GSEA and ORA, respectively). Class: \code{"ANY"}.}
#'    \item{\code{protein.network}:}{a string plot showing protein networks (\code{clusterProfiler::cnetplot}). Class: \code{"ANY"}.}
#'    \item{\code{pathway.network}:}{a list with clusters and string plot showing pathway/set networks (\code{aPEAR::enrichmentNetwork}). Class: \code{"ANY"}.}
#'    \item{\code{NES.plot}:}{(GSEA only) a bar plot showing the NES scores for each significantly enriched geneSet. Class: \code{"ANY"}.}
#'    \item{\code{dotplot_gene.ratio}:}{a dotplot showing the geneRatios for each significantly enriched geneSet. Class: \code{"ANY"}.}
#'    \item{\code{dotplot_fold.enrichment}:}{(ORA only) a dotplot showing the foldEnrichment for each significantly enriched geneSet. Class: \code{"ANY"}.}
#'    \item{\code{parameters}:}{a list containing the parameters used to run the analyses. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.enrichResult

setClass(Class = "DEprot.enrichResult",
         slots = list(enrichment.discovery = "ANY",
                      protein.network = "ANY",
                      pathway.network = "ANY",
                      NES.plot = "ANY",
                      dotplot_gene.ratio = "ANY",
                      dotplot_fold.enrichment = "ANY",
                      parameters = "ANY",
                      affinity.propagation = "ANY"))




#' @title DEprot.pvalues class
#' @name DEprot.pvalues
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{pvalue.distribution}:}{ggplot object depicting the histogram of the distribution of the p-values of the differential expression test of a specific contrast. Class: \code{"ANY"}.}
#'    \item{\code{padjusted.distribution}:}{ggplot object depicting the histogram of the distribution of the adjusted p-values of the differential expression test of a specific contrast. Class: \code{"ANY"}.}
#'    \item{\code{pvalue.rank}:}{ggplot object depicting the curve of the ranked p-values. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.pvalues

setClass(Class = "DEprot.pvalues",
         slots = list(pvalue.distribution = "ANY",
                      padjusted.distribution = "ANY",
                      pvalue.rank = "ANY"))




#' @title DEprot.normality class
#' @name DEprot.normality
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{norm.statement}:}{Logical value indicating whether the samples are normally distributed (\code{TRUE}). Class: \code{"ANY"}.}
#'    \item{\code{norm.AD.tests}:}{List of Anderson-Darling normality test results for each sample. Class: \code{"ANY"}.}
#'    \item{\code{qqplots}:}{List of ggplots objects depicting the Q-Q plots for the Anderson-Darling normality test. Class: \code{"ANY"}.}
#'    \item{\code{densities}:}{List of ggplots objects depicting the destiny distribution of the intensities overlapped to a theoretical normal distribution. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.normality

setClass(Class = "DEprot.normality",
         slots = list(norm.statement = "ANY",
                      norm.AD.tests = "ANY",
                      qqplots = "ANY",
                      densities = "ANY"))





#' @title DEprot.RMSE class
#' @name DEprot.RMSE
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{original.DEprot.object}:}{Object of class \code{DEprot} used to compute the RMSE. Class: \code{"ANY"}.}
#'    \item{\code{percentage.test}:}{Percentage of the total proteins that should be used to perform the comparisons. Class: \code{"ANY"}.}
#'    \item{\code{seed}:}{Seed used for the randomization. Class: \code{"ANY"}.}
#'    \item{\code{fraction.missing.values}:}{Fraction of missing values in the original table. Class: \code{"ANY"}.}
#'    \item{\code{test.dataset}:}{Subset of the original table used for the comparisons. Class: \code{"ANY"}.}
#'    \item{\code{imputed.objects}:}{List of the output of \link{impute.counts} (class \code{DEprot}) using the different imputation methods. Class: \code{"ANY"}.}
#'    \item{\code{RMSE.tables}:}{List of data.frames, one per tested imputation method, containing the coordinates of the imputed values and including the following columns (Class: \code{"ANY"}.):
#'      \itemize{
#'        \item{\code{row.id}: }{id of the row (protein)}
#'        \item{\code{col.id}: }{id of the column (sample)}
#'        \item{\code{expected.values}: }{value measured in the experiment}
#'        \item{\code{imputation.method}: }{id of the method used for the imputation}
#'        \item{\code{imputed.values}: }{value imputed by \code{DEprot}}
#'        \item{\code{residuals}: }{difference of the the values, \code{imputed - expected}}
#'        \item{\code{sq.residuals}: }{the squared value of the residuals}}}
#'    \item{\code{contrasts}:}{List of contrasts. each contrast is a vector indicating, in the order: metadata.table.column - groupA - groupB; (groupA / group B). Class: \code{"ANY"}.}
#'    \item{\code{RMSE.scores}:}{. Class: \code{"ANY"}.}
#'    \item{\code{correlation.plots}:}{. Class: \code{"ANY"}.}
#'    \item{\code{density.residuals}:}{. Class: \code{"ANY"}.}
#'  }
#'
#' @exportClass DEprot.RMSE

setClass(Class = "DEprot.RMSE",
         slots = list(original.DEprot.object = "ANY",
                      percentage.test = "ANY",
                      seed = "ANY",
                      fraction.missing.values = "ANY",
                      test.dataset = "ANY",
                      imputed.objects = "ANY",
                      RMSE.tables = "ANY",
                      RMSE.scores = "ANY",
                      correlation.plots = "ANY",
                      density.residuals = "ANY"))


################# METHODS #################

#' @title DEprot show-method
#' @export
setMethod(f = "show",
          signature = "DEprot",
          definition =
            function(object) {
              if (!is.null(object@raw.counts)) {
                tb = object@raw.counts
                tb.type.raw = "raw"
              } else {tb.type.raw = NULL}

              if (!is.null(object@norm.counts)) {
                tb = object@norm.counts
                tb.type.norm = "normalized"
              } else {tb.type.norm = NULL}

              if (!is.null(object@imputed.counts)) {
                tb = object@imputed.counts
                tb.type.imputed = "imputed"
              } else {tb.type.imputed = NULL}

              cnt.avilable = c(tb.type.raw, tb.type.norm, tb.type.imputed)
              cnt.avilable = cnt.avilable[!is.null(cnt.avilable)]


              cat("DEprot object:")
              cat("\n           Samples: ", ncol(tb))
              cat("\n          Proteins: ", nrow(tb))
              cat("\n  Counts available: ", paste(cnt.avilable, collapse = ", "))
              cat("\nLog transformation: ", ifelse(is.na(object@log.base), yes = "none (linear)", no = paste0("log", object@log.base)))
              cat("\n  Metadata columns: ", paste(colnames(object@metadata), collapse = ", "), "\n")
            }#end definition
) #end method



#' @title DEprot.analyses show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.analyses",
          definition =
            function(object) {

              ## Summarise differential expression results
              recap = data.frame()

              for (i in 1:length(object@analyses.result.list)) {
                recap = rbind(recap,
                              cbind(data.frame(contrast.id = rep(paste0(object@contrasts[[i]]$metadata.column, ": ", object@contrasts[[i]]$var.1, " vs ", object@contrasts[[i]]$var.2), 4),
                                               group.factor = rep(object@contrasts[[i]]$metadata.column, 4),
                                               group1 = rep(object@contrasts[[i]]$var.1, 4),
                                               group2 = rep(object@contrasts[[i]]$var.2, 4),
                                               paired.test = object@contrasts[[i]]$paired.test),
                                    object@analyses.result.list[[i]]$n.diff))

              }

              cat("DEprot.analyses object:")
              cat("\n          Counts used: ", object@differential.analyses.params$counts.used)
              cat("\nFold Change threshold: ", object@differential.analyses.params$linear.FC.th, "(linear)")
              cat("\nFC unresponsive range:  [", object@differential.analyses.params$linear.FC.unresp.range[1],",",object@differential.analyses.params$linear.FC.unresp.range[2],"] (linear)", sep = "")
              cat("\n       padj threshold: ", object@differential.analyses.params$padj.th, "(linear)")
              cat("\n          padj method: ", object@differential.analyses.params$padj.method)
              cat("\n")
              cat("\n")
              cat("\nDifferential results summary:\n")

              print(recap)
            }#end definition
) #end method



#' @title DEprot.analyses summary-method
#' @export
setMethod(f = "summary",
          signature = "DEprot.analyses",
          definition =
            function(object) {

              ## Summarise differential expression results
              recap = data.frame()

              for (i in 1:length(object@analyses.result.list)) {
                recap = rbind(recap,
                              cbind(data.frame(contrast.id = rep(paste0(object@contrasts[[i]]$metadata.column, ": ", object@contrasts[[i]]$var.1, " vs ", object@contrasts[[i]]$var.2), 4),
                                               group.factor = rep(object@contrasts[[i]]$metadata.column, 4),
                                               group1 = rep(object@contrasts[[i]]$var.1, 4),
                                               group2 = rep(object@contrasts[[i]]$var.2, 4),
                                               paired.test = object@contrasts[[i]]$paired.test),
                                    object@analyses.result.list[[i]]$n.diff))

              }
              recap
            }#end definition
) #end method



#' @title DEprot.PCA show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.PCA",
          definition =
            function(object) {
              cat("DEprot.PCA object:")
              cat("\n  Samples analyzed: ", object@sample.subset)
              cat("\n         Data used: ", paste(object@data.used, "(log2)\n"))
              cat("\n")
              cat("\n")
              cat("PCs vectors:\n")
              print(object@PCs)})



#' @title DEprot.correlation show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.correlation",
          definition = function(object) {print(object@heatmap)})



#' @title DEprot.upset show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.upset",
          definition = function(object) {print(object@upset)})


#' @title DEprot.cluster.heatmap show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.contrast.heatmap",
          definition = function(object) {print(object@heatmap)})


#' @title DEprot.counts.heatmap show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.counts.heatmap",
          definition = function(object) {print(object@heatmap)})



#' @title DEprot.enrichResult show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.enrichResult",
          definition = function(object) {print(object@enrichment.discovery@result)})


#' @title DEprot.pvalues show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.pvalues",
          definition =
            function(object) {
              require(patchwork, quietly = TRUE)
              plot = (object@pvalue.distribution / object@padjusted.distribution) | object@pvalue.rank
              print(plot)
            })


#' @title DEprot.normality show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.normality",
          definition =
            function(object) {
              if (object@norm.statement == TRUE) {
                message("All samples display a normal distribution.")

              } else {
                normality = sapply(X = object@norm.AD.tests, FUN = function(x){x$p.value < p.threshold}, USE.NAMES = TRUE)
                message(paste0("The following samples do not display a normal distribution: ",
                               paste0(names(normality)[isFALSE(normality)], collapse = ", "), "."))
              }
            })



#' @title DEprot.RMSE show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.RMSE",
          definition =
            function(object) {
              plot = patchwork::wrap_plots(object@correlation.plots)
              print(plot)
            })


