#' @title DEprot class
#' @name DEprot
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
#' @exportClass DEprot.upset

setClass(Class = "DEprot.upset",
         slots = list(upset = "ANY",
                      obs.matrix = "ANY"))



#' @title DEprot.contrast.heatmap class
#' @name DEprot.contrast.heatmap
#' @exportClass DEprot.contrast.heatmap

setClass(Class = "DEprot.contrast.heatmap",
         slots = list(heatmap = "ANY",
                      cluster = "ANY"))


#' @title DEprot.counts.heatmap class
#' @name DEprot.counts.heatmap
#' @exportClass DEprot.counts.heatmap

setClass(Class = "DEprot.counts.heatmap",
         slots = list(heatmap = "ANY",
                      row.cluster = "ANY",
                      column.cluster = "ANY"))



#' @title DEprot.enrichResult class
#' @name DEprot.enrichResult
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
#' @exportClass DEprot.pvalues

setClass(Class = "DEprot.pvalues",
         slots = list(pvalue.distribution = "ANY",
                      padjusted.distribution = "ANY",
                      pvalue.rank = "ANY"))




#' @title DEprot.normality class
#' @name DEprot.normality
#' @exportClass DEprot.normality

setClass(Class = "DEprot.normality",
         slots = list(norm.statement = "ANY",
                      norm.AD.tests = "ANY",
                      qqplots = "ANY",
                      densities = "ANY"))





#' @title DEprot.RMSE class
#' @name DEprot.RMSE
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
                      correlation.plots = "ANY"))


################# METHODS ################# "DEprot.analyses"

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
              require(patchwork, quietly = T)
              plot = (object@pvalue.distribution / object@padjusted.distribution) | object@pvalue.rank
              print(plot)
            })


#' @title DEprot.normality show-method
#' @export
setMethod(f = "show",
          signature = "DEprot.normality",
          definition =
            function(object) {
              if (object@norm.statement == T) {
                message("All samples display a normal distribution.")

              } else {
                normality = sapply(X = object@norm.AD.tests, FUN = function(x){x$p.value < p.threshold}, USE.NAMES = T)
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

