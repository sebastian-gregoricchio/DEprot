#' @title filter.proteins
#'
#' @description Function that allows for removing or keeping specific proteins from any \code{DEprot} object class.
#'
#' @param DEprot.object Object of class \code{DEprot}.
#' @param proteins String vector indicating a list of proteins to be use as filtering criteria. It must correspond to the counts rown names and 'prot.id' column of the analyses object.
#' @param mode String indicating the mode to use for the filtering. If 'keep' the proteins are selected and kept, while if 'remove' the proteins are discarded. Default: \code{"keep"}.
#'
#' @return An object of class \code{DEprot} (S4 vector).
#'
#' @importFrom patchwork wrap_plots
#' @importFrom purrr pmap
#' @import dplyr
#'
#' @export filter.proteins

filter.proteins =
  function(DEprot.object,
           proteins,
           mode = "keep") {

    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
      #return(DEprot.object)
    }

    ### check mode
    if (!(mode %in% c("keep", "k", "remove", "rm", "r", "rmv"))) {
      stop("The 'mode' must be one among: 'keep', 'remove'.")
      #return(DEprot.object)
    }


    # filter raw counts
    if (!is.null(DEprot.object@raw.counts)) {
      if (tolower(mode) %in% c("keep","k")) {
        DEprot.object@raw.counts = DEprot.object@raw.counts[(rownames(DEprot.object@raw.counts) %in% proteins),]
      } else {
        DEprot.object@raw.counts = DEprot.object@raw.counts[!(rownames(DEprot.object@raw.counts) %in% proteins),]
      }

      ## replot counts
      DEprot.object@boxplot.raw =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "raw",
                            violin.color = "darkorange",
                            title = DEprot.object@boxplot.raw$labels$title,
                            convert.log2 = T)
    }



    # filter normalized counts
    if (!is.null(DEprot.object@norm.counts)) {
      if (tolower(mode) %in% c("keep","k")) {
        DEprot.object@norm.counts = DEprot.object@norm.counts[(rownames(DEprot.object@norm.counts) %in% proteins),]
      } else {
        DEprot.object@norm.counts = DEprot.object@norm.counts[!(rownames(DEprot.object@norm.counts) %in% proteins),]
      }

      ## replot counts
      DEprot.object@boxplot.norm =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "normalized",
                            violin.color = "purple",
                            title = DEprot.object@boxplot.norm$labels$title,
                            subtitle = DEprot.object@boxplot.norm$labels$subtitle,
                            convert.log2 = TRUE)
    }


    # filter imputed counts
    if (!is.null(DEprot.object@imputed.counts)) {
      if (tolower(mode) %in% c("keep","k")) {
        DEprot.object@imputed.counts = DEprot.object@imputed.counts[(rownames(DEprot.object@imputed.counts) %in% proteins),]
      } else {
        DEprot.object@imputed.counts = DEprot.object@imputed.counts[!(rownames(DEprot.object@imputed.counts) %in% proteins),]
      }

      ## replot counts
      DEprot.object@boxplot.imputed =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "imputed",
                            violin.color = "forestgreen",
                            title = DEprot.object@boxplot.imputed$labels$title,
                            subtitle = DEprot.object@boxplot.imputed$labels$subtitle,
                            convert.log2 = TRUE)
    }



    ### FILTER ANALYSES ##################
    if ("DEprot.analyses" %in% class(DEprot.object)) {
      ## Filter the protein names
      DEprot.object@analyses.result.list =
        lapply(DEprot.object@analyses.result.list,
               function(x) {
                 ### Filtering analyses results
                 if (tolower(mode) %in% c("keep","k")) {
                   x$results = dplyr::filter(.data = x$results, prot.id %in% proteins)
                 } else {
                   x$results = dplyr::filter(.data = x$results, !(prot.id %in% proteins))
                 }


                 ### Recompute ndiff
                 diff.tb = x$results
                 colnames(diff.tb)[5] = "log2.Fold_group1.vs.group2"

                 x$n.diff =
                   data.frame(diff.tb %>%
                                dplyr::group_by(diff.status, .drop = FALSE) %>%
                                dplyr::summarise(n = n(),
                                                 median.FoldChange = median(log2.Fold_group1.vs.group2)))
                 return(x)
               })


      ## Recompute PCA and corr
      DEprot.object@analyses.result.list =
        purrr::pmap(.l = list(res = DEprot.object@analyses.result.list,
                              contr = DEprot.object@contrasts),
                    .f = function(res, contr) {
                      samples = unique(c(contr$group.1, contr$group.2))

                      res$PCA.data = DEprot::perform.PCA(DEprot.object = DEprot.object,
                                                         sample.subset = samples,
                                                         which.data = "imputed")

                      pca.123 = DEprot::plot.PC.scatter.123(DEprot.PCA.object = res$PCA.data, color.column = contr$metadata.column)
                      pca.cumulative = DEprot::plot.PC.cumulative(DEprot.PCA.object = res$PCA.data)
                      res$PCA.plots = patchwork::wrap_plots(pca.123, pca.cumulative, ncol = 1)

                      pearson = DEprot::plot.correlation.heatmap(DEprot.object = DEprot.object,
                                                                 correlation.method = "pearson",
                                                                 sample.subset = samples,
                                                                 which.data = "imputed")

                      spearman = DEprot::plot.correlation.heatmap(DEprot.object = DEprot.object,
                                                                  correlation.method = "spearman",
                                                                  sample.subset = samples,
                                                                  which.data = "imputed")

                      res$correlations = patchwork::wrap_plots(pearson@heatmap, spearman@heatmap, nrow = 1)

                      return(res)
                    })


      ### replot volcanos
      for (i in 1:length(DEprot.object@analyses.result.list)) {
        DEprot.object@analyses.result.list[[i]]$volcano = DEprot::plot.volcano(DEprot.object, contrast = i)
        DEprot.object@analyses.result.list[[i]]$MA.plot = DEprot::plot.MA(DEprot.object, contrast = i)
      }
    } # end analyses


    ### Export the filtered object
    return(DEprot.object)
  } #END function
