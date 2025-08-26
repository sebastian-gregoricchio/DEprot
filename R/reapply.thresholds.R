#' @title reapply.thresholds
#'
#' @description Allows for the re-computation of the differential status and re-plotting volcano and MA plot for each contrast.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param linear.FC Number indicating the (absolute) fold change threshold (linear scale) to use to define differential proteins. Default: \code{2}.
#' @param linear.FC.unresp.range A numeric 2-elements vector indicating the range (linear scale) used to define the unresponsive fold changes. Default: \code{c(1/1.1, 1.1)}.
#' @param p.adjusted Numeric value indicating the p.adjusted threshold to apply to the differential analyses. Default: \code{0.05}.
#' @param up.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"indianred"}.
#' @param down.color String indicating the color to use for up-regulated proteins in the plots. Default: \code{"steelblue"}.
#' @param unresponsive.color String indicating the color to use for unresponsive proteins in the plots. Default: \code{"purple"}.
#' @param null.color String indicating the color to use for null proteins in the plots. Default: \code{"gray"}.
#'
#' @import dplyr
#' @importFrom purrr pmap
#'
#' @author Sebastian Gregoricchio
#'
#' @export reapply.thresholds



reapply.thresholds =
  function(DEprot.analyses.object,
           linear.FC = 2,
           p.adjusted = 0.05,
           linear.FC.unresp.range = c(1/1.1, 1.1),
           up.color = "indianred",
           down.color = "steelblue",
           unresponsive.color = "purple",
           null.color = "gray") {

    # ### Packages
    # require(dplyr)

    ## Define signif status function
    de.status =
      function(FC, p) {
        ifelse(p < p.adjusted,
               yes = ifelse(abs(FC) >= log2(linear.FC),
                            yes = ifelse(sign(FC) == 1,
                                         yes = DEprot.analyses.object@contrasts[[i]]$var.1,
                                         no = DEprot.analyses.object@contrasts[[i]]$var.2),
                            no = "null"),
               no = ifelse(FC >= log2(linear.FC.unresp.range[1]) & FC <= log2(linear.FC.unresp.range[2]),
                           yes = "unresponsive",
                           no = "null"))
      }


    ####################################################

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return(DEprot.analyses.object)
    }


    ## Avoid running the function if the thresholds are the same
    if (all(DEprot.analyses.object@differential.analyses.params$linear.FC.th == linear.FC,
            all(DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range == linear.FC.unresp.range),
            DEprot.analyses.object@differential.analyses.params$padj.th == p.adjusted)) {
      stop("All the fold change and p.adjusted thresholds are the same as the original ones.\nNo modifcations have been applied to the object.")
      #return(DEprot.analyses.object)
    }


    ####################################################


    ### Modify parameters
    DEprot.analyses.object@differential.analyses.params$linear.FC.th = linear.FC
    DEprot.analyses.object@differential.analyses.params$linear.FC.unresp.range = linear.FC.unresp.range
    DEprot.analyses.object@differential.analyses.params$padj.th = p.adjusted


    ### Reapply thresholds to the tables
    for (i in 1:length(DEprot.analyses.object@analyses.result.list)) {

      ## Redefine differential status
      DEprot.analyses.object@analyses.result.list[[i]]$results =
        DEprot.analyses.object@analyses.result.list[[i]]$results %>%
        dplyr::mutate(diff.status = unlist(purrr::pmap(.l = list(FC = DEprot.analyses.object@analyses.result.list[[i]]$results[,5],
                                                                 p = DEprot.analyses.object@analyses.result.list[[i]]$results$padj),
                                                       .f = function(FC,p){de.status(FC,p)}))) %>%
        dplyr::mutate(diff.status = factor(diff.status,
                                           levels = c(DEprot.analyses.object@contrasts[[i]]$var.2,
                                                      DEprot.analyses.object@contrasts[[i]]$var.1,
                                                      "unresponsive",
                                                      "null")))

      ## Recount the n of diff genes
      old_log2FC_col_name = colnames(DEprot.analyses.object@analyses.result.list[[i]]$results)[5]
      colnames(DEprot.analyses.object@analyses.result.list[[i]]$results)[5] = "log2.Fold_group1.vs.group2"

      DEprot.analyses.object@analyses.result.list[[i]]$n.diff =
        data.frame(DEprot.analyses.object@analyses.result.list[[i]]$results %>%
                     dplyr::group_by(diff.status, .drop = FALSE) %>%
                     dplyr::summarise(n = n(),
                                      median.FoldChange = median(log2.Fold_group1.vs.group2)))

      colnames(DEprot.analyses.object@analyses.result.list[[i]]$results)[5] = old_log2FC_col_name



      ## Define color params
      colors.plots = c(up.color, down.color, null.color, unresponsive.color)
      names(colors.plots) = c(DEprot.analyses.object@contrasts[[i]]$var.1, DEprot.analyses.object@contrasts[[i]]$var.2, "null", "unresponsive")


      ## replot volcano
      DEprot.analyses.object@analyses.result.list[[i]]$volcano =
        DEprot::plot.volcano(DEprot.analyses.object = DEprot.analyses.object,
                             contrast = i,
                             up.color = up.color,
                             down.color = down.color,
                             unresponsive.color = unresponsive.color,
                             null.color = null.color)


      ## replot MA
      DEprot.analyses.object@analyses.result.list[[i]]$MA.plot =
        DEprot::plot.MA(DEprot.analyses.object = DEprot.analyses.object,
                        contrast = i,
                        up.color = up.color,
                        down.color = down.color)
    }

    ## return object
    return(DEprot.analyses.object)

  } # END FUCNTION
