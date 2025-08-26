#' @title contrast.scatter
#'
#' @description Plots a scatter plot of the log2(fold change expression) derived from two differential analyses contrasts.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast.x Integer indicating the position of the contrast to use for the x-axis of the plot.
#' @param contrast.y Integer indicating the position of the contrast to use for the y-axis of the plot.
#' @param regression.line.color String indicating any R-supported color to use for the regression line and error. Default: \code{"firebrick"}.
#' @param correlation.method String indicating the clustering method to use to generate the correlation matrix. Possible options: 'pearson', 'spearman', 'kendall'. Default: \code{"pearson"}.
#' @param add.foldchange.threshold Logical value to indicate whether two gray rectangles should be used to highlight the non-differential area in the plot, based on the foldchange threshold indicated during the differential analyses. Default: \code{TRUE}.
#' @param symmetric.axes Logical value indicating whether the axes should be forces to be symmetric between x and y. Default: \code{TRUE}.
#'
#' @return A scatter plot of class ggplot2.
#'
#' @import dplyr
#' @import ggplot2
#' @import ggtext
#' @importFrom stringr str_to_title
#' @importFrom ggpubr stat_cor
#'
#' @export contrast.scatter




contrast.scatter =
  function(DEprot.analyses.object,
           contrast.x,
           contrast.y,
           regression.line.color = "firebrick",
           correlation.method = "pearson",
           add.foldchange.threshold = TRUE,
           symmetric.axes = TRUE) {

    # ### Libraries
    # require(dplyr)
    # require(ggplot2)


    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
      #return(invisible())
    }

    ### check and collect contrast
    if (!(contrast.x %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
      stop("The 'contrast.x' is not available in the DEprot.analyses object.")
      #return(invisible())
    } else if (!(contrast.y %in% 1:length(DEprot.analyses.object@analyses.result.list))) {
      stop("The 'contrast.y' is not available in the DEprot.analyses object.")
      #return(invisible())
    }


    ### get FC threshold
    fc_th = log2(DEprot.analyses.object@differential.analyses.params$linear.FC.th)



    ### collect the fc table
    fc.x = DEprot.analyses.object@analyses.result.list[[contrast.x]]$results
    fc.x = fc.x[,c(1,5)]
    colnames(fc.x)[2] = "contrast.x"

    fc.y = DEprot.analyses.object@analyses.result.list[[contrast.y]]$results
    fc.y = fc.y[,c(1,5)]
    colnames(fc.y)[2] = "contrast.y"

    fc.tb = dplyr::full_join(x = fc.x, y = fc.y, by = "prot.id")



    ### generate scatter
    scatter =
      ggplot(data = fc.tb,
             aes(x = contrast.x,
                 y = contrast.y)) +
      geom_hline(yintercept = 0, colour = "gray", linetype = 2) +
      geom_vline(xintercept = 0, colour = "gray", linetype = 2) +
      geom_point(stroke = NA,
                 size = 3,
                 alpha = 0.25,
                 show.legend = FALSE) +
      geom_rug(alpha = 0.5)

    if (add.foldchange.threshold == TRUE) {
      scatter =
        scatter +
        geom_rect(data = data.frame(xmin = c(-Inf, -fc_th, -fc_th),
                                    xmax = c(+Inf, fc_th, fc_th),
                                    ymin = c(-fc_th, fc_th, -fc_th),
                                    ymax = c(fc_th, +Inf, -Inf)),
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "gray",
                  color = NA,
                  alpha = 0.25)
    }

    scatter =
      scatter +
      geom_smooth(method = glm, formula = y ~ x, color = regression.line.color, fill = regression.line.color) +
      ggtitle(paste0("**",stringr::str_to_title(correlation.method)," correlation**")) +
      ggpubr::stat_cor(method = correlation.method,
                       cor.coef.name = dplyr::case_when(tolower(correlation.method) == "pearson" ~ "R",
                                                        tolower(correlation.method) == "spearman" ~ "rho",
                                                        tolower(correlation.method) == "kendall" ~ "tau")) +
      xlab(paste0("log~2~(Fold Change<sub>",DEprot.analyses.object@contrasts[[contrast.x]]$var.1,"</sup>&frasl;<sub>",DEprot.analyses.object@contrasts[[contrast.x]]$var.2,"</sub></sub>)")) +
      ylab(paste0("log~2~(Fold Change<sub>",DEprot.analyses.object@contrasts[[contrast.y]]$var.1,"</sup>&frasl;<sub>",DEprot.analyses.object@contrasts[[contrast.y]]$var.2,"</sub></sub>)")) +
      theme_classic() +
      theme(axis.title = ggtext::element_markdown(color = "black"),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            aspect.ratio = 1)



    if (symmetric.axes == TRUE) {
      build = ggplot_build(scatter)
      max = max(c(abs(build$layout$panel_params[[1]]$x.range), abs(build$layout$panel_params[[1]]$y.range)))
      scatter = scatter + scale_x_continuous(limits = c(-1,1)*max) + scale_y_continuous(limits = c(-1,1)*max)
    }



    ### export plot
    return(scatter)

  } # END function
