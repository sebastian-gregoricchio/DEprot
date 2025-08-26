#' @title plot.PC.cumulative
#'
#' @description Function to plot the cumulative variance of all the principal components of a PCA.
#'
#' @param DEprot.PCA.object An object of class \code{DEprot.PCA}.
#' @param bar.color String indicating the color to use for the bar fill. Default: \code{"steelblue"}.
#' @param line.color String indicating the color to use for the line and the dots of the cumulative curve. Default: \code{"navyblue"}.
#' @param title String indicating the title of the plot (markdown annotation supported).
#'
#'
#' @return A barplot in ggplot format.
#'
#' @name plot.PC.cumulative
#'
#' @import ggplot2
#' @import ggtext
#'
#' @author Sebastian Gregoricchio
#'
#' @export plot.PC.cumulative

plot.PC.cumulative =
  function(DEprot.PCA.object,
           bar.color = "steelblue",
           line.color = "navyblue",
           title = NULL) {

    # ### Libraries
    # require(ggplot2)


    ### check object
    if (!("DEprot.PCA" %in% class(DEprot.PCA.object))) {
      stop("The input must be an object of class 'DEprot.PCA'.")
    }


    ### Generate plot
    cumulative.plot =
      ggplot(data = DEprot.PCA.object@importance,
             aes(x = PC)) +
      geom_bar(mapping = aes(y = Proportion.of.Variance),
               stat = "identity",
               fill = bar.color) +
      geom_line(mapping = aes(y = Cumulative.Proportion,
                              group = 1),
                color = line.color,
                linetype = 1) +
      geom_point(mapping = aes(y = Cumulative.Proportion,
                               group = 1),
                 stroke = NA,
                 size = 2,
                 color = line.color) +
      ylab("Proportion of variance") +
      xlab("Principal Component (PC)") +
      ggtitle(title) +
      theme_classic() +
      scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
      theme(axis.ticks.x = element_blank(),
            axis.text = element_text(color = "black"),
            axis.line.x = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
            panel.grid.minor.y = element_line(color = "gray", linewidth = 0.1),
            plot.title = ggtext::element_markdown(color = "black", hjust = 0.5),
            aspect.ratio = 0.5)

    return(cumulative.plot)
  } #END function
