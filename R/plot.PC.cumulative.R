#' @title plot.PC.cumulative
#'
#' @description Function to plot the cumulative variance of all the principal components of a PCA.
#'
#' @param DEprot.PCA.object An object of class \code{DEprot.PCA}.
#' @param bar.color String indicating the color to use for the bar fill. Default: \code{"steelblue"}.
#' @param line.color String indicating the color to use for the line and the dots of the cumulative curve. Default: \code{"navyblue"}.
#' @param title String indicating the title of the plot (markdown annotation supported).
#'
#' @return A barplot in ggplot format.
#'
#' @seealso \link{perform.PCA}, \link{plot.PCoA.cumulative}
#'
#' @name plot.PC.cumulative
#'
#' @import ggplot2
#' @import ggtext
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Perform Principal Component Analyses (PCA)
#' pca <- perform.PCA(DEprot.object = DEprot::test.toolbox$dpo.imp)
#'
#' # Plot cumulative variance barplot
#' plot.PC.cumulative(pca)
#'
#' @export plot.PC.cumulative

plot.PC.cumulative =
  function(DEprot.PCA.object,
           bar.color = "steelblue",
           line.color = "navyblue",
           title = NULL) {


    ### check object
    if (!("DEprot.PCA" %in% class(DEprot.PCA.object))) {
      stop("The input must be an object of class 'DEprot.PCA'.")
    }

    ### harmonized access to the ordination slots
    ord = .ordination.slots(DEprot.PCA.object)


    ### Generate plot (shared body with plot.PCoA.cumulative)
    cumulative.plot =
      .cumulative.plot(importance = ord$importance,
                       axis.prefix = ord$axis.prefix,
                       axis.title = ord$axis.title,
                       bar.color = bar.color,
                       line.color = line.color,
                       title = title,
                       broken.stick = FALSE,
                       filter.positive = FALSE)

    return(cumulative.plot)
  } #END function
