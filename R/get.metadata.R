#' @title get.metadata
#'
#' @description Function to extract the metadata from a \code{DEprot} object.
#'
#' @param DEprot.object Any object of class \code{DEprot}.
#'
#' @return Data.frame corresponding to the metadata of the provided object.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Extract metadata for DEprot oeejct
#' get.metadata(DEprot::dpo.imputed.counts)
#'
#' # Extract metadata for DEprot object
#' get.metadata(DEprot::test.toolbox$diff.exp.limma)
#'
#' @export get.metadata



get.metadata =
  function(DEprot.object) {

    ### check object
    if ("DEprot" %in% class(DEprot.object)) {
      return(DEprot.object@metadata)
    } else if ("DEprot.analyses" %in% class(DEprot.object)) {
      return(DEprot.object@metadata)
    } else if ("DEprot.PCA" %in% class(DEprot.object)) {
      return(DEprot.object@PCA.metadata)
    } else if ("DEprot.correlation" %in% class(DEprot.object)) {
      return(DEprot.object@corr.metadata)
    } else {
      stop("The input must be an object of class 'DEprot', 'DEprot.analyses', 'DEprot.PCA', or 'DEprot.correlation'.")
    }

  } # END function
