#' @title get.protein.info
#'
#' @description Function to extract the protein annotation table from a \code{DEprot} object.
#'
#' @param DEprot.object Any object of class \code{DEprot} or \code{DEprot.analyses}.
#'
#' @return Data.frame corresponding to the \code{protein.info} slot of the provided object (row names: protein IDs), or \code{NULL} when no annotation is stored.
#'
#' @seealso \link{add.protein.info}, \link{load.counts2}, \link{get.results}
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Extract the protein annotation from a DEprot object
#' get.protein.info(DEprot::dpo.imputed.counts)
#'
#' # Extract the protein annotation from a DEprot.analyses object
#' get.protein.info(DEprot::test.toolbox$diff.exp.limma)
#'
#' @export get.protein.info



get.protein.info =
  function(DEprot.object) {

    ### check object
    if (!methods::is(DEprot.object, "DEprot")) {
      stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
    }

    protein.info = .get.protein.info(DEprot.object)

    if (is.null(protein.info)) {
      message("No 'protein.info' table is stored in this object. An annotation can be added with `add.protein.info()`.")
      return(invisible(NULL))
    }

    return(protein.info)

  } # END function
