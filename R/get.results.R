#' @title get.results
#'
#' @description Simplifies the access to the differential expression results table
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#'
#' @return A data.frame.
#'
#' @export get.results

get.results =
  function(DEprot.analyses.object,
           contrast = 1) {

    ### check object
    if (!("DEprot.analyses" %in% class(DEprot.analyses.object))) {
      stop("The input must be an object of class 'DEprot.analyses'.")
    }

    ### check and collect contrast
    if (is.numeric(contrast)) {
      if (contrast <= length(DEprot.analyses.object@analyses.result.list)) {
        data = DEprot.analyses.object@analyses.result.list[[contrast]]$results
      } else {
        stop("The 'contrast' indicated is not available.")
      }
    } else {
      stop("The 'contrast' must be a numeric value.")
    }

    ## Add contrast to data attributes
    attributes(data)$contrast = names(DEprot.analyses.object@analyses.result.list)[contrast]


    ### Export plot
    return(data)
  } # END of function

