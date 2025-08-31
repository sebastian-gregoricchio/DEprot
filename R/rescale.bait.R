#' @title rescale.bait
#'
#' @description Rescaling of the imputed counts based on the mean of a specific protein. This is useful for comparing the enrichment of interactors in the RIME experiments.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses} (must contain imputed counts).
#' @param bait.id String indicating the id of a protein of the dataset present in the DEprot object.
#'
#' @return An object of class \code{DEprot} or \code{DEprot.analyses} (same class of the input).
#'
#' @import dplyr
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' recaled_dpo <- rescale.bait(DEprot.object = DEprot::test.toolbox$dpo.imp, bait.id = "protein.29")
#'
#' @export rescale.bait


rescale.bait =
  function(DEprot.object,
           bait.id){

    # ### Libraries
    # require(dplyr)


    ### check object and extract metadata table
    if (!("DEprot" %in% class(DEprot.object))) {
      if (!("DEprot.analyses" %in% class(DEprot.object))) {
        stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.")
        #return(DEprot.object)
      }
    }


    ### Extract the imputed data
    if (isFALSE(DEprot.object@imputed)) {
      stop("The DEprot.object provided does not contain imputed data. The bai rescaling can be applied only to imputed data.")
      #return(DEprot.object)
    } else {
      mat = DEprot.object@imputed.counts
    }


    ### Check bait
    if (!(bait.id %in% rownames(mat))) {
      stop(paste0("The bait.id '", bait.id, "' is not available in the dataset."))
      #return(DEprot.object)
    }


    ############## PERFORM RESCALING ##############

    ### Convert table in linear values
    mat.linear = DEprot.object@log.base^mat

    ### Compute scaling factor
    bait.scores = mat.linear[bait.id,]
    bait.mean = mean(bait.scores)
    scaling.factors = bait.scores / bait.mean

    ### applying scaling factor
    mat.scaled.linear = sweep(x = mat.linear,
                              MARGIN = 2,
                              STATS = scaling.factors,
                              FUN = "/")

    ### Convert to 1 (0 in log) values below 1, and re-transform in log
    mat.scaled.linear[mat.scaled.linear < 1] = 1
    mat.scaled.log = log(x = mat.scaled.linear, base = DEprot.object@log.base)


    ### make scaling method
    scaling.method =
      list(original.imputation = DEprot.object@imputation,
           bait.rescaling = list(source.matrix = mat,
                                 bait.id = bait.id,
                                 bait.scores = data.frame(sample.id = names(bait.scores),
                                                          logarithmic.intensity = log(bait.scores, DEprot.object@log.base),
                                                          linear.intensity = bait.scores,
                                                          scaling.factor = scaling.factors,
                                                          rescaled.linear = log(bait.mean, DEprot.object@log.base),
                                                          rescaled.linear = bait.mean)))


    ### Return new rescaled object
    rescaled.dpo = DEprot.object
    rescaled.dpo@imputed.counts = mat.scaled.log
    rescaled.dpo@imputation = scaling.method

    return(rescaled.dpo)

  } # END function





