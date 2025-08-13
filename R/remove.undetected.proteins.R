#' @title remove.undetected.proteins
#'
#' @description This function allows for the removal of all the proteins that occur below a certain threshold (detected in at least n samples).
#'
#' @param DEprot.object A \code{DEprot} object.
#' @param min.n.samples Numeric value indicating the minimum number of sample in which a protein must be detected to be kept. Default: \code{3}.
#' @param which.data String indicating which type of counts should be used to define the proteins to remove. One among: 'raw', 'normalized', 'norm', 'imputed', 'imp'. Default: \code{"normalized"}.
#'
#' @return A \code{DEprot} object.
#'
#' @export remove.undetected.proteins


remove.undetected.proteins =
  function(DEprot.object,
           min.n.samples = 3,
           which.data = "normalized") {

    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
      #return(DEprot.object)
    }


    ### Check and extract table
    if (tolower(which.data) == "raw") {
      if (!is.null(DEprot.object@raw.counts)) {
        mat = DEprot.object@raw.counts
        data.used = "raw"
      } else {
        stop(paste0("Use of RAW counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("norm", "normalized", "normal")) {
      if (!is.null(DEprot.object@norm.counts)) {
        mat = DEprot.object@norm.counts
        data.used = "normalized"
      } else {
        stop(paste0("Use of NORMALIZED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else if (tolower(which.data) %in% c("imputed", "imp", "impute")) {
      if (!is.null(DEprot.object@imputed.counts)) {
        mat = DEprot.object@imputed.counts
        data.used = "imputed"
      } else {
        stop(paste0("Use of IMPUTED counts was required, but not available.\n",
                    "       Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
        #return(DEprot.object)
      }
    } else {
      stop(paste0("The 'which.data' value is not recognized.\n",
                  "Please indicated a count type among 'raw', 'normalized', 'imputed', using the option 'which.data'."))
      #return(DEprot.object)
    }



    ### Define the proteins to remove
    n.present.values = rowSums(!is.na(mat))

    if (TRUE %in% (n.present.values < min.n.samples)) {
      prot.ids.to.remove = names(n.present.values[n.present.values < min.n.samples])
      message(paste0("These proteins (n=", length(unique(prot.ids.to.remove)),") will be removed:\n", paste0(unique(prot.ids.to.remove), collapse = ", ")))
      return(DEprot::filter.proteins(DEprot.object = DEprot.object, proteins = unique(prot.ids.to.remove), mode = "remove"))
    } else {
      return(DEprot.object)
    }

  } # END function
