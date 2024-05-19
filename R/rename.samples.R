#' @title rename.samples
#'
#' @description Function that allows for the renaming of the columns of the counts table in a \code{DEprot} object.
#'
#' @param DEprot.object An object of class \code{DEprot}.
#' @param metadata.column A string indicating any column from the \code{metadata} table to use to rename the \code{counts} table columns. Default \code{"column.id"} (no renaming).
#'
#' @return An object of class \code{DEprot}. A column called \code{old.column.id} will be added to the \code{metadata} in order to keep track of the original names.
#'
#' @export rename.samples

rename.samples =
  function(DEprot.object,
           metadata.column = "column.id") {

    ### check object
    if (!("DEprot" %in% class(DEprot.object))) {
      warning("The input must be an object of class 'DEprot'.")
      ibmAcousticR:::stop_quietly("")
    }

    ### Check column provided
    if (!(metadata.column %in% colnames(DEprot.object$metadata))) {
      warning(paste0("The column '", metadata.column, "' is not among the column names of the metadata table.\n",
                     "Please chose one among: ", paste0(colnames(DEprot.object$metadata), collapse = ", ")))
      ibmAcousticR:::stop_quietly("")
    } else if (length(DEprot.object$metadata[,metadata.column]) != length(unique(DEprot.object$metadata[,metadata.column]))) {
      warning(paste0("The column '", metadata.column, "' does not contain unique identifiers."))
      ibmAcousticR:::stop_quietly("")
    }


    ### Rename columns
    # Add a column with the old names to metadata
    DEprot.object$metadata$old.column.id = DEprot.object$metadata$column.id
    DEprot.object$metadata$column.id = DEprot.object$metadata[,metadata.column]

    # rename raw counts
    if (!is.null(DEprot.object$raw.counts)) {
      for (i in 1:nrow(DEprot.object$metadata)) {
        colnames(DEprot.object$raw.counts)[which(colnames(DEprot.object$raw.counts) == DEprot.object$metadata$old.column.id[i])] = DEprot.object$metadata$column.id[i]
      }
    }

    # rename normalized counts
    if (!is.null(DEprot.object$norm.counts)) {
      for (i in 1:nrow(DEprot.object$metadata)) {
        colnames(DEprot.object$norm.counts)[which(colnames(DEprot.object$norm.counts) == DEprot.object$metadata$old.column.id[i])] = DEprot.object$metadata$column.id[i]
      }
    }

    return(DEprot.object)
  } #END function
