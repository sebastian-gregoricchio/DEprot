#' @title rename.samples
#'
#' @description Function that allows for the renaming of the columns of the counts table in a \code{DEprot} object.
#'
#' @param DEprot.object An object of class \code{DEprot}.
#' @param metadata.column A string indicating any column from the \code{metadata} table to use to rename the \code{counts} table columns. Default \code{"column.id"} (no renaming).
#'
#' @return An object of class \code{DEprot} (S4 vector). A column called \code{old.column.id} will be added to the \code{metadata} in order to keep track of the original names.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' dpo_renamed <- rename.samples(DEprot.object = DEprot::test.toolbox$dpo.imp, metadata.column = "sample.id")
#'
#' @export rename.samples

rename.samples =
  function(DEprot.object,
           metadata.column = "column.id") {

    ### check object
    if (!("DEprot" %in% class(DEprot.object)) & !("DEprot.analyses" %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot'.")
      #return(DEprot.object)
    }

    ### Check column provided
    if (!(metadata.column %in% colnames(DEprot.object@metadata))) {
      stop(paste0("The column '", metadata.column, "' is not among the column names of the metadata table.\n",
                  "       Please chose one among: ", paste0(colnames(DEprot.object@metadata), collapse = ", ")))
      #return(DEprot.object)
    } else if (length(DEprot.object@metadata[,metadata.column]) != length(unique(DEprot.object@metadata[,metadata.column]))) {
      stop(paste0("The column '", metadata.column, "' does not contain unique identifiers."))
      #return(DEprot.object)
    }


    ### Rename columns
    # Add a column with the old names to metadata
    DEprot.object@metadata$old.column.id = DEprot.object@metadata$column.id
    DEprot.object@metadata$column.id = DEprot.object@metadata[,metadata.column]

    # rename raw counts
    if (!is.null(DEprot.object@raw.counts)) {
      for (i in 1:nrow(DEprot.object@metadata)) {
        colnames(DEprot.object@raw.counts)[which(colnames(DEprot.object@raw.counts) == DEprot.object@metadata$old.column.id[i])] = DEprot.object@metadata$column.id[i]
      }
      ## replot counts
      DEprot.object@boxplot.raw =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "raw",
                            violin.color = "darkorange",
                            title = DEprot.object@boxplot.raw$labels$title,
                            convert.log2 = TRUE)

    }



    # rename normalized counts
    if (!is.null(DEprot.object@norm.counts)) {
      for (i in 1:nrow(DEprot.object@metadata)) {
        colnames(DEprot.object@norm.counts)[which(colnames(DEprot.object@norm.counts) == DEprot.object@metadata$old.column.id[i])] = DEprot.object@metadata$column.id[i]
      }

      ## replot counts
      DEprot.object@boxplot.norm =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "normalized",
                            violin.color = "purple",
                            title = DEprot.object@boxplot.norm$labels$title,
                            subtitle = DEprot.object@boxplot.norm$labels$subtitle,
                            convert.log2 = TRUE)
    }


    # rename imputed counts
    if (!is.null(DEprot.object@imputed.counts)) {
      for (i in 1:nrow(DEprot.object@metadata)) {
        colnames(DEprot.object@imputed.counts)[which(colnames(DEprot.object@imputed.counts) == DEprot.object@metadata$old.column.id[i])] = DEprot.object@metadata$column.id[i]
      }

      ## replot counts
      DEprot.object@boxplot.imputed =
        DEprot::plot.counts(DEprot.object = DEprot.object,
                            which.data = "imputed",
                            violin.color = "forestgreen",
                            title = DEprot.object@boxplot.imputed$labels$title,
                            subtitle = DEprot.object@boxplot.imputed$labels$subtitle,
                            convert.log2 = TRUE)
    }

    return(DEprot.object)
  } #END function
