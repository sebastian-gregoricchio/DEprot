#' @title add.protein.info
#'
#' @description Attaches (or removes) a protein annotation table to an existing
#' \code{DEprot} or \code{DEprot.analyses} object. The table is validated and
#' re-ordered to match the proteins of the object, so that the \code{protein.info}
#' slot stays row-by-row aligned with the counts tables. The annotation can then be
#' appended to the differential expression results by \link{get.results}.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param protein.info A \code{data.frame} (or \code{matrix}) with one row per protein and any number of annotation columns (e.g., gene symbol, description, peptide counts, etc.). The protein IDs must be provided as row names, in a column called \code{prot.id}, or in the column indicated by \code{id.column}; they must correspond to the row names of the counts table. Proteins missing from the annotation are filled with \code{NA}, while annotations of proteins absent from the object are discarded. Use \code{NULL} to remove an annotation already stored in the object.
#' @param id.column String indicating the name of the column of \code{protein.info} containing the protein IDs. If \code{NULL} (default), the row names are used (or a column called \code{prot.id} when the table has no explicit row names).
#' @param overwrite Logical value indicating whether an annotation table already stored in the object should be replaced. Default: \code{FALSE}.
#'
#' @return An object of the same class as the input (\code{DEprot} or \code{DEprot.analyses}).
#'
#' @seealso \link{load.counts2}, \link{get.results}
#'
#' @importFrom methods slot
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' dpo <- load.counts2(counts = DEprot::unimputed.counts,
#'                     metadata = DEprot::sample.config,
#'                     data.type = "raw",
#'                     log.base = 2)
#'
#' # A minimal annotation table
#' info <- data.frame(gene.name = toupper(rownames(DEprot::unimputed.counts)),
#'                    row.names = rownames(DEprot::unimputed.counts))
#'
#' dpo <- add.protein.info(DEprot.object = dpo,
#'                         protein.info = info)
#'
#'
#' @export add.protein.info

add.protein.info =
  function(DEprot.object,
           protein.info,
           id.column = NULL,
           overwrite = FALSE) {

    ### check object
    if (!methods::is(DEprot.object, "DEprot")) {
      stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.", call. = FALSE)
    }

    if (missing(protein.info)) {
      stop(paste0("Provide a 'protein.info' table: a data.frame with one row per protein, ",
                  "the protein IDs as row names (or in a dedicated column) and any number of annotation columns.\n",
                  "Use `protein.info = NULL` to remove an annotation already stored in the object."),
           call. = FALSE)
    }


    ### Current annotation
    current.info = .get.protein.info(DEprot.object)


    ### Removal of the annotation
    if (is.null(protein.info)) {
      if (is.null(current.info)) {
        message("No 'protein.info' table is stored in this object: nothing to remove.")
      } else {
        message("The 'protein.info' table has been removed from the object.")
      }
      methods::slot(DEprot.object, "protein.info") = NULL
      return(DEprot.object)
    }


    ### Check overwriting
    if (!is.null(current.info) & !isTRUE(overwrite)) {
      stop(paste0("This object already contains a 'protein.info' table (columns: ",
                  paste0(colnames(current.info), collapse = ", "), ").\n",
                  "Set `overwrite = TRUE` to replace it, or combine the tables manually, e.g.:\n",
                  "  info <- cbind(<object>@protein.info, <new.columns>)"),
           call. = FALSE)
    }


    ### Collect the protein IDs of the object
    protein.ids = .get.protein.ids(DEprot.object)

    if (is.null(protein.ids)) {
      stop("No counts table with protein row names is available in this object: the 'protein.info' cannot be aligned.",
           call. = FALSE)
    }


    ### Validate and align
    methods::slot(DEprot.object, "protein.info") =
      .check.protein.info(protein.info = protein.info,
                          protein.ids = protein.ids,
                          id.column = id.column,
                          arg.name = "protein.info",
                          verbose = TRUE)

    return(DEprot.object)
  } # END function
