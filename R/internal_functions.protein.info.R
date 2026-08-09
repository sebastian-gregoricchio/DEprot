#' @title .get.protein.info
#'
#' @description Safe accessor for the \code{protein.info} slot. Returns \code{NULL}
#' both when the slot is empty and when the object was created with a version of
#' \code{DEprot} predating the introduction of the slot (backward compatibility).
#'
#' @param object An object of class \code{DEprot} or \code{DEprot.analyses}.
#'
#' @return A \code{data.frame} or \code{NULL}.
#'
#' @importFrom methods .hasSlot slot
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.get.protein.info =
  function(object) {

    if (!methods::.hasSlot(object, "protein.info")) {return(NULL)}

    info = methods::slot(object, "protein.info")

    if (is.null(info)) {return(NULL)}
    if (length(info) == 0) {return(NULL)}
    if (is.atomic(info) && length(info) == 1 && is.na(info)) {return(NULL)}

    return(as.data.frame(info))
  } # END function




#' @title .get.protein.ids
#'
#' @description Extracts the protein IDs (row names) from the first available
#' counts table of a \code{DEprot}/\code{DEprot.analyses} object. All the counts
#' tables of an object share the same rows, hence the first available one is
#' representative of the whole object.
#'
#' @param object An object of class \code{DEprot} or \code{DEprot.analyses}.
#'
#' @return A character vector, or \code{NULL} if no counts table is available.
#'
#' @importFrom methods slot
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.get.protein.ids =
  function(object) {

    count.slots = c("raw.counts", "norm.counts", "random.counts", "imputed.counts")

    for (s in count.slots) {
      mat = methods::slot(object, s)
      if (!.deprot_slot_is_empty(mat)) {
        if (!is.null(rownames(mat))) {
          return(rownames(mat))
        }
      }
    }

    return(NULL)
  } # END function




#' @title .check.protein.info
#'
#' @description Validates a protein annotation table and aligns it to a set of
#' protein IDs. The table is re-ordered to strictly follow the order of
#' \code{protein.ids}: proteins absent from the annotation are filled with
#' \code{NA}, while annotations of proteins absent from the counts are dropped.
#' This guarantees that the \code{protein.info} slot stays row-by-row aligned
#' with the counts tables at any moment.
#'
#' @param protein.info A \code{data.frame} (or \code{matrix}) with one row per protein, or \code{NULL}.
#' @param protein.ids String vector of the protein IDs to which the table must be aligned (usually the row names of the counts).
#' @param id.column String indicating the name of the column of \code{protein.info} containing the protein IDs. If \code{NULL} (default) the row names are used; when the table has no explicit row names a column named \code{prot.id} is looked up.
#' @param arg.name String used in the error/warning messages to refer to the table. Default: \code{"protein.info"}.
#' @param verbose Logical value indicating whether to print messages about the alignment. Default: \code{TRUE}.
#'
#' @return A \code{data.frame} with \code{length(protein.ids)} rows and \code{protein.ids} as row names, or \code{NULL}.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.check.protein.info =
  function(protein.info,
           protein.ids,
           id.column = NULL,
           arg.name = "protein.info",
           verbose = TRUE) {

    ### Nothing to do
    if (is.null(protein.info)) {return(NULL)}
    if (length(protein.info) == 0) {return(NULL)}
    if (is.atomic(protein.info) && length(protein.info) == 1 && is.na(protein.info)) {return(NULL)}


    ### Coerce to data.frame
    if (is.matrix(protein.info)) {
      protein.info = as.data.frame(protein.info, stringsAsFactors = FALSE)
    } else if (is.data.frame(protein.info)) {
      protein.info = as.data.frame(protein.info, stringsAsFactors = FALSE)
    } else {
      stop(paste0("The '", arg.name, "' must be a data.frame (or a matrix) with one row per protein.\n",
                  "The protein IDs must be provided as row names or in a dedicated column (see 'id.column')."),
           call. = FALSE)
    }

    if (ncol(protein.info) == 0) {
      stop(paste0("The '", arg.name, "' table does not contain any column."), call. = FALSE)
    }


    ### Recover the protein IDs: dedicated column > row names > 'prot.id' column
    if (!is.null(id.column)) {
      if (!(id.column %in% colnames(protein.info))) {
        stop(paste0("The column '", id.column, "' is not present in the '", arg.name, "' table.\n",
                    "Available columns: ", paste0(colnames(protein.info), collapse = ", "), "."),
             call. = FALSE)
      }
      info.ids = as.character(protein.info[, id.column])
      protein.info = protein.info[, setdiff(colnames(protein.info), id.column), drop = FALSE]

    } else if (!is.null(rownames(protein.info)) & .row_names_info(protein.info) > 0) {
      info.ids = as.character(rownames(protein.info))

    } else if ("prot.id" %in% colnames(protein.info)) {
      info.ids = as.character(protein.info[, "prot.id"])
      protein.info = protein.info[, setdiff(colnames(protein.info), "prot.id"), drop = FALSE]

    } else {
      stop(paste0("The '", arg.name, "' table has no protein IDs: provide them as row names, ",
                  "in a column called 'prot.id', or indicate the column to use through 'id.column'."),
           call. = FALSE)
    }

    if (ncol(protein.info) == 0) {
      stop(paste0("The '", arg.name, "' table contains only the protein IDs and no annotation column."), call. = FALSE)
    }


    ### Check the IDs
    if (any(is.na(info.ids)) | any(info.ids == "")) {
      stop(paste0("The protein IDs of the '", arg.name, "' table contain missing values ('') or NAs. ",
                  "Replace them with actual values."), call. = FALSE)
    }

    if (any(duplicated(info.ids))) {
      dup.ids = unique(info.ids[duplicated(info.ids)])
      stop(paste0("The '", arg.name, "' table contains duplicated protein IDs (n=", length(dup.ids), "): ",
                  paste0(dup.ids[1:min(10, length(dup.ids))], collapse = ", "),
                  ifelse(length(dup.ids) > 10, ", ...", ""), ".\n",
                  "Only one row per protein is allowed."),
           call. = FALSE)
    }


    ### Align to the protein IDs of the object
    protein.ids = as.character(protein.ids)

    missing.ids = setdiff(protein.ids, info.ids)
    extra.ids = setdiff(info.ids, protein.ids)

    aligned.info = protein.info[match(protein.ids, info.ids), , drop = FALSE]
    rownames(aligned.info) = protein.ids


    ### Report
    if (length(missing.ids) == length(protein.ids)) {
      warning(paste0("None of the protein IDs of the '", arg.name, "' table matches the proteins of the object: ",
                     "the annotation columns will contain only NAs.\n",
                     "Check that the IDs used in the annotation correspond to the row names of the counts table."),
              call. = FALSE, immediate. = TRUE)

    } else if (isTRUE(verbose)) {
      if (length(missing.ids) > 0) {
        message(paste0("[", arg.name, "] ", length(missing.ids), " protein(s) of the object ",
                       ifelse(length(missing.ids) == 1, "is", "are"),
                       " not annotated and will be filled with NAs: ",
                       paste0(missing.ids[1:min(5, length(missing.ids))], collapse = ", "),
                       ifelse(length(missing.ids) > 5, ", ...", ""), "."))
      }

      if (length(extra.ids) > 0) {
        message(paste0("[", arg.name, "] ", length(extra.ids), " annotated protein(s) ",
                       ifelse(length(extra.ids) == 1, "is", "are"),
                       " not present in the object and ",
                       ifelse(length(extra.ids) == 1, "has", "have"), " been discarded."))
      }
    }

    return(aligned.info)
  } # END function
