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
      stop("The '", arg.name, "' must be a data.frame (or a matrix) with one row per protein.\n",
           "The protein IDs must be provided as row names or in a dedicated column (see 'id.column').",
           call. = FALSE)
    }

    if (ncol(protein.info) == 0) {
      stop("The '", arg.name, "' table does not contain any column.", call. = FALSE)
    }


    ### Recover the protein IDs: dedicated column > row names > 'prot.id' column
    if (!is.null(id.column)) {
      if (!(id.column %in% colnames(protein.info))) {
        stop("The column '", id.column, "' is not present in the '", arg.name, "' table.\n",
             "Available columns: ", paste0(colnames(protein.info), collapse = ", "), ".",
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
      stop("The '", arg.name, "' table has no protein IDs: provide them as row names, ",
           "in a column called 'prot.id', or indicate the column to use through 'id.column'.",
           call. = FALSE)
    }

    if (ncol(protein.info) == 0) {
      stop("The '", arg.name, "' table contains only the protein IDs and no annotation column.", call. = FALSE)
    }


    ### Check the IDs
    if (any(is.na(info.ids)) | any(info.ids == "")) {
      stop("The protein IDs of the '", arg.name, "' table contain missing values ('') or NAs. ",
           "Replace them with actual values.", call. = FALSE)
    }

    if (any(duplicated(info.ids))) {
      dup.ids = unique(info.ids[duplicated(info.ids)])
      stop("The '", arg.name, "' table contains duplicated protein IDs (n=", length(dup.ids), "): ",
           paste0(dup.ids[1:min(10, length(dup.ids))], collapse = ", "),
           ifelse(length(dup.ids) > 10, ", ...", ""), ".\n",
           "Only one row per protein is allowed.",
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
      warning("None of the protein IDs of the '", arg.name, "' table matches the proteins of the object: ",
              "the annotation columns will contain only NAs.\n",
              "Check that the IDs used in the annotation correspond to the row names of the counts table.",
              call. = FALSE, immediate. = TRUE)

    } else if (isTRUE(verbose)) {
      if (length(missing.ids) > 0) {
        message("[", arg.name, "] ", length(missing.ids), " protein(s) of the object ",
                ifelse(length(missing.ids) == 1, "is", "are"),
                " not annotated and will be filled with NAs: ",
                paste0(missing.ids[1:min(5, length(missing.ids))], collapse = ", "),
                ifelse(length(missing.ids) > 5, ", ...", ""), ".")
      }

      if (length(extra.ids) > 0) {
        message("[", arg.name, "] ", length(extra.ids), " annotated protein(s) ",
                ifelse(length(extra.ids) == 1, "is", "are"),
                " not present in the object and ",
                ifelse(length(extra.ids) == 1, "has", "have"), " been discarded.")
      }
    }

    return(aligned.info)
  } # END function




# ----------------------------------------------------------------------------------------


#' @title .append.protein.info
#'
#' @description Internal. Appends the columns of the \code{protein.info} table to a results table,
#' matching the rows by protein ID rather than by position, since the results might have been
#' re-ordered or filtered.
#'
#' @param data data.frame containing a \code{prot.id} column.
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}, or \code{NULL}.
#' @param protein.info.columns String vector of the columns to append, or one of the keywords \code{"none"} and \code{"all"}.
#' @param protein.info.prefix String to prepend to the names of the appended columns, or \code{NULL}.
#'
#' @return The input data.frame, with the annotation columns appended when available.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.append.protein.info =
  function(data,
           DEprot.object = NULL,
           protein.info.columns = "none",
           protein.info.prefix = NULL) {

    if (is.null(protein.info.columns)) {protein.info.columns = "none"}

    if (!is.character(protein.info.columns) | length(protein.info.columns) == 0) {
      stop("The 'protein.info.columns' must be a string vector indicating the annotation columns to append, ",
           "or one of the keywords 'none' (default) and 'all'.",
           call. = FALSE)
    }

    keyword = "custom"

    if (length(protein.info.columns) == 1) {
      if (tolower(protein.info.columns[1]) %in% c("none", "all")) {
        keyword = tolower(protein.info.columns[1])
      }
    }

    if (keyword == "none") {return(data)}

    if (is.null(DEprot.object)) {
      warning("No 'DEprot.object' was provided: the protein annotation cannot be appended.",
              call. = FALSE, immediate. = TRUE)
      return(data)
    }

    protein.info = .get.protein.info(DEprot.object)

    if (is.null(protein.info)) {
      warning("No 'protein.info' table is stored in this object: the requested columns cannot be appended.\n",
              "An annotation can be added with `add.protein.info()`.",
              call. = FALSE, immediate. = TRUE)
      return(data)
    }


    ### Subset the columns, if required
    if (keyword == "custom") {
      missing.columns = setdiff(protein.info.columns, colnames(protein.info))

      if (length(missing.columns) > 0) {
        stop("The following 'protein.info.columns' are not available in the annotation table: ",
             paste(missing.columns, collapse = ", "), ".\n",
             "Available columns: ", paste(colnames(protein.info), collapse = ", "), ".",
             call. = FALSE)
      }

      protein.info = protein.info[, protein.info.columns, drop = FALSE]
    }


    ### Add the prefix, if required
    if (!is.null(protein.info.prefix)) {
      colnames(protein.info) = paste0(protein.info.prefix, colnames(protein.info))
    }


    ### Resolve the collisions with the columns of the results table
    # only the annotation columns are renamed, the results columns are left untouched
    colliding.columns = intersect(colnames(protein.info), colnames(data))

    if (length(colliding.columns) > 0) {
      unique.names = make.unique(c(colnames(data), colnames(protein.info)), sep = ".")
      colnames(protein.info) = unique.names[(ncol(data) + 1):length(unique.names)]

      warning("Some columns of the 'protein.info' table have the same name as columns of the results table (",
              paste(colliding.columns, collapse = ", "), ").\n",
              "The annotation columns have been renamed. Use 'protein.info.prefix' to control their naming.",
              call. = FALSE, immediate. = TRUE)
    }


    ### Bind by protein ID and not by position
    aligned.info = protein.info[match(as.character(data$prot.id), rownames(protein.info)), , drop = FALSE]
    rownames(aligned.info) = NULL

    original.rownames = rownames(data)
    data = cbind(data, aligned.info)
    rownames(data) = original.rownames

    return(data)
  } # END function
