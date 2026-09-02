#' @title get.results
#'
#' @description Simplifies the access to the differential expression results table.
#' When the object contains a protein annotation table (\code{protein.info} slot,
#' see \link{load.counts2} and \link{add.protein.info}), the annotation columns can
#' be appended to the results.
#'
#' @param DEprot.analyses.object An object of class \code{DEprot.analyses}.
#' @param contrast Number indicating the position of the contrast to use for the plotting.
#' @param protein.info.columns String vector indicating which columns of the \code{protein.info} slot should be appended to the results table. Two keywords are available: \code{"none"}, to return the results alone, and \code{"all"}, to append the whole annotation. Any other value is interpreted as the names of the columns to append. Ignored when no annotation is stored in the object. Default: \code{"none"}.
#' @param protein.info.prefix String to prepend to the names of the appended annotation columns (e.g., \code{"info."}), useful to distinguish them from the columns of the results. If \code{NULL} (default), the original names are kept and only the ones colliding with the results columns are made unique.
#'
#' @return A data frame.
#'
#' @seealso \link{add.protein.info}, \link{load.counts2}
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' # Results alone
#' get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
#'             contrast = 1)
#'
#' # Results with the whole protein annotation (when available)
#' get.results(DEprot.analyses.object = DEprot::test.toolbox$diff.exp.limma,
#'             contrast = 1,
#'             protein.info.columns = "all")
#'
#'
#' @export get.results

get.results =
  function(DEprot.analyses.object,
           contrast = 1,
           protein.info.columns = "none",
           protein.info.prefix = NULL) {

    ### check object
    if (!(methods::is(DEprot.analyses.object, "DEprot.analyses"))) {
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


    ### Define the annotation columns to append
    protein.info = .get.protein.info(DEprot.analyses.object)

    if (is.null(protein.info.columns)) {protein.info.columns = "none"}

    if (!is.character(protein.info.columns) | length(protein.info.columns) == 0) {
      stop(paste0("The 'protein.info.columns' must be a string vector indicating the annotation columns to append, ",
                  "or one of the keywords 'none' (default) and 'all'."),
           call. = FALSE)
    }

    # keyword ('none'/'all') or explicit column names?
    if (length(protein.info.columns) == 1 & tolower(protein.info.columns[1]) %in% c("none", "all")) {
      keyword = tolower(protein.info.columns[1])

      # a column of the annotation named literally 'none' or 'all' would be ambiguous
      if (!is.null(protein.info)) {
        if (protein.info.columns[1] %in% colnames(protein.info)) {
          warning(paste0("The 'protein.info' table contains a column called '", protein.info.columns[1],
                         "', which is also a keyword of 'protein.info.columns': the keyword takes precedence.\n",
                         "Rename the column to be able to select it."),
                  call. = FALSE, immediate. = TRUE)
        }
      }
    } else {
      keyword = "custom"
    }


    ### Append the protein annotation, when requested and available
    if (keyword != "none" & !is.null(protein.info)) {

      ## subset the columns, if required
      if (keyword == "custom") {
        missing.columns = setdiff(protein.info.columns, colnames(protein.info))

        if (length(missing.columns) > 0) {
          stop(paste0("The following 'protein.info.columns' are not available in the annotation table: ",
                      paste0(missing.columns, collapse = ", "), ".\n",
                      "Available columns: ", paste0(colnames(protein.info), collapse = ", "), "."),
               call. = FALSE)
        }

        protein.info = protein.info[, protein.info.columns, drop = FALSE]
      }

      ## add the prefix, if required
      if (!is.null(protein.info.prefix)) {
        colnames(protein.info) = paste0(protein.info.prefix, colnames(protein.info))
      }

      ## resolve the collisions with the columns of the results table
      # only the annotation columns are renamed, the results columns are left untouched
      colliding.columns = intersect(colnames(protein.info), colnames(data))

      if (length(colliding.columns) > 0) {
        unique.names = make.unique(c(colnames(data), colnames(protein.info)), sep = ".")
        colnames(protein.info) = unique.names[(ncol(data) + 1):length(unique.names)]

        warning(paste0("Some columns of the 'protein.info' table have the same name as columns of the results table (",
                       paste0(colliding.columns, collapse = ", "), ").\n",
                       "The annotation columns have been renamed. Use 'protein.info.prefix' to control their naming."),
                call. = FALSE, immediate. = TRUE)
      }

      ## bind by protein ID (and not by position: the results might have been re-ordered/filtered)
      aligned.info = protein.info[match(as.character(data$prot.id), rownames(protein.info)), , drop = FALSE]
      rownames(aligned.info) = NULL

      original.rownames = rownames(data)
      data = cbind(data, aligned.info)
      rownames(data) = original.rownames

    } else if (keyword == "custom" & is.null(protein.info)) {
      # explicit columns were requested, but the object carries no annotation at all
      warning(paste0("No 'protein.info' table is stored in this object: the requested columns (",
                     paste0(protein.info.columns, collapse = ", "), ") cannot be appended.\n",
                     "An annotation can be added with `add.protein.info()`."),
              call. = FALSE, immediate. = TRUE)
    }


    ## Add contrast to data attributes
    attributes(data)$contrast = names(DEprot.analyses.object@analyses.result.list)[contrast]


    ### Export table
    return(data)
  } # END of function
