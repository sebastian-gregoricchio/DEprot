#' @title export.external
#'
#' @description Converts a \code{DEprot}/\code{DEprot.analyses} object into the container used by
#'   other proteomics and Bioconductor workflows: \code{SummarizedExperiment}, \code{QFeatures},
#'   \code{MSnSet}, \code{limma::EList}, or a plain list of tables. This is the counterpart of
#'   \link{import.external}: all the count matrices available in the object become assays, the
#'   metadata table becomes the column annotation, and the protein information (together with the
#'   differential results, when present) becomes the row annotation.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param format String indicating the output container. One among: \code{"SummarizedExperiment"}
#'   (default), \code{"QFeatures"}, \code{"MSnSet"}, \code{"EList"}, \code{"list"}.
#' @param counts.type String indicating which counts must be used as primary assay. One among
#'   \code{"auto"} (default), \code{"raw"}, \code{"normalized"}, \code{"randomized"},
#'   \code{"imputed"}. When \code{"auto"} the most processed matrix available is used.
#' @param assays String or character vector indicating which count matrices should be stored as
#'   assays. Default \code{"all"} (all the matrices available in the object). Ignored by the
#'   formats holding a single matrix (\code{MSnSet}, \code{EList}).
#' @param add.protein.info Logical, whether the \code{protein.info} table of the object, when
#'   available, should be added to the row annotation. Default \code{TRUE}.
#' @param add.results Logical, whether the differential results should be added to the row
#'   annotation (only for \code{DEprot.analyses} objects). Default \code{TRUE}.
#' @param contrast.subset Numeric vector indicating the contrasts to add to the row annotation.
#'   Default \code{NULL} (all contrasts).
#' @param assay.name String indicating the name of the assay set in the \code{QFeatures} object.
#'   Default \code{"proteins"}.
#' @param keep.object Logical, whether the original \code{DEprot} object should be stored in the
#'   metadata of the exported container, allowing a lossless round-trip. Default \code{FALSE}.
#' @param install.missing String defining the behaviour when an optional package is missing:
#'   \code{"ask"} (default, prompts in interactive sessions), \code{"always"}, \code{"never"}.
#' @param verbose Logical value to indicate whether progress messages should be printed.
#'   Default \code{TRUE}.
#'
#' @details
#' The count matrices are stored as assays named \code{raw}, \code{normalized}, \code{randomized}
#' and \code{imputed}, with the one indicated by \code{counts.type} placed first so that it becomes
#' the default assay. Since a \code{SummarizedExperiment} requires all the assays to share the same
#' dimensions, the matrices are re-indexed on the proteins and samples of the primary assay: proteins
#' removed at a later step of the pipeline (imputation, filtering) are therefore dropped, and any
#' protein missing from a secondary matrix is filled with \code{NA}.
#'
#' The experimental parameters that have no equivalent in the destination class (log base,
#' normalization and imputation methods, contrasts, thresholds) are written in the metadata list of
#' the object, so that nothing is silently lost.
#'
#' Except for \code{"EList"} and \code{"list"}, which need no additional package, the destination
#' classes come from Bioconductor and are declared as optional dependencies:
#' \describe{
#'   \item{\code{SummarizedExperiment}}{\code{format = "SummarizedExperiment"} and \code{"QFeatures"}}
#'   \item{\code{QFeatures}}{\code{format = "QFeatures"}}
#'   \item{\code{MSnbase}}{\code{format = "MSnSet"}}
#' }
#' They are requested only when needed and never installed without an explicit confirmation
#' (see \code{install.missing}).
#'
#' @return An object of the class indicated by \code{format}.
#'
#' @import dplyr
#' @import methods
#' @importFrom utils packageVersion
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' # SummarizedExperiment with all the count matrices available
#' se <- export.external(DEprot.object = DEprot::test.toolbox$dpo.imp,
#'                       format = "SummarizedExperiment")
#'
#' # differential results stored in the rowData
#' se <- export.external(DEprot.object = DEprot::test.toolbox$diff.exp.limma,
#'                       format = "SummarizedExperiment",
#'                       add.results = TRUE)
#'
#' # QFeatures, keeping the original object for a round-trip
#' qf <- export.external(DEprot.object = dpo,
#'                       format = "QFeatures",
#'                       keep.object = TRUE)
#'
#' # single-matrix containers
#' ms <- export.external(DEprot.object = dpo, format = "MSnSet")
#' el <- export.external(DEprot.object = dpo, format = "EList")
#' }
#'
#' @seealso \link{import.external}, \link{export.analyses}
#'
#' @export export.external

export.external =
  function(DEprot.object,
           format = "SummarizedExperiment",
           counts.type = "auto",
           assays = "all",
           add.protein.info = TRUE,
           add.results = TRUE,
           contrast.subset = NULL,
           assay.name = "proteins",
           keep.object = FALSE,
           install.missing = "ask",
           verbose = TRUE) {

    ### ------------------------------------------------------------------- checks
    if (missing(DEprot.object)) {
      stop("Provide 'DEprot.object': an object of class 'DEprot' or 'DEprot.analyses'.", call. = FALSE)
    }

    if (!any(c("DEprot", "DEprot.analyses") %in% class(DEprot.object))) {
      stop("The input must be an object of class 'DEprot' or 'DEprot.analyses'.", call. = FALSE)
    }

    format = match.arg(tolower(as.character(format)[1]),
                       choices = c("summarizedexperiment", "qfeatures", "msnset", "elist", "list"))

    counts.type = match.arg(tolower(as.character(counts.type)[1]),
                            choices = c("auto", "raw", "normalized", "randomized", "imputed"))

    install.missing = match.arg(tolower(as.character(install.missing)[1]),
                                choices = c("ask", "always", "never"))


    ### ------------------------------------------------------------------- assays
    collected = .collect.assays(DEprot.object = DEprot.object,
                                counts.type = counts.type,
                                assays = assays,
                                verbose = verbose)

    features = rownames(collected$assays[[1]])
    samples  = colnames(collected$assays[[1]])

    if (isTRUE(verbose)) {
      message("Exporting ", base::format(length(features), big.mark = ","), " proteins x ",
              length(samples), " samples (primary assay: '", collected$reference, "').")
    }


    ### ------------------------------------------------------------------- annotations
    row.data = .build.row.annotation(DEprot.object = DEprot.object,
                                     features = features,
                                     add.protein.info = add.protein.info,
                                     add.results = add.results,
                                     contrast.subset = contrast.subset,
                                     verbose = verbose)

    col.data = .build.col.annotation(DEprot.object = DEprot.object,
                                     samples = samples)

    obj.metadata = .build.export.metadata(DEprot.object = DEprot.object,
                                          keep.object = keep.object)


    ### ------------------------------------------------------------------- build the container
    exported =
      switch(format,
             "summarizedexperiment" = .build.summarized.experiment(assay.list = collected$assays,
                                                                   row.data = row.data,
                                                                   col.data = col.data,
                                                                   obj.metadata = obj.metadata,
                                                                   install.missing = install.missing),

             "qfeatures" = .build.qfeatures(assay.list = collected$assays,
                                            row.data = row.data,
                                            col.data = col.data,
                                            obj.metadata = obj.metadata,
                                            assay.name = assay.name,
                                            install.missing = install.missing),

             "msnset" = .build.msnset(assay.list = collected$assays,
                                      row.data = row.data,
                                      col.data = col.data,
                                      obj.metadata = obj.metadata,
                                      install.missing = install.missing,
                                      verbose = verbose),

             "elist" = .build.elist(assay.list = collected$assays,
                                    row.data = row.data,
                                    col.data = col.data,
                                    verbose = verbose),

             "list" = list(counts = collected$assays,
                           protein.info = row.data,
                           metadata = col.data,
                           parameters = obj.metadata))

    return(exported)
  } # END function




#' @title as.SummarizedExperiment
#'
#' @description Shortcut for \code{export.external(..., format = "SummarizedExperiment")}.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param ... Any other parameter accepted by \link{export.external}.
#'
#' @return A \code{SummarizedExperiment} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' se <- as.SummarizedExperiment(DEprot::test.toolbox$diff.exp.limma)
#' }
#'
#' @seealso \link{export.external}
#'
#' @export as.SummarizedExperiment

as.SummarizedExperiment =
  function(DEprot.object, ...) {
    export.external(DEprot.object = DEprot.object, format = "SummarizedExperiment", ...)
  } # END function




#' @title as.QFeatures
#'
#' @description Shortcut for \code{export.external(..., format = "QFeatures")}.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param ... Any other parameter accepted by \link{export.external}.
#'
#' @return A \code{QFeatures} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' qf <- as.QFeatures(dpo, assay.name = "proteins")
#' }
#'
#' @seealso \link{export.external}
#'
#' @export as.QFeatures

as.QFeatures =
  function(DEprot.object, ...) {
    export.external(DEprot.object = DEprot.object, format = "QFeatures", ...)
  } # END function




#' @title as.MSnSet
#'
#' @description Shortcut for \code{export.external(..., format = "MSnSet")}.
#'
#' @param DEprot.object An object of class \code{DEprot} or \code{DEprot.analyses}.
#' @param ... Any other parameter accepted by \link{export.external}.
#'
#' @return An \code{MSnSet} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @examples
#' \dontrun{
#' ms <- as.MSnSet(dpo, counts.type = "imputed")
#' }
#'
#' @seealso \link{export.external}
#'
#' @export as.MSnSet

as.MSnSet =
  function(DEprot.object, ...) {
    export.external(DEprot.object = DEprot.object, format = "MSnSet", ...)
  } # END function
